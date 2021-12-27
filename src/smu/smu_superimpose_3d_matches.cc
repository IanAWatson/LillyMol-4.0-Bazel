#include <iostream>
#include <memory>

#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/u3b.h"

// Processing the output of topology_from_geometry, we have
// a list of conformer_id's and bond_topology_id's where a
// different BondTopology matches the conformer.
// For each conformer_id, fetch all conformers of the corresponding
// BondTopology and try to superimpose them.

namespace superimpose_3d_matches {
using std::cerr;

int verbose = 0;

char output_separator = ' ';

void
Usage(int rc) {
  cerr << " -R <fname>      output from reconcile_dupes\n";
  cerr << " -v              verbose output\n";

  exit(rc);
}


// The task is to find the best match for a given ConformerID.
// This class 
class PossibleMatches {
  private:
    Molecule * _target;

    resizable_array_p<Molecule> _mols;
  public:
    PossibleMatches();
    ~PossibleMatches();

    void SetTarget(Molecule * m);

    bool HasTarget() const {
      return _target != nullptr;
    }

    int Add(Molecule * m) {
      return _mols.add(m);
    }

    int NumberCandidates() const {
      return _mols.number_elements();
    }

    int ProcessMatches(std::ostream& output);

    int Write(std::ostream& output);
};

PossibleMatches::PossibleMatches() {
  _target = nullptr;
}

PossibleMatches::~PossibleMatches() {
  if (_target != nullptr) {
    delete _target;
  }
}

void
PossibleMatches::SetTarget(Molecule* m) {
  _target = m;
}

int
PossibleMatches::Write(std::ostream & output) {
  if (_target == nullptr || _mols.number_elements() == 0) {
    return 0;
  }
  output << _target->smiles() << output_separator << _target->name() << output_separator << _mols.number_elements() << output_separator << " T\n";
  IWString tname(_target->name());
  tname.truncate_at_first(' ');
  for (Molecule* m : _mols) {
    output << m->smiles() << output_separator << m->name() << output_separator << tname << '\n';
  }

  return 1;
}

// Take the output from u3b and orient `m`.
void
Orient(Molecule& m,
       const std::unique_ptr<double[]>& t,
       const std::unique_ptr<double[]>& u) {

  const double rotmat11 = u[0];
  const double rotmat12 = u[1];
  const double rotmat13 = u[2];
  const double rotmat21 = u[3];
  const double rotmat22 = u[4];
  const double rotmat23 = u[5];
  const double rotmat31 = u[6];
  const double rotmat32 = u[7];
  const double rotmat33 = u[8];

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++)
  {
    const Atom * a = m.atomi(i);

    const double x0 = a->x() - t[0];
    const double y0 = a->y() - t[1];
    const double z0 = a->z() - t[2];

    double xx = rotmat11 * x0 + rotmat12 * y0 + rotmat13 * z0;
    double yy = rotmat21 * x0 + rotmat22 * y0 + rotmat23 * z0;
    double zz = rotmat31 * x0 + rotmat32 * y0 + rotmat33 * z0;

    m.setxyz(i, static_cast<coord_t>(xx), static_cast<coord_t>(yy), static_cast<coord_t>(zz));
  }
}

int
PossibleMatches::ProcessMatches(std::ostream& output) {
  const int matoms = _target->natoms();

  Molecule mcopy(*_target);
  mcopy.set_all_bonds_to_type(SINGLE_BOND);
  for (int i = 0; i < matoms; ++i) {
    mcopy.set_formal_charge(i, 0);
  }

  Molecule_to_Query_Specifications mqs;
  mqs.set_just_atomic_number_and_connectivity(1);
  mqs.set_make_embedding(1);

  Substructure_Query target_query;
  if (! target_query.create_from_molecule(mcopy, mqs)) {
    cerr << "PossibleMatches::ProcessMatches:cannot convert target to query\n";
    return 0;
  }

#ifdef NOT_NEEDEDJ
  queries.resize(_mols.number_elements());
  for (Molecule * m : _mols) {
    cerr << _target->unique_smiles() << " cmp " << m->unique_smiles() << '\n';
    std::unique_ptr<Substructure_Query> qry = std::make_unique<Substructure_Query>();
    if (! qry->create_from_molecule(*m, mqs)) {
      cerr << "PossibleMatches::ProcessMatches:cannot create query from molecule\n";
      return 0;
    }
    queries.add(qry.release());
  }
#endif

  std::unique_ptr<double[]> weight(new_double(matoms, 1.0));
  std::unique_ptr<double[]> x(new_double(matoms * 3));
  std::unique_ptr<double[]> y(new_double(matoms * 3));
  std::unique_ptr<double[]> u(new_double(matoms * matoms));
  std::unique_ptr<double[]> t(new_double(matoms));
  double rms = 0.0;
  int mode = 1;
  int ier = 0;

  for (int i = 0; i < matoms; ++i) {
    const Atom& a = _target->atom(i);
    x[3 * i] = a.x();
    x[3 * i + 1] = a.y();
    x[3 * i + 2] = a.z();
  }

  cerr << _target->name() << " begin search across " << _mols.number_elements() << " matches\n";

  // Ambiguous as to whether we should get the best match for each possible match
  // or globally across all.
  Molecule best_match;
  float min_rms = std::numeric_limits<float>::max();
  for (Molecule * m : _mols) {
    Substructure_Results sresults;
    const int nhits = target_query.substructure_search(*m, sresults);
    if (nhits == 0) {
      cerr << "No hits for " << _target->unique_smiles() << " in " << m->unique_smiles() << '\n';
      continue;
    }
    cerr << " match has " << nhits << " substructure matches\n";
    for (int i = 0; i < nhits; ++i) {
      const Set_of_Atoms & embedding = *sresults.embedding(i);
      for (int j = 0; j < matoms; ++j) {
        int k = embedding[j];
        if (_target->atomic_number(j) != m->atomic_number(k)) {
          cerr << "ATomic number mismatch\n";
        }
        const Atom & a = *m->atomi(embedding[j]);
        y[3 * j] = a.x();
        y[3 * j + 1] = a.y();
        y[3 * j + 2] = a.z();
      }
      u3b_(weight.get(), x.get(), y.get(), &matoms, &mode, &rms, u.get(), t.get(), &ier);

      if (0 == ier)
        ;
      else if (-1 == ier)
        cerr << "PossibleMatches::ProcessMatches:superposition not unique, but optimal\n";
      else {
        cerr << "PossibleMatches::ProcessMatches:u3b failed\n";
        continue;
      }

      if (rms < min_rms) {
        min_rms = rms;
        best_match = *m;
        Orient(best_match, t, u);
      }
    }
  }
  cerr << "Min rms " << min_rms << '\n';
  output << _target->smiles() << output_separator << _target->name() << '\n';
  output << best_match.smiles() << output_separator << best_match.name() << ' ' << min_rms << '\n';

  return 1;
}

IW_STL_Hash_Map<IWString, PossibleMatches*>
ReadPreFoundMatches(data_source_and_type<Molecule> & input) {
  IW_STL_Hash_Map<IWString, PossibleMatches*> result;

  Molecule * m;
  IWString current_target;
  while ((m = input.next_molecule()) != NULL) {
    if (m->name().ends_with(" T")) {
      IWString tmp(m->name());
      tmp.truncate_at_first(' ');
      cerr << "Is target '" << tmp << "'\n";

      PossibleMatches* p = new PossibleMatches;
      p->SetTarget(m);
      std::pair<IWString, PossibleMatches*> mypair(tmp, p);
      result.insert(mypair);
      current_target = tmp;
    } else {
      result[current_target]->Add(m);
      cerr << "Added as possible match\n";
    }
  }

  return result;
}


IW_STL_Hash_Map<IWString, PossibleMatches*>
ReadPreFoundMatches(const char * fname) {
  data_source_and_type<Molecule> input(FILE_TYPE_SMI, fname);
  if (! input.good()) {
    cerr << "ReadPreFoundMatches:cannot open " << fname << '\n';
    return IW_STL_Hash_Map<IWString, PossibleMatches*>();
  }

  return ReadPreFoundMatches(input);
}

int
ProcessPreviouslyExtracted(const char * fname, std::ostream& output) {
  IW_STL_Hash_Map<IWString, PossibleMatches*> possible_matches = ReadPreFoundMatches(fname);
  if (possible_matches.empty()) {
    cerr << "ProcessPreviouslyExtracted:no data\n";
    return 0;
  }

  int has_target = 0;
  extending_resizable_array<int> num_candidates;
  for (const auto& [confid, matches] : possible_matches) {
    if (! matches->HasTarget()) {
      cerr << "NO targat\n";
      continue;
    }
    has_target++;

    num_candidates[matches->NumberCandidates()]++;
    if (matches->NumberCandidates() == 0) {
      cerr << "no matches\n";
      continue;
    }
    matches->ProcessMatches(output);
  }

  return 1;
}

// return 1 if `m` is needed by the calculation.
// It is needed if `cid_to_process` contains m->name.
// It is needed if `btid_to_cid`  contains m->name
int
MaybeAddToHash(Molecule* m,
               const IW_STL_Hash_Map_String& btid_to_cid,
               IW_STL_Hash_Map<IWString, PossibleMatches*> possible_matches) {
  IWString name = m->name();
  name.truncate_at_first(' ');
  // First check to see if this is needed as a target.
  auto iter1 = possible_matches.find(name);
  if (iter1 != possible_matches.end()) {
    iter1->second->SetTarget(m);
    return 1;
  }

  // Not a target, see if bond topology id needed as a possible match.
  name.chop(3);
  auto iter2 = btid_to_cid.find(name);
  if (iter2 == btid_to_cid.end()) {
    return 0;
  }

  const IWString cid = iter2->second;

  // Great, needed as a match

  return possible_matches[cid]->Add(m);
}

int
ReadMolecules(data_source_and_type<Molecule>& input,
              const IW_STL_Hash_Map_String& btid_to_cid,
              IW_STL_Hash_Map<IWString, PossibleMatches*> possible_matches) {
  Molecule * m;
  while ((m = input.next_molecule()) != NULL) {
    if (MaybeAddToHash(m, btid_to_cid, possible_matches)) {
    } else {
      delete m;
    }
  }

  return 1;
}

int
ReadMolecules(const char * fname,
              const IW_STL_Hash_Map_String& btid_to_cid,
              IW_STL_Hash_Map<IWString, PossibleMatches*> possible_matches) {
  data_source_and_type<Molecule> input(FILE_TYPE_SMI, fname);
  if (! input.good()) {
    cerr << "Cannot open " << fname << '\n';
    return 0;
  }

  return ReadMolecules(input, btid_to_cid, possible_matches);
}

IW_STL_Hash_Map_String
ReadIdsToProcess(iwstring_data_source& input) {
  IW_STL_Hash_Map_String result;

  const_IWSubstring buffer;
  input.next_record(buffer);
  while (input.next_record(buffer)) {
    IWString cid, btid;
    buffer.split(cid, ' ', btid);
    result.emplace(std::make_pair(btid, cid));
  }

  return result;
}

IW_STL_Hash_Map_String
ReadIdsToProcess(const char * fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "ReadIdsToProcess:cannot open " << fname << '\n';
    return IW_STL_Hash_Map_String();
  }

  return ReadIdsToProcess(input);
}

int
Superimpose3dMatches(int argc, char ** argv) {
  Command_Line cl(argc, argv, "vi:A:g:R:P:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  verbose = cl.option_count('v');

  set_append_coordinates_after_each_atom(1);

  if (cl.option_present('P')) {
    const char * fname = cl.option_value('P');
    cerr << "Extract prev from '" << fname << "'\n";
    ProcessPreviouslyExtracted(fname, std::cout);
    return 0;
  }

  // The set of conformer id's and bond_topology ids to be processed.
  // For each conformer id, we need to find the corresponding Molecule.
  // For each conformer derived from the bond topology id, we need to
  // accumulate those molecules.
  IW_STL_Hash_Map_String btid_to_cid;

  if (! cl.option_present('R')) {
    cerr << "MUst specify the output from reconcile_dupes via the -R option\n";
    Usage(1);
  } else {
    IWString fname = cl.string_value('R');
    btid_to_cid = ReadIdsToProcess(fname);
    if (btid_to_cid.empty()) {
      cerr << "Cannot read reconcile_dupes ids to process from " << fname << "\n";
      return 1;
    }

    if (verbose) {
      cerr << "Read " << btid_to_cid.size() << " ids to process from " << fname << '\n';
    }
  }

  // For each cid, we will need a set of possible matches.
  IW_STL_Hash_Map<IWString, PossibleMatches*> possible_matches;
  for (auto& [btid, cid] : btid_to_cid) {
    PossibleMatches* p = new PossibleMatches;
    std::pair<IWString, PossibleMatches*> mypair(cid, p);
    possible_matches.insert(mypair);
  }

  if (verbose) {
    cerr << "Looking for structures to match " << possible_matches.size() << " conformers\n";
  }

  if (cl.empty()) {
    cerr << "Must specify input file(s)\n";
    Usage(1);
  }

  for (const char * fname : cl) {
    if (! ReadMolecules(fname, btid_to_cid, possible_matches)) {
      cerr << "Error reading " << fname << '\n';
      return 1;
    }
  }

  if (verbose) {
    for (const auto& [confid, matches]: possible_matches) {
      cerr << "for " << confid << " has target " << matches->HasTarget() << " matches " << matches->NumberCandidates() << '\n';
    }
  }

  int has_target = 0;
  extending_resizable_array<int> num_candidates;
  for (const auto& [confid, matches] : possible_matches) {
    if (! matches->HasTarget()) {
      continue;
    }
    has_target++;

    num_candidates[matches->NumberCandidates()]++;
    if (matches->NumberCandidates() == 0) {
      continue;
    }
//  matches->ProcessMatches(std::cout);
  }

  for (const auto& [confid, matches] : possible_matches) {
    matches->Write(std::cout);
  }

  return 0;
}

}  // namespace superimpose_3d_matches

int
main(int argc, char ** argv) {
  return superimpose_3d_matches::Superimpose3dMatches(argc, argv);
}
