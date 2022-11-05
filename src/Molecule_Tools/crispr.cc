// A generalization of ring replacement. Given a set of replacement
// fragments, look for each of those fragments in a target molecule and
// try to replace what is there.
// The fragments are converted to queries via mol2qry, so necessarily
// the atom types will be discarded.
// I called this molecular_crispr, but I really don't understand CRISPR
// so I am not sure the name is appropriate.

#include <iostream>
#include <memory>

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

namespace molecular_crispr {

using std::cerr;

void
Usage(int rc) {
  cerr << "Create new molecules by merging in replacement fragments\n";
  cerr << " -F <fname>       file containing molecular fragments\n";
  cerr << " -Y <query>       queries that the product molecules must     contain\n";
  cerr << " -N <query>       queries that the product molecules must not contain\n";
  cerr << " -p               write the parent molecule\n";

  cerr << " -c                remove chirality\n";
  cerr << " -l                strip to largest fragment\n";
  cerr << " -v                verbose output\n";

  ::exit(rc);
}

// Fragments are read and converted to query form. We retain
// the original molecule so it can be inserted into the parent.
class MoleculeAndQuery {
  private:
    Molecule _mol;
    Substructure_Query _query;

    int _molecules_generated;

  public:
    MoleculeAndQuery();

    int Build(const Molecule& m, Molecule_to_Query_Specifications& mqs);

    Substructure_Query& query() {
      return _query;
    }

    int molecules_generated() const {
      return _molecules_generated;
    }

    int substructure_search(Molecule_to_Match& target, Substructure_Results& sresults) {
      return _query.substructure_search(target, sresults);
    }

    int Process(Molecule& m, const Set_of_Atoms& embedding);

    void SetName(Molecule& m) const;
};

MoleculeAndQuery::MoleculeAndQuery() {
  _molecules_generated = 0;
}

int
MoleculeAndQuery::Build(const Molecule& m, Molecule_to_Query_Specifications& mqs) {
  if (! _query.create_from_molecule(m, mqs)) {
    cerr << "MoleculeAndQuery::Build:Cannot convert to query '" << m.name() << '\n';
    return 0;
  }

  _mol = m;

  return 1;
}

void
UnsetImplicitHydrogenInformation(Molecule& m) {
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    m.unset_all_implicit_hydrogen_information(i);
  }
}

//#define DEBUG_MQ_PROCESS

// Code copied from Duplicate/ring_replacement.cc
// Note that this may invert chiral centres.
// If the chiral centre has an implicit Hydrogen,
// adding a bond will remove the lone pair and make
// the new atom part of the chiral centre. But then
// the old atom is removed, and that position
// then becomes an implicit hydrogen.
int
MoleculeAndQuery::Process(Molecule& m,
                          const Set_of_Atoms& embedding) {
#ifdef DEBUG_MQ_PROCESS
  cerr << "start smiles " << m.unique_smiles() << '\n';
  cerr << "Embedding " << embedding << '\n';
#endif
  const int initial_natoms = m.natoms();
  const int initial_aromatic_count = m.aromatic_atom_count();
  const IWString initial_smiles = m.unique_smiles();

  m.add_molecule(&_mol);  
#ifdef DEBUG_MQ_PROCESS
  cerr << "After adding fragment\n";
  write_isotopically_labelled_smiles(m, false, cerr);
#endif

  std::unique_ptr<int[]> to_remove(new_int(m.natoms()));
  embedding.set_vector(to_remove.get(), 1);

  // Add pairs of atoms to this vector.
  Set_of_Atoms replacements;

  const int esize = embedding.number_elements();
  for (int i = 0; i < esize; ++i) {
    const atom_number_t a = embedding[i];
    const Atom& atom = m.atom(a);
    for (const Bond* b : atom) {
      atom_number_t o = b->other(a);
      if (to_remove[o]) {
        continue;
      }
      // cerr << "From atom " << a << " make bond " << (initial_natoms + i) << " to " << o << '\n';

      atom_number_t join_to = initial_natoms + i;
      if (m.implicit_hydrogens(join_to) == 0) {
        return 0;
      }
      replacements << a;
      replacements << join_to;
#ifdef DEBUG_MQ_PROCESS
      cerr << "   adding bond btw " << o << " and " << join_to << '\n';
#endif
      continue;
    }
  }

#ifdef DEBUG_MQ_PROCESS
  cerr << "Replacements " << replacements << '\n';
#endif

  for (int i = 0; i < replacements.number_elements(); i += 2) {
    m.stereo_preserving_substitute(replacements[i], replacements[i + 1]);
  }

#ifdef DEBUG_MQ_PROCESS
  cerr << "After bonds added " << m.unique_smiles() << '\n';
#endif

  m.remove_atoms(to_remove.get());

  m.reduce_to_largest_fragment();

  UnsetImplicitHydrogenInformation(m);

  if (int final_aromatic_atoms = m.aromatic_atom_count(); final_aromatic_atoms < initial_aromatic_count) {
    cerr << initial_aromatic_count << " vs " << final_aromatic_atoms << " aromatic atoms\n";
    cerr << "Warning loss of aromaticty initial " << initial_smiles << " now "  <<  m.unique_smiles() << '\n';
  }

  if (! m.valence_ok()) {
    cerr << "Invalid valence '" << m.smiles() << ' ' << m.name() << ' ' << _mol.name() << '\n';
    for (int i = 0; i < m.natoms(); ++i) {
      if (! m.valence_ok(i)) {
        cerr << " atom " << i << ' ' << m.smarts_equivalent_for_atom(i) << '\n';
      }
    }
    return 0;
  }

  ++_molecules_generated;

  return 1;
}

void
MoleculeAndQuery::SetName(Molecule& m) const {
  IWString new_name = m.name();
  new_name << " %% " << _mol.name();
  m.set_name(new_name);
}

class CoreReplacement {
  private:
    int _verbose;

    int _molecules_read;

    int _reduce_to_largest_fragment;

    int _remove_chirality;

    // Queries formed by converting the -F file to queries.
    resizable_array_p<MoleculeAndQuery> _fragment;

    int _ignore_fragments_not_meeting_requirements;

    int _write_parent_molecule;

    Molecule_to_Query_Specifications _mol2qry;

    IW_STL_Hash_Set _seen;

    extending_resizable_array<int> _matches_per_molecule;

    int _molecules_created;
    int _molecules_written;

    FileType _input_type;

    Chemical_Standardisation _chemical_standardisation;

    resizable_array_p<Substructure_Query> _results_must_have;
    resizable_array_p<Substructure_Query> _results_must_not_have;

    int _duplicates_discarded;
    int _failed_constraints;

  // Private functions.

    int ReadFragments(IWString& fname);
    int ReadFragments(data_source_and_type<Molecule>& input);
    int AddFragment(Molecule& m);

    int ReadMol2Qry(IWString& fname);

    int MatchesMustHave(Molecule_to_Match& target);
    int MatchesMustNotHave(Molecule_to_Match& target);

    int PassesConstraints(Molecule& m);

  public:
    CoreReplacement();

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);

    // For now, only smiles output is supported. But smiles can hold
    // 3d information.
    int Process(Molecule& m, IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;

    FileType input_type() const {
      return _input_type;
    }
};

CoreReplacement::CoreReplacement() {
  _verbose = 0;
  _molecules_read = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _ignore_fragments_not_meeting_requirements = 0;
  _write_parent_molecule = 0;
  _molecules_created = 0;
  _molecules_written = 0;
  _failed_constraints = 0;
  _duplicates_discarded = 0;
  _input_type = FILE_TYPE_INVALID;
}

int
CoreReplacement::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('c')) {
    _remove_chirality = 1;
    if (_verbose) {
      cerr << "Will remove chirality from input molecules\n";
    }
  }

  if (cl.option_present('g')) {
    if (! _chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      cerr << "Cannot initialise chemical standardisation\n";
      return 0;
    }
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce molecules to largest fragment\n";
    }
  }

  if (cl.option_present('P')) {
    IWString fname = cl.string_value('P');
    if (! ReadMol2Qry(fname)) {
      cerr << "Cannot initialise mol2qry proto '" << fname << "'\n";
      return 0;
    }
  }

  if (! cl.option_present('F')) {
    cerr << "Must specify a file of replacement fragments via the -F option\n";
    return 0;
  }

  if (cl.option_present('F')) {
    IWString fname;
    for (int i = 0; cl.value('F', fname, i); ++i) {
      if (! ReadFragments(fname)) {
        cerr << "Cannot read replacement fragments from '" << fname << "'\n";
        return 0;
      }
    }
  }

  //for (MoleculeAndQuery* q : _fragment) {
  //  q->query().set_find_unique_embeddings_only(1);
  //}

  for (int i = 0; i < _fragment.number_elements(); ++i) {
    IWString fname;
    fname << "/tmp/q" << i << ".qry";
    _fragment[i]->query().write_msi(fname);
  }

  if (_verbose) {
    cerr << "CoreReplacement::initialise:read " << _fragment.size() << " replacement fragments\n";
  }

  if (cl.option_present('p')) {
    _write_parent_molecule = 1;
    if (_verbose) {
      cerr << "Will write the parent molecule\n";
    }
  }

  if (cl.option_present('Y')) {
    if (! process_queries(cl, _results_must_have, _verbose, 'Y')) {
      cerr << "Cannot process results must have queries (-Y)\n";
      return 0;
    }
  }

  if (cl.option_present('N')) {
    if (! process_queries(cl, _results_must_not_have, _verbose, 'N')) {
      cerr << "Cannot process results must not have queries (-N)\n";
      return 0;
    }
  }

  if (1 == cl.number_elements() && 0 == strcmp("-", cl[0])) { // reading a pipe, assume smiles
    _input_type = FILE_TYPE_SMI;
  } else if (!all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot discern all file types, use the -i option\n";
    return 0;
  } else if (!process_input_type(cl, _input_type)) {
    return 0;
  }

  return 1;
}

int
CoreReplacement::ReadMol2Qry(IWString& fname) {
  std::optional<molecule_to_query::MoleculeToQuery> proto =
       iwmisc::ReadTextProto<molecule_to_query::MoleculeToQuery>(fname);
  if (!proto) {
    cerr << "CoreReplacement::ReadMol2Qry:cannot read '" << fname << "'\n";
    return 0;
  }
  cerr << "Just read mol2qry proto " << fname << "'\n";

  return _mol2qry.Build(*proto);
}

int
CoreReplacement::ReadFragments(IWString& fname) {
  FileType input_type = _input_type;
  if (_input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    if (input_type == FILE_TYPE_INVALID) {
      cerr << "CoreReplacement::ReadFragments:cannot discern fragment file type '" << fname << "'\n";
      return 0;
    }
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok()) {
    cerr << "CoreReplacement::ReadFragments:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadFragments(input);
}

int
CoreReplacement::ReadFragments(data_source_and_type<Molecule>& input) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (AddFragment(*m)) {
      continue;
    }
    if (_ignore_fragments_not_meeting_requirements) {
      continue;
    }

    cerr << "CoreReplacement::ReadFragments:cannot process " << m->name() << '\n';
    return 0;
  }

  return 1;
}

int
CoreReplacement::AddFragment(Molecule& m) {
  std::unique_ptr<MoleculeAndQuery> q = std::make_unique<MoleculeAndQuery>();
  if (! q->Build(m, _mol2qry)) {
    cerr << "CoreReplacement::AddFragment:cannot create query from " << m.name() << '\n';
    return 0;
  }

  _fragment << q.release();

  return 1;
}

int
CoreReplacement::Process(Molecule& m, IWString_and_File_Descriptor& output) {
  ++_molecules_read;
  // Avoid molecules recreating themselves.
  _seen.insert(m.unique_smiles());
  m.invalidate_smiles();
  int parent_molecule_written = 0;

  Molecule_to_Match target(&m);
  Substructure_Results sresults;
  int matches = 0;
  for (MoleculeAndQuery* q : _fragment) {
    // cerr << " Trying search " << q->query().comment() << '\n';
    if (! q->substructure_search(target, sresults)) {
      continue;
    }
    ++matches;
    // cerr << "Query " << q->query().comment() << " got " << sresults.number_embeddings() << " embeddings\n";
    for (const Set_of_Atoms* e : sresults.embeddings()) {
      Molecule mcopy(m);
      // cerr << "Smiles " << m.unique_smiles() << " and " << mcopy.unique_smiles() << '\n';
      ++_molecules_created;
      if (! q->Process(mcopy, *e)) {
        continue;
      }
      // cerr << "After process " << mcopy.unique_smiles() << '\n';
      if (! PassesConstraints(mcopy)) {
        ++_failed_constraints;
        continue;
      }

      if (_write_parent_molecule && ! parent_molecule_written) {
        output << m.smiles() << ' ' << m.name() << '\n';
        parent_molecule_written = 1;
      }

      ++_molecules_written;
      q->SetName(mcopy);
      mcopy.invalidate_smiles();
      output << mcopy.smiles() << ' ' << mcopy.name() << '\n';
      output.write_if_buffer_holds_more_than(4096);
    }
  }

  ++_matches_per_molecule[matches];

  return 1;
}

int
CoreReplacement::PassesConstraints(Molecule& m) {
  if (_seen.contains(m.unique_smiles())) {
    ++_duplicates_discarded;
    return 0;
  }

  if (_results_must_have.empty() && _results_must_not_have.empty()) {
    return 1;
  }

  Molecule_to_Match target(&m);
  if (! MatchesMustHave(target)) {
    return 0;
  }

  if (MatchesMustNotHave(target)) {
    return 0;
  }

  return 1;
}

int
MatchesAny(Molecule_to_Match& target,
           const resizable_array_p<Substructure_Query>& queries) {
  Substructure_Results sresults;
  for (Substructure_Query* q : queries) {
    if (q->substructure_search(target, sresults)) {
      return 1;
    }
  }

  return 0;
}

int
CoreReplacement::MatchesMustHave(Molecule_to_Match& target) {
  if (_results_must_have.empty()) {
    return 1;
  }

  return MatchesAny(target, _results_must_have);
}

int
CoreReplacement::MatchesMustNotHave(Molecule_to_Match& target) {
  if (_results_must_not_have.empty()) {
    return 0;
  }
  return MatchesAny(target, _results_must_not_have);
}

int
CoreReplacement::Preprocess(Molecule& m) {
  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
  }

  if (m.natoms() == 0) {
    return 0;
  }

  return 1;
}

int
CoreReplacement::Report(std::ostream& output) const {
  output << "CoreReplacement:read " << _molecules_read << " molecules\n";
  if (_molecules_read == 0) {
    return 1;
  }

  output << "Created " << _molecules_created << " molecules\n";
  output << _duplicates_discarded << " duplicate molecules discarded\n";
  output << _failed_constraints << " molecules failed structural and uniqueness constraints\n";
  output << _molecules_written << " molecules written\n";

  return 1;
}

int
ReplaceCore(CoreReplacement& core_replacement,
            Molecule& m,
            IWString_and_File_Descriptor& output) {
  return core_replacement.Process(m, output);
}

int
ReplaceCore(CoreReplacement& core_replacement,
            data_source_and_type<Molecule>& input,
            IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! core_replacement.Preprocess(*m)) {
      return 0;
    }

    if (! ReplaceCore(core_replacement, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
ReplaceCore(CoreReplacement& core_replacement,
            const char * fname,
            IWString_and_File_Descriptor& output) {
  FileType input_type = core_replacement.input_type();
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return ReplaceCore(core_replacement, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:i:g:lcF:Y:N:pP:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  set_copy_name_in_molecule_copy_constructor(1);

  const int verbose = cl.option_count('v');

  if (! process_standard_aromaticity_options(cl, verbose, 'A')) {
    cerr << "Cannot process aromaticity options\n";
    return 1;
  }

  if (! process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process standard elements options (-E)\n";
    return 1;
  }

  CoreReplacement core_replacement;
  if (! core_replacement.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);
  for (const char* fname : cl) {
    if (! ReplaceCore(core_replacement, fname, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    core_replacement.Report(cerr);
  }

  return 0;
}

}  // namespace molecular_crispr

int
main(int argc, char** argv) {
  int rc = molecular_crispr::Main(argc, argv);

  return rc;
}
