// Scan a set of molecules and extract the bond angles specified by a query

#include <stdlib.h>
#include <cstdint>
#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwbits/fixed_bit_vector.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

namespace get_bond_angles {

using std::cerr;

void
Usage(int rc) {
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << "\n";
  cerr << "Accumulates bond angles specified via a set of queries\n";
  cerr << " -s <smarts>      specify queries via smarts\n";
  cerr << " -q <query >      specify queries via query files and directives\n";
  cerr << " -x <2,3,4>       specify number of matched atoms needed, 2->bond length, 3->bond_angle, 4->torsion\n";
  cerr << " -A ...           aromaticity options\n";
  cerr << " -E ...           element options\n";
  cerr << " -g ...           chemical standardisation options\n";
  cerr << " -v               verbose output\n";

  ::exit(rc);
}

constexpr int k180 = 180;

struct AngleHistogram {
  private:
    uint32_t count[k180];

  public:
    AngleHistogram();

    void Extra(float angle);

    int MakeHistogram(IWString& fname) const;

    int Report(std::ostream& output) const;
};

AngleHistogram::AngleHistogram() {
  std::fill_n(count, k180, 0);
}

void
AngleHistogram::Extra(float angle) {
  if (angle < 0.0) {
    angle = 180.0 - angle;
  }

  int ndx = static_cast<int>(angle + 0.49999);

  if (ndx >= k180) {
    ndx = ndx % k180;
  }

  ++count[ndx];
}

int
AngleHistogram::Report(std::ostream& output) const {
  Accumulator_Int<uint64_t> acc;
  for (int i = 0; i < k180; ++i) {
    if (count[i] > 0) {
      acc.extra(i, count[i]);
    }
  }

  if (acc.n() == 0) {
    output << "No query matches\n";
    return 1;
  }

  output << acc.n() << " values btw " << acc.minval() << " and " << acc.maxval() << " mean " << acc.average() << " std " << sqrt(acc.variance()) << '\n';

  return 1;
}

int
AngleHistogram::MakeHistogram(IWString& fname) const {
  IWString_and_File_Descriptor output;
  if (! output.open(fname)) {
    cerr << "AngleHistogram::MakeHistogram:cannot open '" << fname << "'\n";
    return 0;
  }

  output << "angle,count\n";
  for (int i = 0; i < k180; ++i) {
    output << i << ',' << count[i] << '\n';
  }

  return 1;
}

class BondAngles {
  private:
    int _verbose;

    int _reduce_to_largest_fragment;

    Chemical_Standardisation _chemical_standardisation;

    int _molecules_read;

    resizable_array_p<Substructure_Query> _queries;

    // Even though this tool was written for bond angles, it quickly became
    // apparent that it could be used for distances or torsions as well, so
    // make the required number of atoms in the query matches variable.
    uint32_t _required_number_of_atoms;

    // For each query, an accumulator of results.
    AngleHistogram * _acc_per_query;
    AngleHistogram _acc;

    // We can write a histogram for each query.
    IWString _stem_for_histogram;

  // Private functions.

    int Preprocess(Molecule& m);
    int Process(Molecule& m, int query_number, const Set_of_Atoms& embedding,
                IWString_and_File_Descriptor& output);
    int BondAngle(Molecule& m, int query_number,
                    const Set_of_Atoms& embedding, IWString_and_File_Descriptor& output);
    int BondLength(Molecule& m, int query_number,
                    const Set_of_Atoms& embedding, IWString_and_File_Descriptor& output);
    int Torsion(Molecule& m, int query_number,
                    const Set_of_Atoms& embedding, IWString_and_File_Descriptor& output);

  public:
    BondAngles();
    ~BondAngles();

    int verbose() const {
      return _verbose;
    }

    int Initialise(Command_Line& cl);

    int Process(Molecule& m, IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

BondAngles::BondAngles() {
  _verbose = 0;
  _molecules_read = 0;
  _reduce_to_largest_fragment = 0;
  // hard coded for bond angles for now.
  _required_number_of_atoms = 3;
  _acc_per_query = nullptr;
}

BondAngles::~BondAngles() {
  if (_acc_per_query != nullptr) {
    delete [] _acc_per_query;
  }
}

int
BondAngles::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('g')) {
    if (! _chemical_standardisation.construct_from_command_line(cl, _verbose, 'g')) {
      cerr << "Cannot initialise chemical standardisation (-g)\n";
      return 0;
    }
  }

  if (cl.option_present('s')) {
    const_IWSubstring s;
    for (int i = 0; cl.value('s', s, i); ++i) {
      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
      if (! q->create_from_smarts(s)) {
        cerr << "Options::Initialise:cannot parse smarts '" << s << "'\n";
        return 0;
      }
      _queries.add(q.release());
    }
  }

  if (cl.option_present('q')) {
    if (! process_queries(cl, _queries, _verbose, 'q')) {
      cerr << "Cannot process queries (-q)\n";
      return 0;
    }
  }

  if (_queries.empty()) {
    cerr << "BondAngles::Initialise:no queries, use -s and/or -q\n";
    return 0;
  }

  if (_verbose) {
    cerr << "Defined " << _queries.size() << " bond angle queries\n";
  }

  for (Substructure_Query* q : _queries) {
    q->set_find_unique_embeddings_only(1);
  }

  for (uint32_t i = 0; i < _queries.size(); ++i) {
    IWString fname;
    fname << "/tmp/q" << i << ".qry";
    _queries[i]->write_msi(fname);
  }

  _acc_per_query = new AngleHistogram[_queries.size()];

  if (cl.option_present('x')) {
    if (! cl.value('x', _required_number_of_atoms) || _required_number_of_atoms < 2 || _required_number_of_atoms > 4) {
      cerr << "The number of atms matched by the query (-x) must be one of 2,3,4\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will require " << _required_number_of_atoms << " query atoms in each match\n";
    }
  }

  if (cl.option_present('S')) {
    cl.value('S', _stem_for_histogram);
    if (_verbose) {
      cerr << "Histograms written to '" << _stem_for_histogram << '\n';
    }
  }

  return 1;
}

int
BondAngles::Preprocess(Molecule& m) {
  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (m.natoms() == 0) {
    cerr << "Empty molecule\n";
    return 0;
  }

  return 1;
}

int
BondAngles::Process(Molecule& m,
                    IWString_and_File_Descriptor& output) {
  ++_molecules_read;
  if (! Preprocess(m)) {
    return 0;
  }

  Molecule_to_Match target(&m);

  const int nq = _queries.number_elements();
  for (int i = 0; i < nq; ++i) {
    Substructure_Results sresults;
    const int nhits = _queries[i]->substructure_search(target, sresults);
    if (nhits == 0) {
      continue;
    }

    for (const Set_of_Atoms* e : sresults.embeddings()) {
      if (e->size() < _required_number_of_atoms) {
        cerr << "Embedding size mismatch, got " << e->size() << " need " << _required_number_of_atoms << " ignored\n";
        continue;
      }

      if (! Process(m, i, *e, output)) {
        continue;
      }
    }
  }

  return 1;
}

int
BondAngles::Process(Molecule& m,
                    int query_number,
                    const Set_of_Atoms& embedding,
                    IWString_and_File_Descriptor& output) {
  switch (_required_number_of_atoms) {
    case 2:
      return BondLength(m, query_number, embedding, output);
    case 3:
      return BondAngle(m, query_number, embedding, output);
    case 4:
      return Torsion(m, query_number, embedding, output);
  }

  return 1;
}

int
BondAngles::BondAngle(Molecule& m,
                    int query_number,
                    const Set_of_Atoms& embedding,
                    IWString_and_File_Descriptor& output) {
  float angle = m.bond_angle(embedding[0], embedding[1], embedding[2]) * RAD2DEG;
  _acc_per_query[query_number].Extra(angle);
  _acc.Extra(angle);
  return 1;
}

int
BondAngles::BondLength(Molecule& m,
                    int query_number,
                    const Set_of_Atoms& embedding,
                    IWString_and_File_Descriptor& output) {
  cerr << "BondLength not implemented\n";
  return 1;
}

int
BondAngles::Torsion(Molecule& m,
                    int query_number,
                    const Set_of_Atoms& embedding,
                    IWString_and_File_Descriptor& output) {
  float angle = m.dihedral_angle(embedding[0], embedding[1], embedding[2], embedding[3]);
  _acc_per_query[query_number].Extra(angle);
  _acc.Extra(angle);
  return 1;
}
                

int
BondAngles::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  if (_molecules_read == 0) {
    return 1;
  }

  for (int i = 0; i < _queries.number_elements(); ++i) {
    output << "Query " << i << ' ' << _queries[i]->comment() << '\n';
    _acc_per_query[i].Report(output);
    if (_stem_for_histogram.empty()) {
      continue;
    }

    IWString fname;
    fname << _stem_for_histogram << '_' << i << ".csv";
    _acc_per_query[i].MakeHistogram(fname);

  }

  if (_queries.size() > 1) {
    output << "For all queries\n";
    _acc.Report(output);
    if (! _stem_for_histogram.empty()) {
      IWString fname;
      fname << _stem_for_histogram << "_all" << ".csv";
      _acc.MakeHistogram(fname);
    }
  }

  return 1;
}

int
GetBondAngles(Molecule& m,
              BondAngles& bond_angles,
              IWString_and_File_Descriptor& output) {
  return bond_angles.Process(m, output);
}

int
GetBondAngles(data_source_and_type<Molecule>& input,
              BondAngles& bond_angles,
              IWString_and_File_Descriptor& output) {
  Molecule * m;
  while (nullptr != (m = input.next_molecule()))
  {
    std::unique_ptr<Molecule> free_m(m);
    if (! GetBondAngles(*m, bond_angles, output)) {
      cerr << "Error processing " << m->name() << '\n';
      return 0;
    }
  }

  return 1;
}

int
GetBondAngles(const char* fname,
              FileType input_type,
              BondAngles& bond_angles,
              IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok()) {
    cerr << "cannot open '" << fname << "'\n";
    return 0;
  }

  if (bond_angles.verbose()) {
    cerr << "Processing '" << fname << "'\n";
  }

  return GetBondAngles(input, bond_angles, output);
}

int
GetBondAngles(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:E:i:q:s:x:S:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (! process_standard_aromaticity_options(cl, verbose, 'A')) {
    cerr << "Cannot process aromaticity options\n";
    return 1;
  }

  if (! process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process standard elements options (-E)\n";
    return 1;
  }

  BondAngles bond_angles;
  if (! bond_angles.Initialise(cl)) {
    cerr << "Cannot initialise\n";
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (!cl.option_present('i')) {
    if (1 == cl.number_elements() && 0 == strcmp("-", cl[0]))  // reading a pipe, assume smiles
      input_type = FILE_TYPE_SMI;
    else if (!all_files_recognised_by_suffix(cl)) {
      cerr << "Cannot discern all file types, use the -i option\n";
      return 4;
    }
  } else if (!process_input_type(cl, input_type)) {
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  for (const char * fname : cl) {
    if (! GetBondAngles(fname, input_type, bond_angles, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    bond_angles.Report(cerr);
  }

  return 0;
}

}  // namespace get_bond_angles

int
main(int argc, char** argv) {
  int rc = get_bond_angles::GetBondAngles(argc, argv);

  return rc;
}
