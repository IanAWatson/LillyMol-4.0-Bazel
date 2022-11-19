//  Fingerprint a subset of a molecule defined by a query.

#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_pair_fingerprint.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/iwmfingerprint.h"
#include "Molecule_Lib/linear_fingerprint.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"
#include "Molecule_Tools/ec_fingerprint.h"
#include "Molecule_Tools/fingerprint_writer.h"

namespace molecule_subset_fingerprint {

using std::cerr;

IWString smiles_tag = "$SMI<";
IWString identifier_tag = "PCN<";

void
Usage(int rc)
{
  cerr << " -s <smarts>       specify query atoms via smarts\n";
  cerr << " -q <query>        specify query atoms via query\n";
  cerr << " -R <radius>       include atoms within <radius> of matched atoms\n";
  cerr << " -J <tag>          generate fingerprints with tag <tag>\n";
  cerr << " -F <fptype>       Fingperint type: linear, ec, ap\n";
  cerr << " -P <atype>        atom type specification\n";
  cerr << " -l                strip to largest fragment\n";
  cerr << " -v                verbose output\n";

  ::exit(rc);
}

class FingerprintSubset {
 private:

  int _verbose;

  int _molecules_read;

  int _reduce_to_largest_fragment;

  int _remove_chirality;

  int _radius;

  int _ignore_molecules_not_matching_any_queries;

  resizable_array_p<Substructure_Query> _queries;

  //  The different fingerprints we can generate.
  fingerprint_writer::FingerprintWriter _fingerprint_writer;
  ec_fingerprint::ECFingerprint _ec_fingerprint;
  ec_fingerprint::JobParameters _ec_parameters;
  atom_pair_fingerprint::AtomPairFingerprint _atom_pair_fingerprint;
  IWMFingerprint _iwfp;

  enum class FingerprintType { kUndef = 0, kLinear = 1, kEc = 2, kAtomPair = 3, kIwfp = 4 };

  FingerprintType _fptype;

  linear_fingerprint::LinearFingerprintGenerator _linear_fp;

  Atom_Typing_Specification _atom_typing;

  FileType _input_type;

  Chemical_Standardisation _chemical_standardisation;

  extending_resizable_array<int> _bits_set;

  //  Private functions.

  int Process(Molecule& m, const Set_of_Atoms& embedding, IWString_and_File_Descriptor& output);
  int Expand(Molecule& m, const Set_of_Atoms& embedding, int* include_atom) const;

  int LinearFingerprint(Molecule& m, const uint32_t* atype, const int* in_subset,
                        IWString_and_File_Descriptor& output);
  int AtomPairFingerprint(Molecule& m, const uint32_t* atype, const int* in_subset,
                          IWString_and_File_Descriptor& output);
  int ECFingerprint(Molecule& m, const uint32_t* atype, const int* in_subset,
                    IWString_and_File_Descriptor& output);
  int IWFP(Molecule& m, const uint32_t* atype, const int* in_subset,
           IWString_and_File_Descriptor& output);
  int DoOutput(Molecule& m, const Sparse_Fingerprint_Creator& sfc,
               IWString_and_File_Descriptor& output);

 public:

  FingerprintSubset();

  int Initialise(Command_Line& cl);

  int Preprocess(Molecule& m);

  int Process(Molecule& m, IWString_and_File_Descriptor& output);

  int Report(std::ostream& output) const;

  FileType input_type() const { return _input_type; }
};

FingerprintSubset::FingerprintSubset()
{
  _verbose = 0;
  _molecules_read = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _radius = 0;
  _fptype = FingerprintType::kLinear;
  _ignore_molecules_not_matching_any_queries = 0;
  _input_type = FILE_TYPE_INVALID;
}

int
FingerprintSubset::Initialise(Command_Line& cl)
{
  _verbose = cl.option_count('v');

  if (cl.option_present('c')) {
    _remove_chirality = 1;
    if (_verbose) {
      cerr << "Will remove chirality from input molecules\n";
    }
  }

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      cerr << "Cannot initialise chemical standardisation\n";
      return 0;
    }
  }

  if (!cl.option_present('J')) {
    cerr << "Must specify fingerprint via the -J option\n";
    return 0;
  }

  if (!_fingerprint_writer.Initialise(cl, 'J', _verbose)) {
    cerr << "Cannot initialise fingerprint writer\n";
    return 0;
  }

  if (cl.option_present('R')) {
    if (!cl.value('R', _radius) || _radius < 0) {
      cerr << "The bond distance radius (-R) option must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will include atoms within " << _radius << " bonds of matched atoms\n";
    }
  }

  if (cl.option_present('s')) {
    const_IWSubstring s;
    for (int i = 0; cl.value('s', s, i); ++i) {
      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
      if (!q->create_from_smarts(s)) {
        cerr << "Options::Initialise:cannot parse smarts '" << s << "'\n";
        return 0;
      }
      _queries << q.release();
    }
  }

  if (cl.option_present('q')) {
    if (!process_queries(cl, _queries, _verbose, 'q')) {
      cerr << "Cannot read queries (-q)\n";
      return 0;
    }
  }

  if (_queries.empty()) {
    cerr << "No queries\n";
    return 0;
  }

  if (cl.option_present('z')) {
    const_IWSubstring z;
    for (int i = 0; cl.value('z', z, i); ++i) {
      if (z == 'i') {
        _ignore_molecules_not_matching_any_queries = 1;
        if (_verbose) {
          cerr << "Will ignore molecules not matching any queries\n";
        }
      }
    }
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce molecules to largest fragment\n";
    }
  }

  if (cl.option_present('P')) {
    const_IWSubstring p = cl.string_value('P');
    if (!_atom_typing.build(p)) {
      cerr << "Cannot parse atom typing '-P " << p << "'\n";
      return 0;
    }
  }
  else {
    const_IWSubstring p = "UST:AY";
    _atom_typing.build(p);
  }

  if (cl.option_present('F')) {
    const_IWSubstring f = cl.string_value('F');
    if (f == "linear") {
      _fptype = FingerprintType::kLinear;
    }
    else if (f == "ec") {
      _fptype = FingerprintType::kEc;
    }
    else if (f == "ap") {
      _fptype = FingerprintType::kAtomPair;
    }
    else {
      cerr << "Unrecognised -F directive " << f << "'\n";
      return 0;
    }
  }
  else {
    _fptype = FingerprintType::kLinear;
  }

  if (1 == cl.number_elements() && 0 == strcmp("-", cl[0])) {  //  reading a pipe, assume smiles
    _input_type = FILE_TYPE_SMI;
  }
  else if (!all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot discern all file types, use the -i option\n";
    return 0;
  }
  else if (!process_input_type(cl, _input_type)) {
    return 0;
  }

  return 1;
}

int
FingerprintSubset::Preprocess(Molecule& m)
{
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
FingerprintSubset::Process(Molecule& m, IWString_and_File_Descriptor& output)
{
  Molecule_to_Match target(&m);
  Substructure_Results sresults;

  int got_match = 0;

  for (Substructure_Query* q : _queries) {
    int nhits = q->substructure_search(target, sresults);
    if (nhits == 0) {
      continue;
    }

    for (const Set_of_Atoms* e : sresults.embeddings()) {
      Process(m, *e, output);
      ++got_match;
    }
  }

  if (got_match) {
    return 1;
  }

  if (_ignore_molecules_not_matching_any_queries) {
    return 1;
  }

  cerr << "FingerprintSubset::Process:No match to " << m.name() << '\n';
  return 0;
}

//  For the atoms in `m`, if they are within _radius bonds
//  of an atom in `embedding`, set the corresponding entry
//  in `remove_atom`.
int
FingerprintSubset::Expand(Molecule& m, const Set_of_Atoms& embedding, int* include_atom) const
{
  const int matoms = m.natoms();
  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (include_atom[i]) {
      continue;
    }
    for (atom_number_t j : embedding) {
      if (m.bonds_between(i, j) <= _radius) {
        include_atom[i] = 1;
        ++rc;
        break;
      }
    }
  }

  return rc;
}

int
FingerprintSubset::Process(Molecule& m, const Set_of_Atoms& embedding,
                           IWString_and_File_Descriptor& output)
{
  ++_molecules_read;

  const int matoms = m.natoms();
  std::unique_ptr<int[]> in_subset(new_int(matoms));
  embedding.set_vector(in_subset.get(), 1);
  if (_radius > 0) {
    Expand(m, embedding, in_subset.get());
  }

  std::unique_ptr<uint32_t[]> atype = std::make_unique<uint32_t[]>(matoms);
  _atom_typing.assign_atom_types(m, atype.get());

  switch (_fptype) {
    case FingerprintType::kUndef:
      cerr << "FingerprintSubset::Process:no fptype\n";
      return 0;
    case FingerprintType::kLinear:
      return LinearFingerprint(m, atype.get(), in_subset.get(), output);
    case FingerprintType::kEc:
      return ECFingerprint(m, atype.get(), in_subset.get(), output);
    case FingerprintType::kAtomPair:
      return AtomPairFingerprint(m, atype.get(), in_subset.get(), output);
    case FingerprintType::kIwfp:
      return IWFP(m, atype.get(), in_subset.get(), output);
  }

  cerr << "Should not come to here\n";
  return 0;
}

int
FingerprintSubset::LinearFingerprint(Molecule& m, const uint32_t* atype, const int* in_subset,
                                     IWString_and_File_Descriptor& output)
{
  Sparse_Fingerprint_Creator sfc;

  const int matoms = m.natoms();

  std::unique_ptr<linear_fingerprint::atom_type_t[]> tmp =
      std::make_unique<linear_fingerprint::atom_type_t[]>(matoms);
  for (int i = 0; i < matoms; ++i) {
    tmp[i] = atype[i];
  }

  _linear_fp.Fingerprint(m, in_subset, tmp.get(), sfc);

  return DoOutput(m, sfc, output);
}

int
FingerprintSubset::ECFingerprint(Molecule& m, const uint32_t* atype, const int* in_subset,
                                 IWString_and_File_Descriptor& output)
{
  ec_fingerprint::ProduceFingerprint produce_fp;
  _ec_fingerprint.Fingerprint(m, in_subset, atype, produce_fp);
  const Sparse_Fingerprint_Creator& sfc = produce_fp.sfc();
  return DoOutput(m, sfc, output);
}

int
FingerprintSubset::AtomPairFingerprint(Molecule& m, const uint32_t* atype, const int* in_subset,
                                       IWString_and_File_Descriptor& output)
{
  Sparse_Fingerprint_Creator sfc;
  _atom_pair_fingerprint.Fingerprint(m, in_subset, atype, sfc);

  return DoOutput(m, sfc, output);
}

int
FingerprintSubset::DoOutput(Molecule& m, const Sparse_Fingerprint_Creator& sfc,
                            IWString_and_File_Descriptor& output)
{
  if (! _fingerprint_writer.IsWritingDescriptors()) {
    output << smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
  }

  _fingerprint_writer.WriteFingerprint(m.name(), sfc, output);

  if (! _fingerprint_writer.IsWritingDescriptors()) {
    output << "|\n";
  }

  ++_bits_set[sfc.nbits()];

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

//  iwfp does not conform to the newer API's, not implemented.
int
FingerprintSubset::IWFP(Molecule& m, const uint32_t* atype, const int* in_subset,
                        IWString_and_File_Descriptor& output)
{
  return 1;
}

int
FingerprintSubset::Report(std::ostream& output) const
{
  output << "FingerprintSubset:read " << _molecules_read << " molecules\n";
  if (_molecules_read == 0) {
    return 1;
  }

  Accumulator_Int<int> acc;

  for (int i = 0; i < _bits_set.number_elements(); ++i) {
    if (_bits_set[i]) {
      acc.extra(i);
      if (_verbose > 1) {
        output << _bits_set[i] << " molecules had " << i << " bits set\n";
      }
    }
  }
  output << "Set btw " << acc.minval() << " and " << acc.maxval() << " ave " << acc.average()
         << '\n';

  return 1;
}

int
MoleculeSubsetFingerprint(FingerprintSubset& fingerprint_subset, Molecule& m,
                          IWString_and_File_Descriptor& output)
{
  if (!fingerprint_subset.Process(m, output)) {
    return 0;
  }

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
MoleculeSubsetFingerprint(FingerprintSubset& fingerprint_subset,
                          data_source_and_type<Molecule>& input,
                          IWString_and_File_Descriptor& output)
{
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (!fingerprint_subset.Preprocess(*m)) {
      return 0;
    }

    if (!MoleculeSubsetFingerprint(fingerprint_subset, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
MoleculeSubsetFingerprint(FingerprintSubset& fingerprint_subset, const char* fname,
                          IWString_and_File_Descriptor& output)
{
  FileType input_type = fingerprint_subset.input_type();
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return MoleculeSubsetFingerprint(fingerprint_subset, input, output);
}

int
Main(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vE:A:i:g:lcJ:R:s:q:P:F:z:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose, 'A')) {
    cerr << "Cannot process aromaticity options\n";
    return 1;
  }

  if (!process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process standard elements options (-E)\n";
    return 1;
  }

  FingerprintSubset fingerprint_subset;
  if (!fingerprint_subset.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);
  for (const char* fname : cl) {
    if (!MoleculeSubsetFingerprint(fingerprint_subset, fname, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    fingerprint_subset.Report(cerr);
  }

  return 0;
}

}  //  namespace molecule_subset_fingerprint

int
main(int argc, char** argv)
{
  int rc = molecule_subset_fingerprint::Main(argc, argv);

  return rc;
}
