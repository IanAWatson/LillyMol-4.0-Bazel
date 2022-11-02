// Identify substituents defined by a substructure query.

#include <iostream>
#include <limits>
#include <memory>
#include <string>

#include "google/protobuf/text_format.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Duplicate/dicer_fragment.pb.h"

namespace get_substituent {

using std::cerr;

void
Usage(int rc) {
  cerr << "Identify substituents defined by a substructure query\n";
  cerr << " -s <smarts>   define the query via one or more smarts\n";
  cerr << " -q <query>    define the query via one or more query files\n";
  cerr << " -m <natoms>   minimum number of atoms in each substituent\n";
  cerr << " -M <natoms>   maximum number of atoms in each substituent\n";
  cerr << " -X <length>   maximum distance from the attachment point\n";
  cerr << " -. <smarts>   queries that must     be present in the substituent\n";
  cerr << " -. <smarts>   queries that must not be present in the substituent\n";
  cerr << " -z i          ignore molecules not matching the query\n";
  cerr << " -I <iso>      isotopically label substituents\n";
  cerr << " -Y <smarts>   substituents must match\n";
  cerr << " -N <smarts>   substituents must not match\n";
  cerr << " -c                remove chirality\n";
  cerr << " -l                strip to largest fragment\n";
  cerr << " -v                verbose output\n";

  ::exit(rc);
}

class GetSubstituent {
  private:
    int _verbose;

    int _molecules_read;

    int _reduce_to_largest_fragment;

    int _remove_chirality;

    int _min_natoms;
    int _max_natoms;

    int _max_length;

    resizable_array_p<Substructure_Query> _query;

    resizable_array_p<Substructure_Query> _must_have;
    resizable_array_p<Substructure_Query> _must_not_have;

    int _ignore_molecules_not_matching_query;
    int _molecules_not_matching_query;

    int _isotope;

    IW_STL_Hash_Map<IWString, Dicer::DicerFragment> _fragments;

    extending_resizable_array<int> _fragment_size;

    FileType _input_type;

    Chemical_Standardisation _chemical_standardisation;

  // Private functions.
    int OkSize(const Molecule& m) const;
    int OkQueries(Molecule& m);
    int OkLength(Molecule& m, atom_number_t anchor) const;

    int Process(Molecule& m,
                const Substructure_Results& sresults,
                IWString_and_File_Descriptor& output);
    int Process(Molecule& m,
                const Set_of_Atoms& embedding,
                int* storage,
                int* xref,
                IWString_and_File_Descriptor& output);

  public:
    GetSubstituent();

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);

    int Process(Molecule& m, IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;

    int WriteFragments(IWString_and_File_Descriptor& output) const;

    FileType input_type() const {
      return _input_type;
    }
};

GetSubstituent::GetSubstituent() {
  _verbose = 0;
  _molecules_read = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;

  _min_natoms = 0;
  _max_natoms = std::numeric_limits<int>::max();

  _max_length = std::numeric_limits<int>::max();

  _ignore_molecules_not_matching_query = 0;
  _molecules_not_matching_query = 0;

  _isotope = 0;
  _input_type = FILE_TYPE_INVALID;
}

int
GetSubstituent::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('c')) {
    _remove_chirality = 1;
    if (_verbose) {
      cerr << "Will remove chirality from input molecules\n";
    }
  }

  if (cl.option_present('g')) {
    if (! _chemical_standardisation.construct_from_command_line(cl, _verbose, 'g')) {
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

  if (cl.option_present('s')) {
    const_IWSubstring s;
    for (int i = 0; cl.value('s', s, i); ++i) {
      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
      if (! q->create_from_smarts(s)) {
        cerr << "Invalid smarts '" << s << "'\n";
        return 0;
      }
      _query << q.release();
    }
  }

  if (cl.option_present('q')) {
    if (! process_queries(cl, _query, _verbose, 'q')) {
      cerr << "Cannot process command line queries (-q)\n";
      return 0;
    }
  }

  if (_query.empty()) {
    cerr << "No queries\n";
    return 0;
  }

  if (cl.option_present('z')) {
    _ignore_molecules_not_matching_query = 1;
    if (_verbose) {
      cerr << "Will ignore molecules not matching the query\n";
    }
  }

  if (cl.option_present('m')) {
    if (! cl.value('m', _min_natoms) || _min_natoms < 1) {
      cerr << "The minimum fragment size (-m) must be a whole +ve number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will only write fragments with at least " << _min_natoms << " atoms\n";
    }
  }

  if (cl.option_present('M')) {
    if (! cl.value('M', _max_natoms) || _max_natoms < _min_natoms) {
      cerr << "The maximum fragment size (-M) must be a whole +ve number greater than " << _min_natoms << '\n';
      return 0;
    }
    if (_verbose) {
      cerr << "Will only write fragments with at most " << _max_natoms << " atoms\n";
    }
  }

  if (cl.option_present('I')) {
    if (! cl.value('I', _isotope) || _isotope < 0) {
      cerr << "The isotopic label to apply (-I) must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will label attachment points with isotope " << _isotope << '\n';
    }
  }

  if (cl.option_present('Y')) {
    const_IWSubstring y;
    for (int i = 0; cl.value('Y', y, i); ++i) {
      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
      if (! q->create_from_smarts(y)) {
        cerr << "Invalid smarts '" << y << "'\n";
        return 0;
      }
      _must_have << q.release();
    }
  }

  if (cl.option_present('N')) {
    const_IWSubstring n;
    for (int i = 0; cl.value('N', n, i); ++i) {
      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
      if (! q->create_from_smarts(n)) {
        cerr << "Invalid smarts '" << n << "'\n";
        return 0;
      }
      _must_not_have << q.release();
    }
  }

  if (cl.option_present('L')) {
    if (! cl.value('L', _max_length) || _max_length < 1) {
      cerr << "The max length option (-L) must be a whole +ve number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will only find fragments shorter than " << _max_length << " bonds from the attachment point\n";
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
GetSubstituent::Preprocess(Molecule& m) {
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
GetSubstituent::Report(std::ostream& output) const {
  output << "GetSubstituent:read " << _molecules_read << " molecules\n";
  if (_molecules_read == 0) {
    return 1;
  }

  output << _molecules_not_matching_query << " molecules did not match the query\n";

  for (int i = 0; i < _fragment_size.number_elements(); ++i) {
    if (_fragment_size[i]) {
      output << _fragment_size[i] << " fragments had " << i << " atoms\n";
    }
  }

  return 1;
}

int
GetSubstituent::WriteFragments(IWString_and_File_Descriptor& output) const {
  static google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);

  std::string buffer;

  for (const auto& [usmi, frag] : _fragments) {
    printer.PrintToString(frag, &buffer);
    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

int
GetSubstituent::Process(Molecule& m,
                        IWString_and_File_Descriptor& output) {
  ++_molecules_read;

  Molecule_to_Match target(&m);
  Substructure_Results sresults;
  for (Substructure_Query* q : _query) {
    const int nhits = q->substructure_search(target, sresults);
    if (nhits == 0) {
      continue;
    }

    return Process(m, sresults, output);
  }

  ++_molecules_not_matching_query;
  return _ignore_molecules_not_matching_query;
}

int
GetSubstituent::Process(Molecule& m,
                const Substructure_Results& sresults,
                IWString_and_File_Descriptor& output) {
  const int matoms = m.natoms();

  std::unique_ptr<int[]> storage = std::make_unique<int[]>(matoms);
  // Xref is only needed if we are applying isotopes or if _max_length used.
  std::unique_ptr<int[]> xref = std::make_unique<int[]>(matoms);

  (void) m.ring_membership();

  for (const Set_of_Atoms* e : sresults.embeddings()) {
    Process(m, *e, storage.get(), xref.get(), output);
  }

  return 1;
}

int
IdentifySubstituent(const Molecule& m,
                    atom_number_t zatom,
                    int * visited) {
  visited[zatom] = 1;
  int rc = 1;
  const Atom& a = m.atom(zatom);
  for (const Bond* b : a) {
    atom_number_t o = b->other(zatom);
    if (visited[o]) {
      continue;
    }

    rc += IdentifySubstituent(m, o, visited);
  }

  return rc;
}

int
GetSubstituent::OkSize(const Molecule& m) const {
  const int matoms = m.natoms();
  if (matoms < _min_natoms) {
    return 0;
  }
  if (matoms > _max_natoms) {
    return 0;
  }

  return 1;
}

int
GetSubstituent::OkQueries(Molecule& m) {
  if (_must_have.empty() && _must_not_have.empty()) {
    return 1;
  }

  Molecule_to_Match target(&m);
  Substructure_Results sresults;
  int got_must_have = 0;
  for (Substructure_Query* q : _must_have) {
    if (q->substructure_search(target, sresults)) {
      got_must_have = 1;
      break;
    }
  }

  if (! _must_have.empty() && ! got_must_have) {
    return 0;
  }

  for (Substructure_Query* q : _must_not_have) {
    if (q->substructure_search(target, sresults)) {
      return 0;
    }
  }

  return 1;
}

int
GetSubstituent::OkLength(Molecule& m, atom_number_t anchor) const {
  if (_max_length == std::numeric_limits<int>::max()) {
    return 1;
  }

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (i == anchor) {
      continue;
    }
    if (m.bonds_between(anchor, i) > _max_length) {
      return 0;
    }
  }

  return 1;
}

int
GetSubstituent::Process(Molecule& m,
                const Set_of_Atoms& embedding,
                int* storage,
                int* xref,
                IWString_and_File_Descriptor& output) {
  if (embedding.size() < 2) {
    return 0;
  }

  const int matoms = m.natoms();

  const atom_number_t a0 = embedding[0];
  const atom_number_t a1 = embedding[1];
  if (m.in_same_ring(a0, a1)) {
    return 0;
  }

  std::fill_n(storage, matoms, 0);
  storage[a0] = 1;
  IdentifySubstituent(m, a1, storage);
  storage[a0] = 0;

  Molecule substituent;
  m.create_subset(substituent, storage, 1, xref);

  if (! OkSize(substituent)) {
    return 0;
  }

  if (_isotope > 0) {
    substituent.set_isotope(xref[a1], _isotope);
  }

  if (! OkQueries(substituent)) {
    return 0;
  }

  if (! OkLength(substituent, xref[a1])) {
    return 0;
  }

  ++_fragment_size[substituent.natoms()];

  const IWString& usmi = substituent.unique_smiles();
  auto iter = _fragments.find(usmi);
  if (iter == _fragments.end()) {
    Dicer::DicerFragment frag;
    frag.set_smi(usmi.AsString());
    frag.set_id(m.name().AsString());
    frag.set_n(1);
    _fragments.emplace(std::make_pair(usmi, std::move(frag)));
  } else {
    const auto n = iter->second.n();
    iter->second.set_n(n + 1) ;
  }

  return 1;
}

int
GetSubstituents(GetSubstituent& get_substituents,
            Molecule& m,
            IWString_and_File_Descriptor& output) {
  return get_substituents.Process(m, output);
}

int
GetSubstituents(GetSubstituent& get_substituents,
            data_source_and_type<Molecule>& input,
            IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! get_substituents.Preprocess(*m)) {
      return 0;
    }

    if (! GetSubstituents(get_substituents, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
GetSubstituents(GetSubstituent& get_substituents,
            const char * fname,
            IWString_and_File_Descriptor& output) {
  FileType input_type = get_substituents.input_type();
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return GetSubstituents(get_substituents, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:i:g:lcz:s:q:m:M:I:Y:N:X:");
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

  GetSubstituent get_substituents;
  if (! get_substituents.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);
  for (const char* fname : cl) {
    if (! GetSubstituents(get_substituents, fname, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  get_substituents.WriteFragments(output);
  output.flush();

  if (verbose) {
    get_substituents.Report(cerr);
  }

  return 0;
}

}  // namespace get_substituent

int
main(int argc, char** argv) {
  int rc = get_substituent::Main(argc, argv);

  return rc;
}
