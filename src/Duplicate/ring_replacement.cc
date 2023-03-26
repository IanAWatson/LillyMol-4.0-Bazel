// Using queries generated by ring_extraction, replace rings in molecules.

#include <iostream>

#include "google/protobuf/io/zero_copy_stream_impl_lite.h"
#include "google/protobuf/text_format.h"

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Duplicate/replacement_ring.pb.h"

namespace ring_replacement {

using std::cerr;

void
Usage(int rc) {
  cerr << " -R <fname>    file of labelled rings created by ring_extraction\n";
  cerr << " -s <smarts>   only replace rings matched by <smarts>\n";
  cerr << " -q <query>    only replace rings matched by <query>\n";
  cerr << " -u            unique molecules only\n";
  cerr << " -a            allow change of aromaticity on ring replacement\n";
  cerr << " -p            write parent molecule\n";
  cerr << " -n <ex>       only process replacement rings with <ex> or more examples\n";
  cerr << " -Y <query>    product molecules MUST     match a query in <query>\n";
  cerr << " -N <query>    product molecules must NOT match any queries in <query>\n";
  cerr << " -D <query>    discard any replacement ring that matches <query>\n";
  cerr << " -I .          remove isotopes from product molecules\n";
  cerr << " -I <n>        change all existing isotopes to <n> (useful if atom types used)\n";
  cerr << " -w            sort output by precedent count\n";
  cerr << " -X ...        more options\n";
  cerr << " -c            remove chirality\n";
  cerr << " -l            strip to largest fragment\n";
  cerr << " -v            verbose output\n";

  ::exit(rc);
}

class Replacement {
  private:
    // Constructed from `smt` in the proto.
    Substructure_Query _query;

    // Constructed from `smi` in the proto.
    Molecule _m;

    // Copied from `id` in the proto.
    IWString _id;

    // the kind of ring replacement this is, 5a6a...
    IWString _label;

    // Copied from `n` in the proto. The number of examplars.
    int _count;

    // save for debugging.
    IWString _smarts;

  public:
    int BuildFromProto(const RplRing::ReplacementRing& proto);

    int count() const {
      return _count;
    }

    int SubstructureSearch(Molecule_to_Match& target, Substructure_Results& sresults);

    int MakeVariant(Molecule& m, const Set_of_Atoms& embedding, const uint32_t* atypes, std::unique_ptr<Molecule>& result) const;
};

int
Replacement::BuildFromProto(const RplRing::ReplacementRing& proto) {
  IWString smt(proto.smt());
  if (! _query.create_from_smarts(smt)) {
    cerr << "Replacement::BuildFromProto:cannot parse '" << proto.smt() << "'\n";
    return 0;
  }

  if (! _m.build_from_smiles(proto.smi())) {
    cerr << "Replacement::BuildFromProto:cannot parse smiles:" << proto.ShortDebugString() << "'n";
    return 0;
  }

  _id = proto.id();
  _label = proto.label();

  _count = proto.n();

  // _smarts = proto.smt();

  return 1;
}

int
Replacement::SubstructureSearch(Molecule_to_Match& target, Substructure_Results& sresults) {
  return _query.substructure_search(target, sresults);
}

void
UnsetImplicitHydrogenInformation(Molecule& m) {
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    m.unset_all_implicit_hydrogen_information(i);
  }
}

int
EnvironmentMatch(const Molecule& m,
                 atom_number_t zatom,
                 atom_number_t to_join,
                 const uint32_t* atypes) {
  if (! atypes) {
    return 1;
  }

  const isotope_t iso = m.isotope(zatom);
  // cerr << " atom " << zatom << " iso " << m.isotope(zatom) << " to_join " << to_join << " atype " << atypes[to_join] << " eq " << (iso == atypes[to_join]) << '\n';
  return iso == atypes[to_join];
}

int
Replacement::MakeVariant(Molecule& parent,
                         const Set_of_Atoms& embedding,
                         const uint32_t* atypes,
                         std::unique_ptr<Molecule>& result) const {
  std::unique_ptr<Molecule> m = std::make_unique<Molecule>(parent);
  const int initial_natoms = m->natoms();
  const int initial_aromatic_count = parent.aromatic_atom_count();

  m->set_name(parent.name());
  m->add_molecule(&_m);
  Molecule mcopy(_m);
  //cerr << mcopy.smiles() << " was added " << _id << " smarts " << _smarts << '\n';
#ifdef DEBUG_MAKE_VARIANT
  cerr << "After adding replacement " << m->smiles() << " initial_natoms " << initial_natoms << '\n';
  write_isotopically_labelled_smiles(*m, false, cerr);
  cerr << '\n';
#endif

  std::unique_ptr<int[]> to_remove(new_int(m->natoms()));
  embedding.set_vector(to_remove.get(), 1);

  resizable_array<const Bond*> remove_bonds;
  remove_bonds.reserve(10);  // approximately.

  const int esize = embedding.number_elements();
  for (int i = 0; i < esize; ++i) {
    const atom_number_t a = embedding[i];
    const Atom& atom = m->atom(a);
    for (const Bond* b : atom) {
      atom_number_t o = b->other(a);
      if (to_remove[o]) {
        continue;
      }
      if (! b->is_single_bond()) {
        continue;
      }
      if (! EnvironmentMatch(_m, i, o, atypes)) {
        return 0;
      }
      // cerr << "From atom " << a << " make bond " << (initial_natoms + i) << " to " << o << '\n';
      remove_bonds << b;

      atom_number_t join_to = initial_natoms + i;
      // This happens with double bonds being removed and related.
      if (m->hcount(join_to) == 0) {
        return 0;
      }

      if (b->is_single_bond()) {
        m->add_bond(o, join_to, SINGLE_BOND);
      } else if (b->is_double_bond()) {
        m->add_bond(o, join_to, DOUBLE_BOND);
      } else {
        m->add_bond(o, join_to, TRIPLE_BOND);
      }
    }
  }

  // Note that `b` is destroyed by remove_bond_between_atoms
  for (const Bond* b: remove_bonds) {
    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();
    // cerr << "Removing bond btw " << a1 << ' ' << a2 << '\n';
    if (m->are_bonded(a1, a2)) {
      m->remove_bond_between_atoms(a1, a2);
    }
  }

  m->remove_atoms(to_remove.get());

  UnsetImplicitHydrogenInformation(*m);

  if (int final_aromatic_atoms = m->aromatic_atom_count(); final_aromatic_atoms < initial_aromatic_count) {
    cerr << "Warning loss of aromaticty " << parent.unique_smiles() << " vs " << m->unique_smiles() << '\n';
    cerr << initial_aromatic_count << " vs " << final_aromatic_atoms << " aromatic atoms\n";
  }

  if (! m->valence_ok()) {
    cerr << "Invalid valence '" << m->smiles() << ' ' << parent.name() << ' ' << _id << '\n';
    for (int i = 0; i < m->natoms(); ++i) {
      if (! m->valence_ok(i)) {
        cerr << " atom " << i << ' ' << m->smarts_equivalent_for_atom(i) << '\n';
      }
    }
    return 0;
  }

  IWString name;
  name << parent.name();
  name << " %% " << _id;
  if (! _label.empty()) {
    name << '.' << _label;
  }
  name << ' ' << _count;
  m->set_name(name);

  result.reset(m.release());

  return 1;
}

class RingReplacement {
  private:
    int _verbose;

    int _molecules_read;

    int _reduce_to_largest_fragment;

    int _remove_chirality;

    int _write_parent_molecule;

    Atom_Typing_Specification _atype;

    // One per -P option.
    int _nreplacements;
    resizable_array_p<Replacement>* _rings;

    resizable_array_p<Substructure_Query> _queries;
    resizable_array_p<Substructure_Query> _products_must_have;
    resizable_array_p<Substructure_Query> _products_must_not_have;

    int _molecules_discarded_by_must_have_queries = 0;
    int _molecules_discarded_by_must_not_have_queries = 0;

    int _molecules_with_invalid_valences_suppressed;

    int _unique_molecules_only = 0;

    int _min_examples_needed = 0;
    int _replacement_rings_discarded_for_count = 0;

    IW_STL_Hash_Set _seen;
    int _duplicates_suppressed;

    int _molecules_formed;

    extending_resizable_array<int> _number_variants;

    // Remove isotopes from product molecules.
    int _remove_isotopes;

    // Translate all isotopes to a single value.
    isotope_t _translate_isotopes;

    int _sort_output_by_precedent = 0;

    FileType _input_type;

    Chemical_Standardisation _chemical_standardisation;

  // Private functions.

    int ReadReplacementRings(IWString& fname, int ndx);
    int ReadReplacementRings(iwstring_data_source& input, int ndx);
    int OkSupport(const Replacement& r);
    int IsUnique(Molecule& m);
    int Process(resizable_array_p<Molecule>& mols, int ndx, const uint32_t* atypes);
    int Write(const resizable_array_p<Molecule>& mols, IWString_and_File_Descriptor& output);

  public:
    RingReplacement();
    ~RingReplacement();

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);

    int Process(Molecule& m, IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;

    FileType input_type() const {
      return _input_type;
    }
};

RingReplacement::RingReplacement() {
  _verbose = 0;
  _molecules_read = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;

  _write_parent_molecule = 0;

  _molecules_discarded_by_must_have_queries = 0;
  _molecules_discarded_by_must_not_have_queries = 0;

  _molecules_with_invalid_valences_suppressed = 0;

  _unique_molecules_only = 0;

  _min_examples_needed = 0;
  _replacement_rings_discarded_for_count = 0;

  _duplicates_suppressed = 0;

  _molecules_formed = 0;

  _remove_isotopes = 0;
  _translate_isotopes = 0;

  _sort_output_by_precedent = 0;

  _nreplacements = 0;
  _rings = nullptr;

  _input_type = FILE_TYPE_INVALID;
}

RingReplacement::~RingReplacement() {
  if (_rings != nullptr) {
    delete [] _rings;
  }
}

int
RingReplacement::Initialise(Command_Line& cl) {
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
    const_IWSubstring smarts;
    for (int i = 0; cl.value('s', smarts, i); ++i) {
      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
      if (! q->create_from_smarts(smarts)) {
        cerr << "Cannot parse smarts '" << smarts << "'\n";
        return 0;
      }
      _queries << q.release();
    }
  }

  if (cl.option_present('q')) {
    if (! process_queries(cl, _queries, _verbose, 'q')) {
      cerr << "Cannot read queries (-q)\n";
      return 0;
    }
  }

  for (Substructure_Query* q : _queries) {
    q->set_embeddings_do_not_overlap(1);
  }

  if (cl.option_present('Y')) {
    if (! process_queries(cl, _products_must_have, _verbose, 'M')) {
      cerr << "Cannot read products must have queries (-M)\n";
      return 0;
    }
  }

  if (cl.option_present('N')) {
    if (! process_queries(cl, _products_must_not_have, _verbose, 'N')) {
      cerr << "Cannot read products must not have queries (-X)\n";
      return 0;
    }
  }

  if (cl.option_present('u')) {
    _unique_molecules_only = 1;
    if (_verbose) {
      cerr << "Will only produce unique molecules\n";
    }
  }

  if (cl.option_present('p')) {
    _write_parent_molecule = 1;
    if (_verbose) {
      cerr << "Will write the parent molecule\n";
    }
  }

  if (cl.option_present('n')) {
    if (! cl.value('n', _min_examples_needed) || _min_examples_needed < 1) {
      cerr << "The min examples needed option (-n) must be a whole +ve number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will only consider rings with at least " << _min_examples_needed << " examples\n";
    }
  }

  // Make sure this is set before reading the replacement rings.
  if (cl.option_present('w')) {
    _sort_output_by_precedent = 1;
    if (_verbose) {
      cerr << "Will write new molecules in order of precedent\n";
    }
  }

  if (! cl.option_present('R')) {
    cerr << "Must specify one of more sets of replacement rings via the -R option\n";
    return 0;
  }

  _nreplacements = cl.option_count('R');

  _rings = new resizable_array_p<Replacement>[_nreplacements];

  if (cl.option_present('R')) {
    IWString r;
    for (int i = 0; cl.value('R', r, i); ++i) {
      if (! ReadReplacementRings(r, i)) {
        cerr << "Cannot read replpcement rings '" << r << "'\n";
        return 0;
      }
    }
  }

  if (cl.option_present('I')) {
    const_IWSubstring i = cl.string_value('I');
    if (i == '.') {
      _remove_isotopes = 1;
      if (_verbose) {
        cerr << "Will remove isotopes from product molecules\n";
      }
    } else if (i.numeric_value(_translate_isotopes)) {
      if (_verbose) {
        cerr << "Will translate all isotopes to " << _translate_isotopes << '\n';
      }
    } else {
      cerr << "Unrecognised -I qualifier '" << i << "'\n";
      return 1;
    }
  }

  if (cl.option_present('P')) {
    const_IWSubstring p = cl.string_value('P');
    if (! _atype.build(p)) {
      cerr << "Invalid atom type specification '" << p  << "'\n";
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
RingReplacement::ReadReplacementRings(IWString& fname, int ndx) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "RingReplacement::ReadReplacementRings:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadReplacementRings(input, ndx);
}

int
RingReplacement::ReadReplacementRings(iwstring_data_source& input, int ndx) {
  const_IWSubstring line;
  while (input.next_record(line)) {
    google::protobuf::io::ArrayInputStream zero_copy_array(line.data(), line.nchars());
    RplRing::ReplacementRing proto;
    if (!google::protobuf::TextFormat::Parse(&zero_copy_array, &proto)) {
      cerr << "RingReplacement:ReadReplacementRings:cannot parse proto " << line << '\n';
      return 0;
    }

    std::unique_ptr<Replacement> r = std::make_unique<Replacement>();
    if (! r->BuildFromProto(proto)) {
      cerr << "RingReplacement::ReadReplacementRings:cannot parse " << proto.ShortDebugString() << "\n";
      return 0;
    }

    if (! OkSupport(*r)) {
      continue;
    }

    _rings[ndx] << r.release();
  }

  if (_rings[ndx].empty()) {
    cerr << "RingReplacement::ReadReplacementRings:no replacement rings\n";
    return 0;
  }

  if (_sort_output_by_precedent) {
    _rings[ndx].iwqsort_lambda([](const Replacement* r1, const Replacement* r2) {
      if (r1->count() < r2->count()) {
        return 1;
      } else if (r1->count() > r2->count()) {
        return -1;
      } else {
        return 0;
      }
    });
  }

  return 1;
}

int
RingReplacement::OkSupport(const Replacement& r) {
  if (r.count() < _min_examples_needed) {
    ++_replacement_rings_discarded_for_count;
    return 0;
  }

  return 1;
}

int
RingReplacement::Preprocess(Molecule& m) {
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
RingReplacement::Report(std::ostream& output) const {
  output << "RingReplacement:read " << _molecules_read << " molecules\n";
  if (_molecules_read == 0) {
    return 1;
  }
  if (_min_examples_needed) {
    output << _replacement_rings_discarded_for_count << " replacement rings discarded for support below " << _min_examples_needed << '\n';
  }
  output << _duplicates_suppressed << " duplicates suppressed\n";
  output << "Formed " << _molecules_formed << " molecules\n";

  for (int i = 0; i < _number_variants.number_elements(); ++i) {
    if (_number_variants[i]) {
      output << _number_variants[i] << " molecules generated " << i << " variants\n";
    }
  }

  return 1;
}

int
RingReplacement::Process(Molecule& m,
                         IWString_and_File_Descriptor& output) {
  ++_molecules_read;


  std::unique_ptr<uint32_t[]> atypes;
  if (_atype.active()) {
    atypes.reset(new uint32_t[m.natoms()]);
    _atype.assign_atom_types(m, atypes.get());
  }

  resizable_array_p<Molecule> generated;
  generated << new Molecule(m);

  if (! Process(generated, 0, atypes.get())) {
    return 0;
  }

  return Write(generated, output);
}

void
TranslateIsotopes(Molecule& m, isotope_t iso) {
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (m.isotope(i)) {
      m.set_isotope(i, iso);
    }
  }
}

int
RingReplacement::Process(resizable_array_p<Molecule>& mols,
                         int ndx,
                         const uint32_t* atypes) {

  resizable_array_p<Molecule> generated_here;
  generated_here.reserve(10 * mols.size());

  for (Molecule* m : mols) {
    Molecule_to_Match target(m);
    Substructure_Results sresults;

    for (Replacement* r : _rings[ndx]) {
      if (! r->SubstructureSearch(target, sresults)) {
        continue;
      }

      for (const Set_of_Atoms* e : sresults.embeddings()) {
        std::unique_ptr<Molecule> variant;
        if (! r->MakeVariant(*m, *e, atypes, variant)) {
          continue;
        }
        if (!IsUnique(*variant)) {
          continue;
        }
        if (_remove_isotopes) {
          variant->transform_to_non_isotopic_form();
        } else if (_translate_isotopes) {
          TranslateIsotopes(*variant, _translate_isotopes);
        }
        variant->reduce_to_largest_fragment();

        generated_here << variant.release();
      }
    }
  }

  ++_number_variants[generated_here.size()];

  mols.transfer_in(generated_here);

  if (ndx < _nreplacements - 1) {
    return Process(mols, ndx + 1, atypes);
  }

  return 1;
}

int
RingReplacement::IsUnique(Molecule& m) {
  if (! _unique_molecules_only) {
    return 1;
  }

  if (_seen.contains(m.unique_smiles())) {
    return 0;
  }

  _seen.insert(m.unique_smiles());

  return 1;
}

int
RingReplacement::Write(const resizable_array_p<Molecule>& mols,
                       IWString_and_File_Descriptor& output) {
  for (Molecule* m : mols) {
    output << m->smiles() << ' ' << m->name() << '\n';
    output.write_if_buffer_holds_more_than(8192);
  }
  output.flush();

  _molecules_formed += mols.number_elements();

  return 1;
}


int
ReplaceCore(RingReplacement& ring_replacement,
            Molecule& m,
            IWString_and_File_Descriptor& output) {
  return ring_replacement.Process(m, output);
}

int
ReplaceCore(RingReplacement& ring_replacement,
            data_source_and_type<Molecule>& input,
            IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! ring_replacement.Preprocess(*m)) {
      return 0;
    }

    if (! ReplaceCore(ring_replacement, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
ReplaceCore(RingReplacement& ring_replacement,
            const char * fname,
            IWString_and_File_Descriptor& output) {
  FileType input_type = ring_replacement.input_type();
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return ReplaceCore(ring_replacement, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:i:g:lcY:N:R:P:s:q:uapn:D:I:wX:");
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

  RingReplacement ring_replacement;
  if (! ring_replacement.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);
  for (const char* fname : cl) {
    if (! ReplaceCore(ring_replacement, fname, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    ring_replacement.Report(cerr);
  }

  return 0;
}

}  // namespace ring_replacement

int
main(int argc, char** argv) {
  int rc = ring_replacement::Main(argc, argv);

  return rc;
}
