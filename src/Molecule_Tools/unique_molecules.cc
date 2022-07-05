/*
  Sometimes we need to determine the unique molecules from a file
  of structures.
*/

#include <iostream>
#include <memory>
#include <stdlib.h>
#include <sys/stat.h>

#include "absl/container/flat_hash_set.h"
#include "absl/hash/hash.h"

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwstring/absl_hash.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/iwreaction.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/numass.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/rmele.h"
#include "Molecule_Lib/smiles.h"

using std::cerr;

const char* prog_name = nullptr;

class UniqueMoleculesImplementation {
  private:
    int _verbose;

    // The number of molecules examined vial IsUnique().
    int _molecules_processed;

    // The number of duplicates found by IsUnique().
    int _duplicates_found;

    // the number of unique molecules found by IsUnique().
    int _unique_molecules;

    // The various standardization actiona available.

    // If set, input molecules are transformed to their largest fragment.
    int _reduce_to_largest_fragment;

    // If set, chirality information is excluded from unique smiles formation.
    int _exclude_chiral_info;

    // If set, cis/trans information is excluded from unique smiles formation.
    int _exclude_cis_trans_bonding_info;

    // We can optionally ignore isotopic information.
    int _ignore_isotopes;

    Chemical_Standardisation _chemical_standardisation;

    // Useful for things like compressing the heavy halogens.
    Element_Transformations _element_transformations;

    // Specific elements can be removed before unique smiles generation.
    Elements_to_Remove _elements_to_remove;

    resizable_array_p<IWReaction> _reaction;
    int _molecules_changed_by_reactions;

    // A unique smiles hash for every atom count.
    resizable_array_p<absl::flat_hash_set<IWString>> _smiles_hash;

    // An array of atom hashes from quick_bond_hash();
    resizable_array_p<absl::flat_hash_set<uint64_t>> _atom_hash;
    int _use_atom_hash;

    // One possible comparison is to compare molecules via their graph form.
    int _compare_as_graph;

    // Timing results show that generally the molecular formula hash is not
    // worth computing.
    resizable_array_p<absl::flat_hash_set<IWString>> _formula_hash;
    int _perform_formula_check;

    // A variant on unique molecule finding is that molecules are unique only
    // if the name field also matches.
    int _only_same_if_structure_and_name_the_same;

    // This may have performance benefits.
    int _default_primary_hash_size;

    // private functions
    int PerformReactions(Molecule& m);
    int Preprocess(Molecule& m);
    int IsUniqueInner(Molecule& m);
    void FormUniqueSmiles(Molecule& m, IWString& usmi) const;
    void FormUniqueSmilesInner(Molecule& m, IWString& usmi) const;

  public:
    UniqueMoleculesImplementation();

    // Setter functions.
    void set_verbose(int s) {
      _verbose = s;
    }
    void set_reduce_to_largest_fragment(int s) {
      _reduce_to_largest_fragment = s;
    }
    void set_exclude_chiral_info(int s) {
      _exclude_chiral_info = s;
    }
    void set_exclude_cis_trans_bonding_info(int s) {
      _exclude_cis_trans_bonding_info = s;
    }
    void set_ignore_isotopes(int s) {
      _ignore_isotopes = s;
    }

    int Initialise(Command_Line& cl);

    // Add the information for `m` to the internal hashes.
    // It adds the unique smiles to the hash, and returns whether or not the molecule
    // is unique, but does not update any counters.
    int IngestPreviousMolecule(Molecule& m);

    // Returns 1 if `m` is unique.
    int IsUnique(Molecule& m);

    int Report(std::ostream& output) const;
};

UniqueMoleculesImplementation::UniqueMoleculesImplementation() {
  _verbose = 0;
  _molecules_processed = 0;
  _duplicates_found = 0;
  _unique_molecules = 0;
  _reduce_to_largest_fragment = 0;
  _exclude_chiral_info = 0;
  _exclude_cis_trans_bonding_info = 0;
  _ignore_isotopes = 0;

  _use_atom_hash = 0;
  _perform_formula_check = 0;
  _compare_as_graph = 0;
  _only_same_if_structure_and_name_the_same = 0;
  _default_primary_hash_size = 5000;

  _molecules_changed_by_reactions = 0;

  return;
}

int
UniqueMoleculesImplementation::Report(std::ostream& output) const {
  output << _molecules_processed << " molecules read, ";
  output << _duplicates_found << " duplicates, ";
  output << _unique_molecules << " unique structures\n";

  if (_use_atom_hash) {
    unsigned int s = 0;
    for (auto & h : _atom_hash) {
      s += h->size();
    }
    cerr << "Atom hash contains " << s << " discrete values\n";
  }

  return output.good();
}

// No checking, no NOTHING!
// Beware if reaction queries overlap!!
// Be very careful!
int
UniqueMoleculesImplementation::PerformReactions(Molecule& m) {
  int rc = 0;

  for (IWReaction* rxn : _reaction) {
    Substructure_Results sresults;

    int nhits = rxn->determine_matched_atoms(m, sresults);
    if (nhits == 0) {
      continue;
    }

    if (_verbose > 2) {
      cerr << nhits << " hits for reaction " << rxn->comment() << '\n';
    }

    rxn->in_place_transformations(m, sresults);

    ++rc;
  }

  if (rc && _verbose > 1) {
    cerr << m.name() << ", " << rc << " reactions changed\n";
  }

  if (rc) {
    _molecules_changed_by_reactions++;
  }

  return rc;
  return 1;
}

int
UniqueMoleculesImplementation::Preprocess(Molecule& m) {
  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (_elements_to_remove.active()) {
    _elements_to_remove.process(m);
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_element_transformations.active()) {
    _element_transformations.process(m);
  }

  if (_elements_to_remove.active()) {
    _elements_to_remove.process(m);
  }

  if (! _reaction.empty()) {
    PerformReactions(m);
  }

  if (_ignore_isotopes) {
    m.transform_to_non_isotopic_form();
  }

  return 1;
}

int
UniqueMoleculesImplementation::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce input molecules to largest fragment before usmi generation\n";
    }
  }

  if (cl.option_present('I')) {
    _ignore_isotopes = 1;
    if (_verbose) {
      cerr << "Will transform molecules to non isotopic form\n";
    }
  }

  if (cl.option_present('z')) {
    _exclude_cis_trans_bonding_info = 1;
    if (_verbose) {
      cerr << "Will exclude cis trans bonding from unique smiles generation\n";
    }
  }

  if (cl.option_present('c')) {
    _exclude_chiral_info = 1;
    if (_verbose) {
      cerr << "Will exclude chirality information from unique smiles generation\n";
    }
  }

  if (cl.option_present('h')) {
    _use_atom_hash = 1;
    if (_verbose) {
      cerr << "Will use an atom hash from quick_bond_hash()\n";
    }
  }

  if (cl.option_present('a')) {
    _compare_as_graph = 1;
    if (_verbose) {
      cerr << "Will compare graph forms\n";
    }
  }

  if (cl.option_present('j')) {
    _only_same_if_structure_and_name_the_same = 1;
    if (_verbose)
      cerr << "Molecules will be considered identical only if both structure and name match\n";
  }

  if (cl.option_present('s')) {
    if (! cl.value('s', _default_primary_hash_size) || _default_primary_hash_size < 1) {
      cerr << "The -s option must be followed by a positive whole number\n";
      return 0;
    }

    if (_verbose)
      cerr << "Primary hash size set to " << _default_primary_hash_size << '\n';
  }

  if (cl.option_present('R')) {
    Sidechain_Match_Conditions sidechain_match_conditions;

    if (! read_reactions(cl, _reaction, sidechain_match_conditions, 'R')) {
      cerr << "Cannot read reaction(s) (-R option)\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Defined " << _reaction.size() << " reactions\n";
    }

    for (IWReaction* r : _reaction) {
      r->set_find_unique_embeddings_only(1);
    }
  }

  if (cl.option_present('g')) {
    if (! _chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      cerr << "Cannot process chemical standardisations (-g option)\n";
      return 0;
    }
  }

  if (cl.option_present('X')) {
    if (! _elements_to_remove.construct_from_command_line(cl, _verbose, 'X'))
    {
      cerr << "Cannot discern elements to remove, -X\n";
      return 0;
    }
  }

  if (cl.option_present('T')) {
    if (! _element_transformations.construct_from_command_line(cl, _verbose, 'T')) {
      cerr << "Cannot process element transformations (-t option)\n";
      return 0;
    }
  }

  for (int i = 0; i < 200; i++) {
    _smiles_hash.add(new absl::flat_hash_set<IWString>());

    if (_perform_formula_check) {
      _formula_hash.add(new absl::flat_hash_set<IWString>());
    }
    if (_use_atom_hash) {
      _atom_hash.add(new absl::flat_hash_set<uint64_t>());
    }
  }

  return 1;
}

// No special provision for the empty molecule.
int
UniqueMoleculesImplementation::IsUnique(Molecule& m) {
  ++_molecules_processed;

  if (IsUniqueInner(m)) {
    ++_unique_molecules;
    return 1;
  } else {
    ++_duplicates_found;
    return 0;
  }
}

int
UniqueMoleculesImplementation::IngestPreviousMolecule(Molecule& m) {
  return IsUniqueInner(m);
}

// Form the key that will be used for the hash lookup.
void
UniqueMoleculesImplementation::FormUniqueSmilesInner(Molecule& m,
                        IWString& usmi) const {
  if (_compare_as_graph) {
    const IWString formula = m.molecular_formula();

    Molecule g(m);
    g.change_to_graph_form();

    usmi = g.unique_smiles();

    usmi << ':' << formula;
  } else {
    usmi = m.unique_smiles();
  }

  if (_only_same_if_structure_and_name_the_same) {
    usmi << m.name();
  }
  //cerr << "Testing unique smiles '" << usmi << "'\n";
}

void
UniqueMoleculesImplementation::FormUniqueSmiles(Molecule& m, IWString& usmi) const {
  const auto save_chiral = include_chiral_info_in_smiles();
  const auto save_cis_trans = include_cis_trans_in_smiles();

  if (_exclude_chiral_info) {
    set_include_chiral_info_in_smiles(0);
  }

  if (_exclude_cis_trans_bonding_info) {
    set_include_cis_trans_in_smiles(0);
  }

  FormUniqueSmilesInner(m, usmi);

  set_include_chiral_info_in_smiles(save_chiral);
  set_include_cis_trans_in_smiles(save_cis_trans);
}

#ifdef FORMULA_CHECK_SLOWS_THINGS_DOWN
// This never worked, seems it is too expensive to compute and check the formula.
// quick_bond_hash works better.
    IWString mformula;
    m.formula_distinguishing_aromatic(mformula);

    //cerr << m.name() << " usmi " << usmi << "'\n";

    IW_STL_Hash_Set* s = _formula_hash[matoms];

    if (s->end() == s->find(mformula)) {
      s->insert(mformula);

      IW_STL_Hash_Map_int* h =
          smiles_hash[matoms];    // very important to update the smiles hash too
      (*h)[usmi] = 1;

      unique_molecules++;
      return 1;
    }
#endif

int
UniqueMoleculesImplementation::IsUniqueInner(Molecule& m) {
  Preprocess(m);

  IWString usmi;  // The key that will be used for hash lookup.
  FormUniqueSmiles(m, usmi);

  if (_compare_as_graph) {
    const IWString formula = m.molecular_formula();

    Molecule g(m);
    g.change_to_graph_form();

    usmi = g.unique_smiles();

    usmi << ':' << formula;
  } else {
    usmi = m.unique_smiles();
  }

  //cerr << "Testing unique smiles '" << usmi << "'\n";

  if (_only_same_if_structure_and_name_the_same) {
    usmi << m.name();
  }

  const int matoms = m.natoms();

  while (matoms >= _smiles_hash.number_elements()) {
    _smiles_hash.add(new absl::flat_hash_set<IWString>);
    if (_use_atom_hash) {
      _atom_hash.add(new absl::flat_hash_set<uint64_t>());
    }
  }

  //cerr << matoms << " smiles_hash " << smiles_hash.number_elements() << " formula_hash " << formula_hash.number_elements() << '\n';

  if (_use_atom_hash) {
    const uint64_t h = m.quick_bond_hash();

    absl::flat_hash_set<uint64_t>* ha = _atom_hash[matoms];

    const auto f = ha->find(h);

    // Never seen before, molecule is unique.
    if (f == ha->end()) {
      ha->insert(h);

      absl::flat_hash_set<IWString>* h = _smiles_hash[matoms];
      h->emplace(usmi);

      ++_unique_molecules;
      return 1;
    }
  }

#ifdef FORMULA_CHECK_SLOWS_THINGS_DOWN
  if (_perform_formula_check) {
    PerformFormulaCheck(usmi);  // not implemented.
  }
#endif

  absl::flat_hash_set<IWString>* h = _smiles_hash[matoms];

  auto f = h->find(usmi);

  if (f == h->end()) {
    h->insert(usmi);
    return 1;
  } else {
    return 0;
  }
}

// Driver class for command line unique molecules determination.
class UniqueMoleculesOptions {
  private:
    int _verbose;

    int _molecules_read;

    Molecule_Output_Object _unique_molecule_stream;

    Molecule_Output_Object _duplicate_molecule_stream;

    int _function_as_filter = 0;

    Report_Progress _report_progress;

    // Tags for TDT input.
    IWString _smiles_tag;
    IWString _identifier_tag;

    // How should element transformations and removals be used?
    // By default, we apply the removal and transformations and the resulting molecule
    // is what gets written. But, what we probably want is for the transformed
    // molecule to be used for comparisons only, and we want the original
    // molecule written
    int _discard_molecule_changes;

  // Private functions
    int BuildPreviousMolecules(data_source_and_type<Molecule>& input,
                       UniqueMoleculesImplementation& um,
                       int& molecules);
  public:
    UniqueMoleculesOptions();

    int Initialise(Command_Line& cl);

    int molecules_read() const {
      return _molecules_read;
    }

    // We can pre-load the `um` caches with other sets of molecules.
    int BuildPreviousMolecules(const const_IWSubstring& fname, FileType input_type,
                       UniqueMoleculesImplementation& um, int& molecules);

    int Process(data_source_and_type<Molecule>& input, UniqueMoleculesImplementation& um);
    int Process(const char* fname, UniqueMoleculesImplementation& um, std::ostream& output);
    int Process(const char* fname,
                UniqueMoleculesImplementation& um,
                FileType input_type);
    int Process(iwstring_data_source& input,
                UniqueMoleculesImplementation& um,
                std::ostream& output);
    int IsUnique(Molecule& m, UniqueMoleculesImplementation& um);
};

UniqueMoleculesOptions::UniqueMoleculesOptions() {
  _verbose = 0;
  _molecules_read = 0;
  _function_as_filter = 0;
  _discard_molecule_changes = 0;
  _smiles_tag = "$SMI<";
  _identifier_tag = "PCN<";
}

int
UniqueMoleculesOptions::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('t')) {
    _discard_molecule_changes = 1;
    if (_verbose) {
      cerr << "Will discard molecule changes\n";
    }
  }

  if (cl.option_present('r')) {
    if (! _report_progress.initialise(cl, 'r', _verbose)) {
      cerr << "The report every option (-r) must be followed by a whole positive number\n";
      return 0;
    }
  }

  if (cl.option_present('G') && ! cl.option_present('f')) {
    cerr << "The identifier tag option (-G) only makes sense with the -f option\n";
    return 0;
  }

  if (cl.option_present('G')) {
    cl.value('G', _identifier_tag);
    if (! _identifier_tag.ends_with('<')) {
      _identifier_tag += '<';
    }

    if (_verbose) {
      cerr << "Molecules identified by '" << _identifier_tag << "' tag\n";
    }
  }

  if (cl.option_present('S')) {
    if (_function_as_filter) {
      cerr << "The -f and -S options are incompatible\n";
      return 0;
    }

    if (! cl.option_present('o')) {
      _unique_molecule_stream.add_output_type(FILE_TYPE_SMI);
    } else if (! _unique_molecule_stream.determine_output_types(cl)) {
      cerr << "Cannot discern output types for unique stream\n";
      return 0;
    }

    const_IWSubstring tmp = cl.string_value('S');

    if (_unique_molecule_stream.would_overwrite_input_files(cl, tmp)) {
      cerr << "Cannot overwrite input file(s) '" << tmp << "'\n";
      return 0;
    }

    if (! _unique_molecule_stream.new_stem(tmp, 1))  {   // causes files to be opened
      cerr << "Could not use stem '" << tmp << "' for ouput\n";
      return 4;
    }

    if (_verbose) {
      cerr << "Unique molecules written to stem '" << tmp << "'\n";
    }
  }

  if (cl.option_present('D')) {

    if (! cl.option_present('o')) {
      _duplicate_molecule_stream.add_output_type(FILE_TYPE_SMI);
    } else if (! _duplicate_molecule_stream.determine_output_types(cl)) {
      cerr << "Cannot discern output types for duplicate stream\n";
      return 0;
    }

    const_IWSubstring tmp = cl.string_value('D');

    if (_duplicate_molecule_stream.would_overwrite_input_files(cl, tmp)) {
      cerr << "Cannot overwrite input file(s) '" << tmp << "'\n";
      return 0;
    }

    if (! _duplicate_molecule_stream.new_stem(tmp, 1))  {   // causes files to be opened
      cerr << "Could not use stem '" << tmp << "' for duplicates\n";
      return 4;
    }

    if (_verbose) {
      cerr << "Duplicate written to stem '" << tmp << "'\n";
    }
  }

  if (! cl.option_present('S') && ! cl.option_present('D')) {
    _verbose = 1;
  }

  return 1;
}

int
UniqueMoleculesOptions::IsUnique(Molecule& m,
                        UniqueMoleculesImplementation& um) {
  int is_unique;
  if (_discard_molecule_changes) {
    Molecule mcopy(m);
    is_unique = um.IsUnique(mcopy);
  } else {
    is_unique = um.IsUnique(m);
  }

  if (is_unique) {
    if (_unique_molecule_stream.active()) {
      m.invalidate_smiles();
      _unique_molecule_stream.write(m);
    }
  } else {
    if (_duplicate_molecule_stream.active()) {
      m.invalidate_smiles();
      _duplicate_molecule_stream.write(m);
    }
  }

  return 1;
}

int
UniqueMoleculesOptions::Process(data_source_and_type<Molecule>& input,
                        UniqueMoleculesImplementation& um) {
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    ++_molecules_read;

    if (_report_progress()) {
      cerr << " processed " << _molecules_read << " molecules\n";
    }

    IsUnique(*m, um);
  }

  return 1;
}

int
UniqueMoleculesOptions::Process(iwstring_data_source& input,
                UniqueMoleculesImplementation& um,
                std::ostream& output)
{
  IW_TDT tdt;
  while (tdt.next(input) && output.good()) {
    const_IWSubstring smi;
    if (! tdt.dataitem_value(_smiles_tag, smi)) {
      cerr << "Yipes, cannot extract smiles from tdt\n";
      cerr << tdt;
      return 0;
    }

    Molecule m;
    if (! m.build_from_smiles(smi)) {
      cerr << "Very bad news, cannot parse smiles '" << smi << "'\n";
      cerr << tdt;
      return 0;
    }

    if (um.IsUnique(m)) {
      output << tdt;
      continue;
    }

    //  If verbose we need to report the ID of the dup.

    if (_verbose || _duplicate_molecule_stream.active()) {
      IWString id;
      tdt.dataitem_value(_identifier_tag, id);

      if (_verbose > 1) {
        cerr << "Is duplicate '" << id << "'\n";
      }

      if (_duplicate_molecule_stream.active()) {
        m.set_name(id);
        m.invalidate_smiles();
        _duplicate_molecule_stream.write(m);
      }
    }
  }

  return output.good();
}

int
UniqueMoleculesOptions::Process(const char* fname,
                UniqueMoleculesImplementation& um, std::ostream& output) {
  iwstring_data_source input(fname);
  if (! input.ok()) {
    cerr << "Cannot open filter input '" << fname << "'\n";
    return 0;
  }

  return Process(input, um, output);
}

int
UniqueMoleculesOptions::Process(const char* fname,
                UniqueMoleculesImplementation& um,
                FileType input_type) {
  if (FILE_TYPE_INVALID == input_type) {
    input_type = discern_file_type_from_name(fname);
    assert(FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok()) {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 0;
  }

  if (_verbose > 1) {
    input.set_verbose(_verbose);
  }

  return Process(input, um);
}

int
UniqueMoleculesOptions::BuildPreviousMolecules(data_source_and_type<Molecule>& input,
                       UniqueMoleculesImplementation& um,
                       int& molecules)
{
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    ++molecules;

    if (_report_progress()) {
      cerr << " processed " << molecules << " previously known molecules\n";
    }

    um.IngestPreviousMolecule(*m);
  }

  return molecules;
}

int
UniqueMoleculesOptions::BuildPreviousMolecules(const const_IWSubstring& fname, FileType input_type,
                       UniqueMoleculesImplementation& um, int& molecules)
{
  if (FILE_TYPE_INVALID == input_type) {
    input_type = discern_file_type_from_name(fname);
    assert(FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok()) {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 0;
  }

  return BuildPreviousMolecules(input, um, molecules);
}

static void
usage(int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
  cerr << "Filters out duplicate structures, based on unique smiles\n";
  cerr << "Usage: " << prog_name << " <options> <file1> <file2> ...\n";
  cerr << "  -l             strip to largest fragment\n";
  cerr << "  -a             compare as tautomers - skeleton and Hcount\n";
  cerr << "  -c             exclude chiral info - optical isomers will be duplicates\n";
  cerr << "  -z             exclude cis/trans bonding information\n";
  cerr << "  -I             ignore isotopic labels\n";
  cerr << "  -f             function as filter (TDT input)\n";
  cerr << "  -G <tag>       identifier tag when working as a filter\n";
  cerr << "  -p <fname>     specify previously collected molecules\n";
  cerr << "  -S <name>      specify output file name stem\n";
  cerr << "  -D <name>      write duplicate structures to <name>\n";
  cerr << "  -R <rxn>       perform reaction(s) on molecules before comparing\n";
  cerr << "  -r <number>    report progress every <number> molecules\n";
  cerr << "  -e             report all molecules together with counts\n";
  cerr << "  -j             items are the same only if both structure and name match\n";
  cerr << "  -T E1=E2       element transformations, enter '-t help' for details\n";
  cerr << "  -i <type>      specify input type\n";
  cerr << "  -o <type>      specify output type(s)\n";
  display_standard_aromaticity_options(cerr);
  cerr << "  -K ...         standard smiles options, enter '-K help' for info\n";
  display_standard_chemical_standardisation_options(cerr, 'g');
  cerr << "  -v             verbose output\n";

  exit(rc);
}

static int
unique_molecule(int argc, char** argv)
{
  Command_Line cl(argc, argv, "T:tag:D:vS:A:E:X:i:o:lczfG:p:Ir:n:K:R:jh");

  const int verbose = cl.option_count('v');

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  if (cl.option_present('E')) {
    if (! process_elements(cl, verbose, 'E')) {
      cerr << "Cannot discern elements, -E\n";
      usage(8);
    }
  }

  if (! process_standard_aromaticity_options(cl, verbose, 'A')) {
    cerr << "Cannot process standard aromaticity options\n";
    usage(2);
  }

  if (cl.option_present('K')) {
    if (! process_standard_smiles_options(cl, verbose, 'K')) {
      cerr << "Cannot initialise smiles options\n";
      return 5;
    }
  }

  UniqueMoleculesImplementation um;
  if (! um.Initialise(cl)) {
    cerr << "Cannot initialise unique molecules\n";
    usage(1);
  }

  int function_as_filter = 0;

  if (cl.option_present('f')) {
    if (cl.option_present('i')) {
      cerr << "The -i and -f options are incompatible\n";
      usage(1);
    }

    function_as_filter = 1;
    if (verbose) {
      cerr << "Will function as a TDT filter\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (function_as_filter) {    // don't check anything about the input type
  } else if (! cl.option_present('i')) {
    if (1 == cl.number_elements() && 0 == strncmp(cl[0], "-", 1))
      input_type = FILE_TYPE_SMI;
    else if (! all_files_recognised_by_suffix(cl)) {
      cerr << "Cannot automatically determine input type(s)\n";
      usage(8);
    }
  } else if (! process_input_type(cl, input_type)) {
    cerr << "Cannot determine input type\n";
    usage(1);
  }

  // Test this before opening any files

  if (cl.empty()) {
    cerr << "No input files specified\n";
    usage(1);
  }

  UniqueMoleculesOptions options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise unique molecules handline\n";
    return 1;
  }

  if (cl.option_present('p')) {
    int molecules = 0;

    const_IWSubstring p;
    for (int i = 0; cl.value('p', p, i); ++i) {
      if (! options.BuildPreviousMolecules(p, input_type, um, molecules)) {
        cerr << "Cannot process the -p option, '" << p << "'\n";
        return 72;
      }
    }

    if (verbose) {
      cerr << "read " << molecules << " molecules from previous set(s)\n";
    }
  }

  for (const char * fname : cl) {
    int rc;
    if (function_as_filter) {
      rc = options.Process(fname, um, std::cout);
    } else {
      rc = options.Process(fname, um, input_type);
    }
    if (rc == 0) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    cerr << "Read " << options.molecules_read() << " molecules\n";
    um.Report(cerr);
  }

  return 0;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = unique_molecule(argc, argv);

  return rc;
}
