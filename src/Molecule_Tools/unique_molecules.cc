/*
  Sometimes we need to determine the unique molecules from a file
  of structures.
*/

#include <stdlib.h>
#include <sys/stat.h>

#include <iostream>
#include <memory>

#define RESIZABLE_ARRAY_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/report_progress.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"

#include "Molecule_Tools/unique_molecules_api.h"

using std::cerr;

const char* prog_name = nullptr;

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
                       unique_molecules::UniqueMoleculesImplementation& um,
                       int& molecules);
  public:
    UniqueMoleculesOptions();

    int Initialise(Command_Line& cl);

    int molecules_read() const {
      return _molecules_read;
    }

    // We can pre-load the `um` caches with other sets of molecules.
    int BuildPreviousMolecules(const const_IWSubstring& fname, FileType input_type,
                       unique_molecules::UniqueMoleculesImplementation& um, int& molecules);

    int Process(data_source_and_type<Molecule>& input, unique_molecules::UniqueMoleculesImplementation& um);
    int Process(const char* fname, unique_molecules::UniqueMoleculesImplementation& um, std::ostream& output);
    int Process(const char* fname,
                unique_molecules::UniqueMoleculesImplementation& um,
                FileType input_type);
    int Process(iwstring_data_source& input,
                unique_molecules::UniqueMoleculesImplementation& um,
                std::ostream& output);
    int IsUnique(Molecule& m, unique_molecules::UniqueMoleculesImplementation& um);
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
                        unique_molecules::UniqueMoleculesImplementation& um) {
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
                        unique_molecules::UniqueMoleculesImplementation& um) {
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    ++_molecules_read;

    if (_report_progress()) {
      cerr << " processed " << _molecules_read << " molecules " << um.unique_molecules() << " unique\n";
    }

    IsUnique(*m, um);
  }

  return 1;
}

int
UniqueMoleculesOptions::Process(iwstring_data_source& input,
                unique_molecules::UniqueMoleculesImplementation& um,
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
                unique_molecules::UniqueMoleculesImplementation& um, std::ostream& output) {
  iwstring_data_source input(fname);
  if (! input.ok()) {
    cerr << "Cannot open filter input '" << fname << "'\n";
    return 0;
  }

  return Process(input, um, output);
}

int
UniqueMoleculesOptions::Process(const char* fname,
                unique_molecules::UniqueMoleculesImplementation& um,
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

// Read molecules from `input` and call um.IngestPreviousMolecule
// on each such molecule.
int
UniqueMoleculesOptions::BuildPreviousMolecules(data_source_and_type<Molecule>& input,
                       unique_molecules::UniqueMoleculesImplementation& um,
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
                       unique_molecules::UniqueMoleculesImplementation& um, int& molecules)
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

  unique_molecules::UniqueMoleculesImplementation um;
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
