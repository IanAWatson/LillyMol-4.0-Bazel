// Tester for unique smiles

#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/is_actually_chiral.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

using std::cerr;

namespace tsmiles {

int verbose = 0;

int permutations = 10;

int molecules_read = 0;

int molecules_containing_unmarked_chirality = 0;

int failures = 0;

int display_error_messages = 1;

int reduce_to_largest_fragment = 0;

int remove_chirality = 0;

int consider_chirality_in_unique_smiles = 1;

int try_explicit_hydrogen_variant = 0;

Molecule_Output_Object stream_for_failures;

Chemical_Standardisation chemical_standardisation;

Report_Progress report_progress;

void
usage(int rc) {
  cerr << "Tester for unique smiles generation\n";
  cerr << "  -p <n>         number of permutations of each molecule to generate\n";
  cerr << "  -c             exclude chirality\n";
  cerr << "  -F <fname>     write failed molecules to <fname>\n";
  cerr << "  -q             do NOT display failure messages\n";
  cerr << "  -l             reduce to largest fragment\n";
  cerr << "  -c             remove chirality from incoming molecules\n";
  cerr << "  -C             exclude chirality from the unique smiles determination\n";
  cerr << "  -h             after testing each molecules, make hydrogens explicit and retest\n";
  cerr << "  -b             use legacy uniqueness determination\n";
  cerr << "  -v             verbose output\n";
  exit(rc);
}

void
preprocess(Molecule& m) {
  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (remove_chirality) {
    m.remove_all_chiral_centres();
  }
}

void
make_report(std::ostream & output) {
  output << "Read " << molecules_read << " molecules, " << failures << " failures\n";
  if (molecules_containing_unmarked_chirality > 0)
    output << "Skipped " << molecules_containing_unmarked_chirality << " molecules containing unmarked chirality\n";
}

int
contains_unmarked_or_wrong_chirality(Molecule& m) {
  const int nchiral = m.chiral_centres();
  if (nchiral == 0) {
    return 0;
  }

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (m.chiral_centre_at_atom(i) != NULL) {
      if (! is_actually_chiral(m, i)) {
        m.remove_chiral_centre_at_atom(i);  // Could get strange, hopefully OK.
      }
      continue;
    }
    const Atom * a = m.atomi(i);
    if (a->atomic_number() != 6) {
      continue;
    }
    if (a->ncon() == a->nbonds()) {
      continue;
    }
    if (m.ring_bond_count(i) > 0) {   // Maybe...
      continue;
    }
    if (a->ncon() <= 2) {
      continue;
    }
    if (is_actually_chiral(m, i)) {
      return 1;
    }
  }

  return 0;
}

int
tsmiles(Molecule& m) {
  if (contains_unmarked_or_wrong_chirality(m)) {
    molecules_containing_unmarked_chirality++;
    return 1;
  }
  const IWString initial_smiles = m.unique_smiles();
  for (int i = 0; i < permutations; ++i) {
    Molecule mcopy;
    const IWString& rsmi = m.random_smiles();
    if (! mcopy.build_from_smiles(rsmi)) {
      cerr << "Cannot parse smiles '" << rsmi << "' from " << m.name() << '\n';
      return 0;
    }
    if (verbose > 2) {
      cerr << "tsmiles i = " << i << ' ';
      write_atom_map_number_labelled_smiles(mcopy, false, cerr) << '\n';
    }
    if (mcopy.unique_smiles() == initial_smiles) {
      continue;
    }
    failures++;
    if (stream_for_failures.active())
      stream_for_failures.write(m);
    if (display_error_messages) {
      cerr << "SMiles mismatch: got " << mcopy.unique_smiles() << '\n';
      cerr << "expected             " << initial_smiles << ' ' << m.name() << '\n';
      write_atom_map_number_labelled_smiles(mcopy, false, cerr) << '\n';
      return 1;
    }
  }
  return 1;
}

int
tsmiles(data_source_and_type<Molecule>& input) {
  Molecule * m;
  while ((m = input.next_molecule()) != NULL) {
    std::unique_ptr<Molecule> free_m(m);
    molecules_read++;
    preprocess(*m);
    tsmiles(*m);
    if (try_explicit_hydrogen_variant) {
      m->make_implicit_hydrogens_explicit();
      tsmiles(*m);
    }
    if (report_progress()) {
      make_report(cerr);
    }
  }

  return 1;
}

int
tsmiles(const char * fname,
        FileType input_type) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert (FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok())
  {
    cerr << "Cannot open '" << fname << "' for input\n";
    return 1;
  }

  return tsmiles(input);
}

int tsmiles(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:l:g:K:i:r:p:qF:hcCs:b");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (! process_elements(cl, verbose))
  {
    cerr << "Cannot parse element specifications\n";
    usage(2);
  }

  if (cl.option_present('K'))
  {
    if (! process_standard_smiles_options(cl, verbose, 'K'))
      usage(5);
  }

  if (! process_standard_aromaticity_options(cl, verbose))
    usage(6);

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type))
      usage(3);
  } else if (! all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot determine input file type\n";
    return 1;
  }

  if (cl.option_present('g'))
  {
    if(! chemical_standardisation.construct_from_command_line(cl, verbose > 2))
    {
      cerr << "Cannot parse -g option\n";
      return 61;
    }
  }

  if (cl.option_present('p')) {
    if (! cl.value('p', permutations) || permutations <= 0) {
      cerr << "The number of permutations (-p) must be a whole +ve number\n";
      usage(1);
    }
    if (verbose)
      cerr << "Will perform " << permutations << " permutations on each molecule\n";
  }

  if (cl.option_present('s')) {
    uint32_t seed;
    if (! cl.value('s', seed)) {
      cerr << "The -s option needs a valid int\n";
      usage(1);
    }
    set_smiles_random_number_seed(seed);
    if (verbose)
      cerr << "Using seed " << seed << '\n';
  }

  if (cl.option_present('r')) {
    if (! report_progress.initialise(cl, 'r', verbose)) {
      cerr << "Cannot initialise progress reporting\n";
      return 1;
    }
  }

  if (cl.option_present('b')) {
    set_unique_smiles_legacy_atom_ordering(1);
    if (verbose)
      cerr << "Will use legacy atom ordering in unique determination\n";
  }

  if (cl.option_present('q')) {
    display_error_messages = 0;

    if (verbose)
      cerr << "Will NOT display error messages\n";
  }

  if (cl.option_present('c')) {
    remove_chirality = 1;
    if (verbose)
      cerr << "Chirality removed from molecules\n";
  }

  if (cl.option_present('C')) {
    consider_chirality_in_unique_smiles = 0;
    if (verbose)
      cerr << "Chirality not included in unique smiles determination\n";
    set_include_chiral_info_in_smiles(0);
  }

  if (cl.option_present('h')) {
    try_explicit_hydrogen_variant = 1;
    if (verbose)
      cerr << "Will test an explicit Hydrogen variant of each input molecule\n";
  }

  if (cl.empty()) {
    cerr << "Must specify one or more input files\n";
    usage(1);
  }

  if (cl.option_present('F')) {
    const IWString f = cl.option_value('F');
    stream_for_failures.add_output_type(FILE_TYPE_SMI);
    if (stream_for_failures.would_overwrite_input_files(cl, f)) {
      cerr << "Cannot overwrite input files '" << f << "'\n";
      return 1;
    }
    if (! stream_for_failures.new_stem(f)) {
      cerr << "Cannot open stream for failed molecules '" << f << "'\n";
      return 1;
    }

    if (verbose)
      cerr << "Failed molecules written to '" << f << "'\n";
  }

  for (const char * fname : cl) {
    if (! tsmiles(fname, input_type)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  make_report(cerr);

  return 0;
}

}  // namespace tsmiles

int
main (int argc, char ** argv)
{
  int rc = tsmiles::tsmiles(argc, argv);

  return rc;
}
