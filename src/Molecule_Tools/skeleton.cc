// Form all possible scaffolds from a molecule.

#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/atom_separations.h"

namespace separated_atoms {

using std::cerr;

void
Usage(int rc) {

  cerr << " -c                remove chirality\n";
  cerr << " -l                strip to largest fragment\n";
  cerr << " -v                verbose output\n";

  ::exit(rc);
}

class CoreReplacement {
  private:
    int _verbose;

    int _molecules_read;

    int _reduce_to_largest_fragment;

    int _remove_chirality;

    FileType _input_type;

    Chemical_Standardisation _chemical_standardisation;

  public:
    CoreReplacement();

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);

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

  if (m.empty()) {
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

  return 1;
}

int
ReplaceCore(CoreReplacement& core_replacement,
            Molecule& m,
            IWString_and_File_Descriptor& output) {
  return 1;
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
  Command_Line cl(argc, argv, "vE:A:i:g:lc");
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

}  // namespace separated_atoms

int
main(int argc, char** argv) {
  int rc = separated_atoms::Main(argc, argv);

  return rc;
}
