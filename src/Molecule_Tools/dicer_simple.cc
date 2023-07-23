// Dicer via dicer_api

#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/dicer_api.h"

namespace dicer_api {

using std::cerr;

void
Usage(int rc) {
  cerr << "Dicer variant using dicer_api\n";

  cerr << " -c                remove chirality\n";
  cerr << " -l                strip to largest fragment\n";
  cerr << " -v                verbose output\n";

  ::exit(rc);
}

class LocalOptions {
  private:
    int _verbose;

    int _molecules_read;

    int _reduce_to_largest_fragment;

    int _remove_chirality;

    Chemical_Standardisation _chemical_standardisation;

  public:
    LocalOptions();

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);

    int Report(std::ostream& output) const;
};

LocalOptions::LocalOptions() {
  _verbose = 0;
  _molecules_read = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 1;
}

int
LocalOptions::Initialise(Command_Line& cl) {
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

  return 1;
}

int
LocalOptions::Preprocess(Molecule& m) {
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
LocalOptions::Report(std::ostream& output) const {
  output << "LocalOptions:read " << _molecules_read << " molecules\n";
  if (_molecules_read == 0) {
    return 1;
  }

  return 1;
}

int
Dicer(DicerApi& dicer,
      Molecule& m,
      IWString_and_File_Descriptor& output) {
  std::cout << m.smiles() << ' ' << m.name() << '\n';
  Sparse_Fingerprint_Creator sfc;
  return dicer.Process(m, sfc);
}

int
Dicer(LocalOptions& local_options,
      DicerApi& dicer,
      data_source_and_type<Molecule>& input,
      IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! local_options.Preprocess(*m)) {
      return 0;
    }

    Dicer(dicer, *m, output);
  }

  return 1;
}

int
Dicer(LocalOptions& local_options, DicerApi& dicer,
            const char * fname,
            FileType input_type,
            IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return Dicer(local_options, dicer, input, output);
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

  LocalOptions local_options;
  if (! local_options.Initialise(cl)) {
    cerr << "Cannot initialise local options\n";
    Usage(1);
  }

  DicerApi dicer;
  if (! dicer.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (all_files_recognised_by_suffix(cl)) {
  } else {
    input_type = FILE_TYPE_SMI;
  }

  IWString_and_File_Descriptor output(1);
  for (const char* fname : cl) {
    if (! Dicer(local_options, dicer, fname, input_type, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    dicer.Report(cerr);
  }

  return 0;
}

}  // namespace dicer_api

int
main(int argc, char** argv) {
  int rc = dicer_api::Main(argc, argv);

  return rc;
}
