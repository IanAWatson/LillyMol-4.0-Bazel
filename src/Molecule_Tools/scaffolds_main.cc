#include <stdlib.h>

#include <iostream>

#include "google/protobuf/text_format.h"

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"

#include "Molecule_Tools/scaffolds.h"
#include "Molecule_Tools/scaffolds.pb.h"

namespace scaffolds_main {

using std::cerr;
using scaffolds::ScaffoldFinder;

class LocalOptions {
  private:
    int _verbose;

    int _molecules_read;

    int _reduce_to_largest_fragment;

    int _remove_chirality;

    int _non_ring_molecules_skipped;

    Chemical_Standardisation _chemical_standardisation;

    int _write_parent;

    // Output can be either as smiles or proto.
    int _write_smiles;

  // private functions
 
    int WriteSmiles(Molecule& parent,
                    const Scaffolds::ScaffoldData& result,
                    IWString_and_File_Descriptor& output) const;

  public:
    LocalOptions();

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);

    void AnotherNonRingMolecule() {
      ++_non_ring_molecules_skipped;
    }

    int Write(Molecule& parent,
              const Scaffolds::ScaffoldData& result,
              IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

LocalOptions::LocalOptions() {
  _molecules_read = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _non_ring_molecules_skipped = 0;
  _write_parent = 0;
  _write_smiles = 1;
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

  if (cl.option_present('p')) {
    _write_parent = 1;
    if (_verbose) {
      cerr << "Will write the parent molecule\n";
    }
  }

  if (cl.option_present('y')) {
    _write_smiles = 0;
    if (_verbose) {
      cerr << "Will write as proto form\n";
    }
  }

  return 1;
}

int
LocalOptions::Preprocess(Molecule& m) {
  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (m.empty()) {
    return 0;
  }

  return 1;
}

int
LocalOptions::Write(Molecule& parent,
                    const Scaffolds::ScaffoldData& result,
                    IWString_and_File_Descriptor& output) {
  static constexpr char kSep = ' ';

  if (_write_smiles) {
    return WriteSmiles(parent, result, output);
  }

  static google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);

  if (_write_parent) {
    output << parent.smiles() << kSep << parent.name() << kSep << "parent\n";
  }

  std::string buffer;

  printer.PrintToString(result, &buffer);
  output << buffer << '\n';
  output.write_if_buffer_holds_more_than(8192);

  return 1;

  for (const Scaffolds::ScaffoldSubset& subset : result.subset()) {
    printer.PrintToString(subset, &buffer);
    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

int
LocalOptions::WriteSmiles(Molecule& parent,
                    const Scaffolds::ScaffoldData& result,
                    IWString_and_File_Descriptor& output) const {
  static constexpr char kSep = ' ';
  if (_write_parent) {
    output << parent.smiles() << kSep << parent.name() << '\n';
  }

  int ndx = 0;
  for (const auto& subset : result.subset()) {
    output << subset.smi() << kSep << parent.name() << '.' << ndx << kSep << subset.ring_sys() << '\n';
    ++ndx;
    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

void
Usage(int rc) {
  cerr << " -G <fname>             config containing Scaffolds::ScaffoldsOptions text proto\n";
  cerr << " -y                     output as Scaffolds::ScaffoldData textproto\n";
  ::exit(rc);
}
  
int
Scaffolds(LocalOptions& local_options,
          ScaffoldFinder& make_scaffolds,
          Molecule& m,
          IWString_and_File_Descriptor& output) {
  if (m.nrings() == 0) {
    local_options.AnotherNonRingMolecule();
    return 1;
  }

  Scaffolds::ScaffoldData result;
  result.set_smi(m.smiles().AsString());
  result.set_par(m.name().AsString());

  if (! make_scaffolds.MakeScaffolds(m, result)) {
    return 0;
  }

  return local_options.Write(m, result, output);
}

int
Scaffolds(LocalOptions& local_options,
          ScaffoldFinder& make_scaffolds,
          data_source_and_type<Molecule>& input,
          IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! local_options.Preprocess(*m)) {
      return 0;
    }

    if (! Scaffolds(local_options, make_scaffolds, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
Scaffolds(LocalOptions& local_options,
          ScaffoldFinder& make_scaffolds,
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

  return Scaffolds(local_options, make_scaffolds, input, output);
}

int
Scaffolds(int argc, char** argv) {
  Command_Line cl(argc, argv, "vi:A:E:g:clpG:y");
  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    return 1;
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
    return 1;
  }

  scaffolds::ScaffoldFinder make_scaffolds;

  if (! make_scaffolds.Initialise(cl, 'G')) {
    cerr << "Cannot initialise calculation\n";
    return 1;
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(6);
    }
  }

  if (FILE_TYPE_INVALID != input_type) {  // great, explicitly specified
  } else if (1 == cl.number_elements() && 0 == strcmp("-", cl[0])) {  // reading from a pipe, assume smiles input
    if (verbose)
      cerr << "Assuming smiles input from pipe read\n";
    input_type = FILE_TYPE_SMI;
  } else if (all_files_recognised_by_suffix(cl)) {
  } else {
    cerr << "Cannot discern file types from names\n";
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  for (const char * fname: cl) {
    const const_IWSubstring as_string(fname);

    if (verbose) {
      cerr << "Processing '" << fname << "'\n";
    }

    if (! Scaffolds(local_options, make_scaffolds, fname, input_type, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  if (verbose) {
    make_scaffolds.Report(cerr);
  }

  return 0;
}

}  // namespace scaffolds_main

int
main(int argc, char **argv) {
  int rc = scaffolds_main::Scaffolds(argc, argv);

  return rc;
}
