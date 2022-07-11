// Scan through SMU and extract the smiles and test BondTopology

#include <iostream>

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

#include "smu/dataset.pb.h"
#include "smu/support.h"

namespace get_smiles_bt {

using std::cerr;
using GoogleSmu::BondTopology;

void Usage(int rc) {
  ::exit(rc);
}

struct Options {
  int verbose = 0;

  int conformers_read = 0;

  int no_molecule = 0;

  int no_bond_topologies = 0;

  int no_openbabel_smiles = 0;

  int babel_rdkit_same = 0;
  int babel_lillymol_same = 0;
  int rdkit_lillymol_same = 0;

  Report_Progress report_progress;

  public:
    int Build(Command_Line& cl);

    int Report(std::ostream& output) const;
};

int
Options::Build(Command_Line& cl) {
  verbose = cl.option_count('v');

  if (cl.option_present('r')) {
    if (! report_progress.initialise(cl, 'r', verbose)) {
      cerr << "Cannot initialise progress reporting\n";
      return 0;
    }
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << conformers_read << " conformers\n";
  output << no_bond_topologies << " had no bond topologies\n";
  output << no_molecule << " had no molecule\n";
  output << no_openbabel_smiles << " had no openbabel smiles\n";
  output << babel_rdkit_same << " openbabel and rdkit the same\n";
  output << babel_lillymol_same << " openbabel and lillymol the same\n";
  output << rdkit_lillymol_same << " rdkit and lillymol the same\n";

  return 1;
}

int
GetSmilesBT(Options& options,
            GoogleSmu::Molecule& conformer,
            IWString_and_File_Descriptor& output) {
  if (conformer.bond_topologies().size() == 0) {
    options.no_bond_topologies++;
    return 1;
  }

  options.report_progress();  // For some reason is not working.

  if (conformer.properties().errors().fate() != GoogleSmu::Properties::FATE_SUCCESS) {
    return 1;
  }

  constexpr char kSep = ' ';

  output << conformer.molecule_id();
  const std::string& openbabel_smiles = conformer.properties().smiles_openbabel();
  if (openbabel_smiles.length() > 0) {
    output << kSep << openbabel_smiles;
  } else {
    ++options.no_openbabel_smiles;
    output << kSep << '*';
  }
  for (int i = 0; i < conformer.bond_topologies_size(); ++i) {
    const BondTopology& bt = conformer.bond_topologies(i);

    std::optional<Molecule> maybe_mol = smu::MoleculeFromBondTopology(bt);
    if (! maybe_mol) {
      ++options.no_molecule;
      continue;
    }
    maybe_mol->remove_all(1);
    output << kSep << i;
    output << kSep << bt.smiles();  // RDkit smiles
    output << kSep << maybe_mol->unique_smiles();  // LillyMol smiles
    output << kSep << bt.is_starting_topology();
    output << kSep << bt.topology_score();
    output << kSep << bt.geometry_score();

    if (bt.is_starting_topology()) {
      if (maybe_mol->unique_smiles() == bt.smiles()) {
        ++options.rdkit_lillymol_same;
      }
      if (! openbabel_smiles.empty()) {
        if (openbabel_smiles == bt.smiles()) {
          ++options.babel_rdkit_same;
        }
        if (maybe_mol->unique_smiles() == openbabel_smiles) {
          ++options.babel_lillymol_same;
        }
      }
    }
  }
  output << '\n';
  output.write_if_buffer_holds_more_than(8192);

  return 1;
}

int
GetSmilesBT2(Options& options,
             const const_IWSubstring& buffer,
             IWString_and_File_Descriptor& output) {
  const std::string as_string(buffer.data(), buffer.length());
  GoogleSmu::Molecule conformer;
  if (! conformer.ParseFromString(as_string)) {
    cerr << "Cannot decode proto\n";
    return 0;
  }
  options.conformers_read++;

  return GetSmilesBT(options, conformer, output);
}

int
GetSmilesBT(Options& options,
            iw_tf_data_record::TFDataReader& reader,
            IWString_and_File_Descriptor& output) {
  while (reader.good() && ! reader.eof()) {
    std::optional<const_IWSubstring> data = reader.Next();
    if (! data) {
      return 1;
    }
    if (! GetSmilesBT2(options, *data, output)) {
      return 0;
    }
  }

  return 1;
}

int
GetSmilesBT(Options& options,
            const char * fname,
            IWString_and_File_Descriptor& output) {
  iw_tf_data_record::TFDataReader reader(fname);
  if (! reader.good()) {
    cerr << "GetSmilesBT:cannto open '" << fname << "'\n";
    return 0;
  }

  return GetSmilesBT(options, reader, output);
}

int
GetSmilesBT(int argc, char** argv) {
  Command_Line cl(argc, argv, "v");

  Options options;
  if (! options.Build(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }
  IWString_and_File_Descriptor output(1);

  for (const char * fname: cl) {
    if (! GetSmilesBT(options, fname, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  options.Report(cerr);

  return 0;
}


}  // namespace get_smiles_bt

int
main(int argc, char ** argv)
{
  int rc = get_smiles_bt::GetSmilesBT(argc, argv);

  return rc;
}
