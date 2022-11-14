// Discern topology from geometry

#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"

#include "Molecule_Tools/topology_from_geometry.h"

namespace topology_from_geometry {

using std::cerr;

void
Usage(int rc) {
  cerr << "Discerns plausible bond topologies from atomic positions\n";

  cerr << " -v                verbose output\n";

  ::exit(rc);
}

int
DoTopologyFromGeometry(TopologyFromGeometry& topology_from_geometry,
            Molecule& m,
            IWString_and_File_Descriptor& output) {
  if (! topology_from_geometry.Process(m)) {
    cerr << "Cannot process '" << m.name() << '\n';
    return 0;
  }

  return 1;
}

int
DoTopologyFromGeometry(TopologyFromGeometry& topology_from_geometry,
            data_source_and_type<Molecule>& input,
            IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    m->remove_all_bonds();

    if (! topology_from_geometry.Preprocess(*m)) {
      return 0;
    }

    if (! DoTopologyFromGeometry(topology_from_geometry, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
DoTopologyFromGeometry(TopologyFromGeometry& topology_from_geometry,
            const char * fname,
            IWString_and_File_Descriptor& output) {
  FileType input_type = topology_from_geometry.input_type();
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return DoTopologyFromGeometry(topology_from_geometry, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:i:B:");
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

  TopologyFromGeometry topology_from_geometry;
  if (! topology_from_geometry.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);
  for (const char* fname : cl) {
    if (! DoTopologyFromGeometry(topology_from_geometry, fname, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    topology_from_geometry.Report(cerr);
  }

  return 0;
}

}  // namespace topology_from_geometry

int
main(int argc, char** argv) {
  int rc = topology_from_geometry::Main(argc, argv);

  return rc;
}
