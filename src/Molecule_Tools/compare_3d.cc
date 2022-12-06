// Quick throwaway tool, not really a robust application.
// Input file consists of consecutive duplicate molecules. Compare
// their 3D structures. We assume that the atom orders are the same.

// Start with two 3D smiles files, each containing unique smiles.
//    fileconv -o usmi -o smi3d file1.sdf
//    fileconv -o usmi -o smi3d file2.sdf
// join them with fetch_smiles_quick
//   fetch_smiles_quick -v -c 2 -C 2 file1.smi file2.smi > both.smi
// That output has smiles1 id smiles2. move smiles2 to a separate line
//   sed -e 's/ /_/' -e 's/ /\n/' -e 's/_/ /' both.smi > bb.smi
// bb.smi is the input to this tool.

#include <algorithm>
#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

namespace compare_3d {

using std::cerr;

extern "C" void u3b_(const double * w, double * c1, double * c2, const int * n, const int * mode, double *rms, double * u, double * t, int  * ier);

void
Usage(int rc) {

  cerr << "Compares the 3D coordinates of adjacent molecules\n";
  cerr << "Consecutive molecules must be the same molecule, just with different coordinates\n";
  cerr << " -S <fname>        write superimposed molecules to <fname>\n";
  cerr << " -c                remove chirality\n";
  cerr << " -l                strip to largest fragment\n";
  cerr << " -v                verbose output\n";

  ::exit(rc);
}

class CoordinateComparison {
  private:
    int _verbose;

    int _molecules_read;

    int _reduce_to_largest_fragment;

    int _remove_chirality;

    Accumulator<double> _acc_initial_rms;
    Accumulator<double> _acc_final_rms;

    Molecule_Output_Object _aligned;

    FileType _input_type;

    Chemical_Standardisation _chemical_standardisation;

  public:
    CoordinateComparison();

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);

    int Process(Molecule& m1, Molecule& m2, IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;

    FileType input_type() const {
      return _input_type;
    }
};

CoordinateComparison::CoordinateComparison() {
  _verbose = 0;
  _molecules_read = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _input_type = FILE_TYPE_INVALID;
}

int
CoordinateComparison::Initialise(Command_Line& cl) {
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

  if (cl.option_present('S')) {
    IWString fname = cl.string_value('S');
    _aligned.add_output_type(FILE_TYPE_SDF);
    _aligned.add_output_type(FILE_TYPE_SMI);
    set_append_coordinates_after_each_atom(1);
    if ( !_aligned.new_stem(fname)) {
      cerr << "Cannot set output file name '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Aligned molecules written to " << fname << "'\n";
    }
  }

  return 1;
}

int
CoordinateComparison::Preprocess(Molecule& m) {
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
CoordinateComparison::Process(Molecule& m1, Molecule& m2,
                         IWString_and_File_Descriptor& output) {
  ++_molecules_read;
  
  const int matoms = m1.natoms();
  m2.set_name(m1.name());

  double * c1 = new double[matoms * 3 * 2]; std::unique_ptr<double[]> free_c1(c1);
  double * c2 = c1 + (matoms * 3);

  for (int i = 0; i < matoms; ++i) {
    c1[i * 3] = m1.x(i);
    c1[i * 3 + 1] = m1.y(i);
    c1[i * 3 + 2] = m1.z(i);
    c2[i * 3] = m2.x(i);
    c2[i * 3 + 1] = m2.y(i);
    c2[i * 3 + 2] = m2.z(i);
  }

  double* weight = new double[matoms]; std::unique_ptr<double[]> free_weight(weight);
  std::fill_n(weight, matoms, 1.0 / matoms);

  int mode = 1;
  double u[9];
  double rms;
  double t[3];
  int ier = 0;

  u3b_(weight, c1, c2, &matoms, &mode, &rms, u, t, &ier);

  if (ier == 0)
    ;
  else if (ier == -1)
    cerr << "Reaction_Rotate_Fragment::process:superposition not unique, but optimal\n";
  else
  {
    cerr << "Reaction_Rotate_Fragment::process:u3b failed, ier " << ier << '\n';
    return 0;
  }
  cerr << "RMS " << rms << ' ' << m1.name() << '\n';
  _acc_initial_rms.extra(rms);

  double rotmat11 = u[0];
  double rotmat12 = u[1];
  double rotmat13 = u[2];
  double rotmat21 = u[3];
  double rotmat22 = u[4];
  double rotmat23 = u[5];
  double rotmat31 = u[6];
  double rotmat32 = u[7];
  double rotmat33 = u[8];

  for (int i = 0; i < matoms; i++)
  {
#ifdef DEBUG_PROCESS_3D_REPLACE
    cerr << "Atom " << i << " '" << m.smarts_equivalent_for_atom(i) << "' is moving\n";
#endif

    const Atom * a = m2.atomi(i);

    double x0 = a->x() - t[0];
    double y0 = a->y() - t[1];
    double z0 = a->z() - t[2];

    double xx = rotmat11 * x0 + rotmat12 * y0 + rotmat13 * z0;
    double yy = rotmat21 * x0 + rotmat22 * y0 + rotmat23 * z0;
    double zz = rotmat31 * x0 + rotmat32 * y0 + rotmat33 * z0;

    m2.setxyz( i, static_cast<coord_t> (xx), static_cast<coord_t> (yy), static_cast<coord_t> (zz) );
  }

  double tot = 0.0;
  for (int i = 0; i < matoms; ++i) {
    const Space_Vector<float> a1 = m1.atom(i);
    const Space_Vector<float> a2 = m2.atom(i);
    double d = a1.distance(a2);
    tot += d * d;
  }

  double final_rms = sqrt(tot/ matoms);
  cerr << "Final rms " << final_rms << " diff " << abs(final_rms - rms) << ' ' << m1.name() << '\n';
  _acc_final_rms.extra(final_rms);

  if (_aligned.active()) {
    _aligned.write(m1);
    _aligned.write(m2);
    Molecule both(m1);
    both.add_molecule(&m2);
    both.set_name(m1.name());
    _aligned.write(both);
  }

  return 1;
}

int
CoordinateComparison::Report(std::ostream& output) const {
  output << "CoordinateComparison:read " << _molecules_read << " molecules\n";
  if (_molecules_read == 0) {
    return 1;
  }

  output << "Initial RMS btw " << _acc_initial_rms.minval() << " and " << _acc_initial_rms.maxval() << " ave " << _acc_initial_rms.average() << '\n';
  output << "Final   RMS btw " << _acc_final_rms.minval() << " and " << _acc_final_rms.maxval() << " ave " << _acc_final_rms.average() << '\n';

  return 1;
}

int
CompareCoordinates(CoordinateComparison& coordinate_comparison,
            Molecule& m1, Molecule& m2,
            IWString_and_File_Descriptor& output) {
  const int matoms = m1.natoms();
  if (m2.natoms() != matoms) {
    cerr << "Atom count mismatch " << m1.name() << " has " << matoms << " other has " << m2.natoms() << '\n';
    return 0;
  }

  for (int i = 0; i < matoms; ++i) {
    const IWString smt1 = m1.smarts_equivalent_for_atom(i);
    const IWString smt2 = m2.smarts_equivalent_for_atom(i);
    if (smt1 == smt2) {
      continue;
    }
    cerr << "Smarts mismatch, atom " << i << " in " << m1.name() << " " << smt1 << " vs " << smt2 << '\n';
    return 0;
  }

  return coordinate_comparison.Process(m1, m2, output);
}

int
CompareCoordinates(CoordinateComparison& coordinate_comparison,
            data_source_and_type<Molecule>& input,
            IWString_and_File_Descriptor& output) {
  Molecule * m1;
  while ((m1 = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m1(m1);
    Molecule* m2 = input.next_molecule();
    std::unique_ptr<Molecule> free_m2(m2);

    if (! CompareCoordinates(coordinate_comparison, *m1, *m2, output)) {
      return 0;
    }
  }

  return 1;
}

int
CompareCoordinates(CoordinateComparison& coordinate_comparison,
            const char * fname,
            IWString_and_File_Descriptor& output) {
  FileType input_type = coordinate_comparison.input_type();
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return CompareCoordinates(coordinate_comparison, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:i:g:lcS:");
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

  CoordinateComparison coordinate_comparison;
  if (! coordinate_comparison.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);
  for (const char* fname : cl) {
    if (! CompareCoordinates(coordinate_comparison, fname, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    coordinate_comparison.Report(cerr);
  }

  return 0;
}

}  // namespace compare_3d

int
main(int argc, char** argv) {
  int rc = compare_3d::Main(argc, argv);

  return rc;
}
