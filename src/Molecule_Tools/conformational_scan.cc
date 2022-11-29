// Scan around a given bond and try to avoid close atoms.

#include <iostream>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

namespace conformational_scan {

using std::cerr;

void
Usage(int rc) {
  cerr << "Perform a conformational scan in order to avoid atom clashes\n";
  cerr << " -T <dist>         distance threshold for acceptance\n";
  cerr << " -r <deg>          rotation angle\n";
  cerr << " -q <query>        query to define the bond to be scanned (first 2 matched atoms)\n";
  cerr << " -s <smarts>       query as smarts\n";
  cerr << " -S <fname>        file for successfully relaxed molecules\n";
  cerr << " -U <fname>        file for molecules not relaxed\n";
  cerr << " -c                remove chirality\n";
  cerr << " -l                strip to largest fragment\n";
  cerr << " -v                verbose output\n";

  ::exit(rc);
}

// The main function needs various optional behaviours.
class Options {
  private:
    int _ignore_molecules_not_matching_queries = 0;
    int _write_molecules_not_matching_queries = 0;
    int _molecules_not_matching_queries;

    FileType _input_type;

    Molecule_Output_Object _output;
    Molecule_Output_Object _failures;

  public:
    Options();

    int Initialise(Command_Line& cl);

    FileType input_type() const {
      return _input_type;
    }

    int FailuresActive() const {
      return _failures.active();
    }

    int WriteSuccess(Molecule& m) {
      return _output.write(m);
    }

    int WriteFailed(Molecule& m) {
      return _failures.write(m);
    }
};

Options::Options() {
  _ignore_molecules_not_matching_queries = 0;
  _write_molecules_not_matching_queries = 0;
  _molecules_not_matching_queries = 0;
  _input_type = FILE_TYPE_SDF;
}

int
Options::Initialise(Command_Line& cl) {
  const int verbose = cl.option_count('v');

  if (cl.option_present('z')) {
    const_IWSubstring z;
    for (int i = 0; cl.value('z', z, i); ++i) {
      if (z == 'i') {
        _ignore_molecules_not_matching_queries = 1;
      } else if (z == 'w') {
        _write_molecules_not_matching_queries = 1;
      } else {
        cerr << "Unrecognised -z qualifier '" << z << "'\n";
        return 0;
      }
    }
  }

  if (! cl.option_present('S')) {
    cerr << "Must specify output file name via the -S option\n";
    return 0;
  }

  if (cl.option_present('o')) {
    if (! _output.determine_output_types(cl, 'o')) {
      cerr << "Cannot determine output type(s) -o\n";
      return 1;
    }
    _failures.determine_output_types(cl, 'o');
  } else {
    _output.add_output_type(FILE_TYPE_SDF);
    _failures.add_output_type(FILE_TYPE_SDF);
  }

  if (cl.option_present('S')) {
    IWString fname = cl.string_value('S');
    if (_output.would_overwrite_input_files(cl, fname)) {
      cerr << "Cannot overwrite input file(s) -S '" << fname << "'\n";
      return 0;
    }

    if (! _output.new_stem(fname)) {
      cerr << "Cannot initialise output stream '" << fname << "'\n";
      return 0;
    }

    if (verbose) {
      cerr << "output to '" << fname << "'\n";
    }
  }

  if (cl.option_present('U')) {
    IWString fname = cl.string_value('U');
    if (_failures.would_overwrite_input_files(cl, fname)) {
      cerr << "Cannot overwrite input file(s) -U '" << fname << "'\n";
      return 0;
    }

    if (! _failures.new_stem(fname)) {
      cerr << "Cannot initialise failure stream '" << fname << "'\n";
      return 0;
    }

    if (verbose) {
      cerr << "failed molecules to '" << fname << "'\n";
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

class CoreReplacement {
  private:
    int _verbose;

    int _molecules_read;

    int _reduce_to_largest_fragment;

    int _remove_chirality;

    // There are two ways in which a molecule can succeed.
    //  1. If we find a conformer where all atoms are separated by _distance_threshold
    //  2. The conformation that has the longest, closest separation.
    float _distance_threshold;
    int _longest_separation;

    // With what resolution do we can the angle(s);
    float _angle_delta;

    extending_resizable_array<int> _rotations_performed;

    // The first two matched atoms of the first query to match is scanned.
    resizable_array_p<Substructure_Query> _query;

    // What to do if the none of the queries match.
    int _ignore_molecules_not_matching_queries;
    int _write_molecules_not_matching_queries;
    int _molecules_not_matching_queries;

    int _report_failed_scans = 0;

    int _successful_calculation = 0;
    int _initial_conformation_ok = 0;

    // There are a great many acceptance criteria.

    FileType _input_type;

    Chemical_Standardisation _chemical_standardisation;


  // Private functions.
    int GetMatchedAtoms(Molecule& m, atom_number_t& a1, atom_number_t& a2);
    int Process(Molecule& m, atom_number_t a1, atom_number_t a2);
    int ShowFailedDistances(const Molecule& m, std::ostream& output) const;

  public:
    CoreReplacement();

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);

    // Convert `m` to a form that passes distance constraints.
    // If it fails, `m` is left with an undetermined geometry.
    int Process(Molecule& m);

    int Report(std::ostream& output) const;
};

CoreReplacement::CoreReplacement() {
  _verbose = 0;
  _molecules_read = 0;
  _successful_calculation = 0;
  _report_failed_scans = 0;
  _initial_conformation_ok = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _angle_delta = iwmisc::Deg2Rad(30.0f);
  _distance_threshold = 0.0;
  _ignore_molecules_not_matching_queries = 0;
  _write_molecules_not_matching_queries = 0;
  _molecules_not_matching_queries = 0;
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

  if (! cl.option_present('T')) {
    cerr << "Must specify threshold via the -T option\n";
    return 0;
  }
  if (cl.option_present('T')) {
    if (! cl.value('T', _distance_threshold) || _distance_threshold <= 0.0f) {
      cerr << "Invalid distance threshold (-T)\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will fail bump checks at separation " << _distance_threshold << '\n';
    }
  }

  if (cl.option_present('r')) {
    if (! cl.value('r', _angle_delta) || _angle_delta < 1.0f) {
      cerr << "Invalid angle (-r)\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Rotation by " << _angle_delta << " degrees\n";
    }
    _angle_delta = iwmisc::Deg2Rad(_angle_delta);
  }

  if (cl.option_present('s')) {
    const_IWSubstring s;
    for (int i = 0; cl.value('s', s, i); ++i) {
      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
      if (! q->create_from_smarts(s)) {
        cerr << "Invalid smarts '" << s << "'\n";
        return 0;
      }
      _query << q.release();
    }
  }

  if (cl.option_present('q')) {
    if (! process_queries(cl, _query, _verbose, 'q')) {
      cerr << "Cannot process command line queries (-q)\n";
      return 0;
    }
  }

  if (_query.empty()) {
    cerr << "No queries\n";
    return 0;
  }

  for (Substructure_Query* q : _query) {
    q->set_find_unique_embeddings_only(1);
    q->set_perceive_symmetry_equivalent_matches(0);
  }

  if (cl.option_present('b')) {
    _report_failed_scans = 1;
    if (_verbose) {
      cerr << "Will report distances in failed scans\n";
    }
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

  if (m.natoms() == 0) {
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

  if (_molecules_not_matching_queries == _molecules_read) {
    output << "All molecules matched a query\n";
  } else {
    output << _molecules_not_matching_queries << " did not match a query\n";
  }

  output << _initial_conformation_ok << " initial conformations passed constraints\n";
  output << _successful_calculation << " successful calculations " << iwmisc::Fraction<float>(_successful_calculation, _molecules_read) << '\n';

  return 1;
}

int
CoreReplacement::GetMatchedAtoms(Molecule& m,
                                 atom_number_t& a1,
                                 atom_number_t& a2) {
  Molecule_to_Match target(&m);
  Substructure_Results sresults;
  for (Substructure_Query* q : _query) {
    if (! q->substructure_search(target, sresults)) {
      continue;
    }
    a1 = sresults.embedding(0)->item(0);
    a2 = sresults.embedding(0)->item(1);
    return 1;
  }

  return 0;
}

int
CoreReplacement::Process(Molecule& m) {
  ++_molecules_read;

  atom_number_t a1, a2;
  if (! GetMatchedAtoms(m, a1, a2)) {
    ++_molecules_not_matching_queries;
    return 0;
  }

  return Process(m, a1, a2);
}

int
UpdateClosestApproach(const Molecule& m,
                      const int* classification,
                      float& global_max) {
  const int matoms = m.natoms();
  float local_min = std::numeric_limits<float>::max();
  for (int i = 0; i < matoms; ++i) {
    for (int j = i + 1; j < matoms; ++j) {
      if (classification[i] == classification[j]) {
        continue;
      }
      float d = m.distance_between_atoms(i, j);
      if (d < local_min) {
        local_min = d;
      }
    }
  }

  if (local_min <= global_max) {
    return 0;
  }

  global_max = local_min;
  return 1;
}

int
CoreReplacement::Process(Molecule& m,
                         atom_number_t a1,
                         atom_number_t a2) {
  assert(_distance_threshold > 0.0f);

  const int matoms = m.natoms();
  std::unique_ptr<int[]>sides(new_int(matoms));
  if (! m.identify_side_of_bond(sides.get(), a1, 1, a2)) {
    cerr << "CoreReplacement::Process:cannot identify atoms on either side of bond " << m.name() << '\n';
    m.ring_membership();
    cerr << " atom " << a1 << ' ' << m.smarts_equivalent_for_atom(a1) << '\n';
    cerr << " atom " << a2 << ' ' << m.smarts_equivalent_for_atom(a2) << '\n';
    return 0;
  }

#ifdef CHECK_COLLISIONS_HERE
  int nshort = 0;
  for (int i = 0; i < matoms; ++i) {
    for (int j = i + 1; j < matoms; ++j) {
      if (sides[i] == sides[j]) {
        continue;
      }
      const float d = m.distance_between_atoms(i, j);
      if (d < _distance_threshold) {
        cerr << "   short " << j << ' ' << d << '\n';
        ++nshort;
      }
    }
  }
#endif

  // If already OK,
  if (m.any_bump_check(sides.get(), _distance_threshold) == 0) {
    ++_successful_calculation;
    ++_initial_conformation_ok;
    return 1;
  }

  Space_Vector<double> axis{m.x(a1) - m.x(a2), m.y(a1) - m.y(a2), m.z(a1) - m.z(a2)};
  axis.normalise();

  // We want to record the conformer with the longest, shortest distance.
  float longest_distance = 0.0f;
  // The conformer that has that shortest distance.
  std::unique_ptr<float[]> best_conf;

  const int nconf = static_cast<int>(2.0 * M_PI / _angle_delta);
  for (int i = 1; i < nconf; ++i) {
    m.rotate_atoms(axis, static_cast<double>(i * _angle_delta), sides.get());
    // If we are looking for the longest separation, check each one.
    if (_longest_separation) {
      if (UpdateClosestApproach(m, sides.get(), longest_distance)) {
        best_conf = m.GetCoords();
      }
    } else {
      // Any conformation that works is accepted.
      if (m.any_bump_check(sides.get(), _distance_threshold) == 0) {
        ++_successful_calculation;
        return 1;
      }
    }
  }

  // if we were looking for anything that passed the bump check, and we get to 
  // here, we have failed.
  if (! _longest_separation) {
    if (_report_failed_scans) {
      ShowFailedDistances(m, cerr);
    }
    return 0;
  }

  // Set the atoms to the best found.
  m.SetXyz(best_conf.get());
  ++_successful_calculation;
  return 1;
}


int
CoreReplacement::ShowFailedDistances(const Molecule& m, std::ostream& output) const {
  cerr << "Failed calculation " << m.name() << " last\n";
  const int matoms = m.natoms();

  for (int i = 0; i < matoms; ++i) {
    for (int j = i + 1; j < matoms; ++j) {
      if (m.are_bonded(i, j)) {
        continue;
      }
      float d = m.distance_between_atoms(i, j);
      if (d >= _distance_threshold) {
        continue;
      }
      output << " btw " << i << " and " << j << " dist " << d << '\n';
    }
  }

  return 1;
}

int
ReplaceCore(CoreReplacement& core_replacement,
            Options& options,
            Molecule& m) {
  std::unique_ptr<float[]> initial_coords;

  if (options.FailuresActive()) {
    initial_coords = m.GetCoords();
  }

  if (core_replacement.Process(m)) {
    return options.WriteSuccess(m);
  }

  if (options.FailuresActive()) {
    m.SetXyz(initial_coords.get());
    return options.WriteFailed(m);
  }

  return 1;
}

int
ReplaceCore(CoreReplacement& core_replacement,
            Options& options,
            data_source_and_type<Molecule>& input) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! core_replacement.Preprocess(*m)) {
      return 0;
    }

    if (! ReplaceCore(core_replacement, options, *m)) {
      return 0;
    }
  }

  return 1;
}

int
ReplaceCore(CoreReplacement& core_replacement,
            Options& options,
            const char * fname) {
  FileType input_type = options.input_type();
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return ReplaceCore(core_replacement, options, input);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:i:g:lcq:s:S:U:r:T:fz:b");
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

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  for (const char* fname : cl) {
    if (! ReplaceCore(core_replacement, options, fname)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    core_replacement.Report(cerr);
  }

  return 0;
}

}  // namespace conformational_scan

int
main(int argc, char** argv) {
  int rc = conformational_scan::Main(argc, argv);

  return rc;
}
