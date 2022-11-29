// Align a set of 3D molecules based on a substructure query.
// One matched atom is placed at the origin, one along the X axis
// and another as close to the Y axis as possible.

#include <functional>
#include <iostream>
#include <memory>
#include <optional>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

namespace align_molecules {

using std::cerr;

void
Usage(int rc) {

  cerr << "Aligns molecules based on query matches\n";
  cerr << " -q <query>        specify query\n";
  cerr << " -s <smarts>       query as smarts\n";
  cerr << " -O <atoms>        the matched atom(s) that define the origin\n";
  cerr << " -X <atoms>        the matched atom(s) that define the X axis\n";
  cerr << " -Y <atoms>        the matched atom(s) that define the Y axis\n";
  cerr << " -C <angle>        rather than rotate an atom to the Y axis, rotate the closest atom to <angle> degrees from Y\n";
  cerr << " -z i              ignore molecules not matching the query\n";
  cerr << " -z w              write molecules not matching the query\n";
  cerr << " -T <x,y,z>        apply a final translation to molecules just before writing\n";
  cerr << " -c                remove chirality\n";
  cerr << " -l                strip to largest fragment\n";
  cerr << " -v                verbose output\n";

  ::exit(rc);
}

class AlignByMatchedAtoms {
  private:
    int _verbose;

    int _molecules_read;

    int _reduce_to_largest_fragment;

    int _remove_chirality;

    resizable_array_p<Substructure_Query> _query;

    int _ignore_molecules_not_matching_queries;
    int _write_molecules_not_matching_queries;

    int _molecules_not_matching_queries;

    // Rather than a single matched atom specifying the location
    // we allow multiple matched atom numbers. Not sure if that is
    // a good idea or not.
    resizable_array<int> _origin;
    resizable_array<int> _xaxis;
    resizable_array<int> _yaxis;

    // Rather than a matched atom aligned to the Y axis, place the atom with the highest X coordinate
    // at an angle relative to Y
    std::optional<float> _yangle;

    // After alignment has been done, we can apply a translation before writing.
    Space_Vector<float> _final_translation;
    // To enable a quick check on whether or not _final_translation is set.
    int _final_translation_active;

    FileType _input_type;

    Chemical_Standardisation _chemical_standardisation;

  // Private functions.

    int RunQueries(Molecule& m, Substructure_Results& sresults);
    int Process(Molecule& m, const Substructure_Results& sresults, Molecule_Output_Object& output);
    int HandleMoleculesNotMatchingQueries(Molecule& m, Molecule_Output_Object& output);
    int RotateClosestToYaxis(Molecule& m, const Set_of_Atoms& embedding,
                        Molecule_Output_Object& output);
    int Write(Molecule& m, Molecule_Output_Object& output);

  public:
    AlignByMatchedAtoms();

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);

    int Process(Molecule& m, Molecule_Output_Object& output);

    int Report(std::ostream& output) const;

    FileType input_type() const {
      return _input_type;
    }
};

AlignByMatchedAtoms::AlignByMatchedAtoms() {
  _verbose = 0;
  _molecules_read = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _ignore_molecules_not_matching_queries = 0;
  _write_molecules_not_matching_queries = 0;
  _molecules_not_matching_queries = 0;
  _final_translation_active = 0;
  _input_type = FILE_TYPE_INVALID;
}

// Fetch the integer values associated with the `flag` option
// and place them in `destination`.
// Note no error or duplicate checking.
// Returns true even if no values are specified - on purpose.
int
GetMatchedAtoms(Command_Line& cl,
                char flag,
                resizable_array<int>& destination) {
  int value;
  for (int i = 0; cl.value(flag, value, i); ++i) {
    destination << value;
  }

  return 1;
}

int
AlignByMatchedAtoms::Initialise(Command_Line& cl) {
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

  if (! cl.option_present('O')) {
    cerr << "Must specify the origin matched atom number(s) via the -O option\n";
    return 0;
  }

  if (! GetMatchedAtoms(cl, 'O', _origin)) {
    cerr << "Cannot get origin matched atoms (-O)\n";
    return 0;
  }

  if (! GetMatchedAtoms(cl, 'X', _xaxis)) {
    cerr << "Cannot get x axis matched atoms (-X)\n";
    return 0;
  }

  if (! GetMatchedAtoms(cl, 'Y', _yaxis)) {
    cerr << "Cannot get Y axis matched atoms (-Y)\n";
    return 0;
  }

  if (cl.option_present('T')) {
    const_IWSubstring t = cl.option_value('T');
    if (! _final_translation.read(t, ',')) {
      cerr << "Invalid final translation specification -T '" << t << "'\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will translate before writing by " << _final_translation << '\n';
    }
    _final_translation_active = 1;
  }

  if (cl.option_present('C')) {
    float angle;
    if (! cl.value('C', angle)) {
      cerr << "Invalid closest atom Y axis angle (-C)\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will place the highest X atom at " << angle << " degrees from the Y axis\n";
    }
    _yangle = angle * DEG2RAD;
  }

  if (1 == cl.number_elements() && 0 == strcmp("-", cl[0])) { // reading a pipe, assume sdf
    _input_type = FILE_TYPE_SDF;
  } else if (!all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot discern all file types, use the -i option\n";
    return 0;
  } else if (!process_input_type(cl, _input_type)) {
    return 0;
  }

  return 1;
}

int
AlignByMatchedAtoms::Preprocess(Molecule& m) {
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
AlignByMatchedAtoms::Write(Molecule& m,
                           Molecule_Output_Object& output) {
  if (! _final_translation_active) {
    return output.write(m);
  }

  m.translate_atoms(_final_translation);

  return output.write(m);
}

Space_Vector<float>
Center(const Molecule& m,
       const Set_of_Atoms& embedding,
       const resizable_array<int>& matched_atoms) {
  Space_Vector<float> result;

  for (int i : matched_atoms) {
    atom_number_t j = embedding[i];
    result += m.atom(j);
  }

  result /= static_cast<float>(matched_atoms.number_elements());

  return result;
}

int
AlignByMatchedAtoms::Process(Molecule& m,
                         const Substructure_Results& sresults,
                         Molecule_Output_Object& output) {
  Space_Vector<float> origin = Center(m, *sresults.embedding(0), _origin);
  m.translate_atoms(-origin);

  if (_xaxis.empty()) {
    return Write(m, output);
  }

  Space_Vector<float> xatoms = Center(m, *sresults.embedding(0), _xaxis);
  xatoms.normalise();

  static const Space_Vector<float> xaxis(1.0f, 0.0f, 0.0f);
  static const Space_Vector<float> yaxis(0.0f, 1.0f, 0.0f);
  static const Space_Vector<float> zaxis(0.0f, 0.0f, 1.0f);

  angle_t theta = xaxis.angle_between_unit_vectors(xatoms);

  // The axis about which we will rotate.
  xatoms.cross_product(xaxis);
  xatoms.normalise();

  m.translate_atoms(- xatoms);
  m.rotate_atoms(xatoms, theta);
  m.translate_atoms(xatoms);

  if (_yaxis.empty() && ! _yangle) {
    return Write(m, output);
  }

  if (_yangle) {
    return RotateClosestToYaxis(m, *sresults.embedding(0), output);
  }

  Space_Vector<float> yatoms = Center(m, *sresults.embedding(0), _yaxis);
  // cerr << "Before X rotation Y at " << yatoms << '\n';

  // We need to set the x coordinate to zero, so make a temp.
  // Not really necessary, but clearer/safer...
  Space_Vector<float> tmp(yatoms);
  tmp.set_x(0.0f);
  tmp.normalise();
  theta = yaxis.angle_between_unit_vectors(tmp);

  // If close enough, just write.
  if (abs(yatoms.z()) < 1.0e-03) {
    return Write(m, output);
  }

  if (yatoms.z() < 0.0f) {
    m.rotate_atoms(xaxis, theta);
  } else {
    m.rotate_atoms(xaxis, -theta);
  }

  Space_Vector<float> after_rotation = Center(m, *sresults.embedding(0), _yaxis);
  // cerr << "After  X rotation Y at " << after_rotation << '\n';
  if (abs(after_rotation.z()) < 1.0e-03) {
    return Write(m, output);
  }

  cerr << "Rotation to Y seems to have failed " << after_rotation << ' ' << m.name() << '\n';

  return Write(m, output);
}

void
Scatter(const resizable_array<int>& from,
        const Set_of_Atoms& embedding,
        int* to,
        int value) {
  for (int x : from) {
    to[embedding[x]] = value;
  }
}

//#define DEBUG_ROTATE_CLOSEST_TO_Y_AXIS
// by definition, the highest X value must come from one of the atoms that did NOT
// define the origin or the X axis.
// Note that Hydrogen atoms are not considered.
int
AlignByMatchedAtoms::RotateClosestToYaxis(Molecule& m, const Set_of_Atoms& embedding,
                        Molecule_Output_Object& output) {

  const int matoms = m.natoms();
  // Keep track of which atoms are to be considered. Note that in the common case
  // of both origin and x axis each being defined by a single matched atom, this
  // is quite inefficient.
  std::unique_ptr<int[]> include_atom(new_int(matoms, 1));
  Scatter(_origin, embedding, include_atom.get(), 0);
  Scatter(_xaxis, embedding, include_atom.get(), 0);

  float highest_x_coord = -std::numeric_limits<int>::max();
  atom_number_t highest_x_atom = INVALID_ATOM_NUMBER;
  for (int i = 0; i < matoms; ++i) {
    if (! include_atom[i]) {
      continue;
    }
    if (m.atomic_number(i) == 1) {
      continue;
    }
    const float x = m.x(i);
    if (x < highest_x_coord) {
      continue;
    }
    highest_x_coord = x;
    highest_x_atom = i;
  }

#ifdef DEBUG_ROTATE_CLOSEST_TO_Y_AXIS
  cerr << "Highest X coord atom " << highest_x_atom << " " << m.smarts_equivalent_for_atom(highest_x_atom) << " at " << highest_x_coord << " " << m.name() << '\n';
#endif

  // The case where all atoms were involved with defining the origin and x axis.
  if (highest_x_atom == INVALID_ATOM_NUMBER) {
    return output.write(m);
  }

  static const Space_Vector<float> xaxis(1.0f, 0.0f, 0.0f);
  static const Space_Vector<float> yaxis(0.0f, 1.0f, 0.0f);

  static const Space_Vector<float> target(0.0f, cos(*_yangle), sin(*_yangle));

  Space_Vector<float> closest(m.x(highest_x_atom), m.y(highest_x_atom), m.z(highest_x_atom));

  Space_Vector<float>tmp(closest);
  tmp.set_x(0.0f);
  tmp.normalise();
  float current_angle = target.angle_between_unit_vectors(tmp);
  float delta = current_angle - *_yangle;
  tmp.cross_product(target);

  if (tmp.x() < 0.0f) {
    m.rotate_atoms(xaxis, -delta);
  } else {
    m.rotate_atoms(xaxis, delta);
  }

  Space_Vector<float> after_rotation(m.x(highest_x_atom), m.y(highest_x_atom), m.z(highest_x_atom));
  after_rotation.set_x(0.0f);
  after_rotation.normalise();
  const float new_angle = after_rotation.angle_between_unit_vectors(target);
#ifdef DEBUG_ROTATE_CLOSEST_TO_Y_AXIS
  cerr << "require " << *_yangle << " currently " << current_angle << '\n';
#endif
  if (abs(new_angle - *_yangle) < 0.0001) {
#ifdef DEBUG_ROTATE_CLOSEST_TO_Y_AXIS
    static int ngood = 0;
    ++ngood;
    cerr << "good " << ngood << '\n';
#endif
    return Write(m, output);
  }

  static int nbad = 0;
  ++nbad;
  cerr << "RotateClosestToYaxis:rotation did not work, expected " << *_yangle << " got " << new_angle << " tmp.x " << tmp.x() << '\n';
  cerr << "Initial angle " << current_angle << " delta " << delta << " bad " << nbad << '\n';
  cerr << "Final coords " << m.get_coords(highest_x_atom) << '\n';

  return Write(m, output);
}

int
AlignByMatchedAtoms::RunQueries(Molecule& m, Substructure_Results& sresults) {
  Molecule_to_Match target(&m);
  for (Substructure_Query* q : _query) {
    if (q->substructure_search(target, sresults)) {
      return 1;
    }
  }

  return 0;
}

int
AlignByMatchedAtoms::Process(Molecule& m,
                         Molecule_Output_Object& output) {
  ++_molecules_read;

  Substructure_Results sresults;
  if (! RunQueries(m, sresults)) {
    return HandleMoleculesNotMatchingQueries(m, output);
  }

  return Process(m, sresults, output);
}

int
AlignByMatchedAtoms::Report(std::ostream& output) const {
  output << "AlignByMatchedAtoms:read " << _molecules_read << " molecules\n";
  if (_molecules_read == 0) {
    return 1;
  }
  if (_molecules_not_matching_queries == 0) {
    cerr << "All molecules matched a query\n";
  } else {
    output << _molecules_not_matching_queries << " molecules matched no queries\n";
  }

  return 1;
}
    
int
AlignByMatchedAtoms::HandleMoleculesNotMatchingQueries(Molecule& m,
                Molecule_Output_Object& output) {
  ++_molecules_not_matching_queries;
  if (_ignore_molecules_not_matching_queries) {
    return 1;
  }

  if (_write_molecules_not_matching_queries) {
    return output.write(m);
  }

  cerr << m.name() << " no matches to any of " << _query.size() << " substructure queries\n";
  return 0;
}

int
AlignMolecules(AlignByMatchedAtoms& core_replacement,
            Molecule& m,
            Molecule_Output_Object& output) {
  return core_replacement.Process(m, output);
}

int
AlignMolecules(AlignByMatchedAtoms& core_replacement,
            data_source_and_type<Molecule>& input,
            Molecule_Output_Object& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! core_replacement.Preprocess(*m)) {
      return 0;
    }

    if (! AlignMolecules(core_replacement, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
AlignMolecules(AlignByMatchedAtoms& core_replacement,
            const char * fname,
            Molecule_Output_Object& output) {
  FileType input_type = core_replacement.input_type();
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return AlignMolecules(core_replacement, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:i:g:lcq:s:O:X:Y:o:S:z:T:C:");
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

  AlignByMatchedAtoms core_replacement;
  if (! core_replacement.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  if (! cl.option_present('S')) {
    cerr << "Must specify output file via the -S option\n";
    Usage(1);
  }

  Molecule_Output_Object output;
  if (cl.option_present('o')) {
    if (! output.determine_output_types(cl, 'o')) {
      cerr << "Cannot determine output types (-o)\n";
      return 1;
    }
  } else {
    output.add_output_type(FILE_TYPE_SDF);
  }

  if (cl.option_present('S')) {
    IWString fname = cl.string_value('S');
    if (output.would_overwrite_input_files(cl, fname)) {
      cerr << "Cannot overwrite input(s) '" << fname << "'\n";
      return 1;
    }

    if (! output.new_stem(fname)) {
      cerr << "Cannot set up output stream '" << fname << "'\n";
      return 1;
    }

    if (verbose) {
      cerr << "Output written to '" << fname << "'\n";
    }
  }

  for (const char* fname : cl) {
    if (! AlignMolecules(core_replacement, fname, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    core_replacement.Report(cerr);
  }

  return 0;
}

}  // namespace align_molecules

int
main(int argc, char** argv) {
  int rc = align_molecules::Main(argc, argv);

  return rc;
}
