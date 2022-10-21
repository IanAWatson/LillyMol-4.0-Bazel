// Scan a series of molecules likely from dicer, and generate
// an index based on inter-attachment point distances.
// If dicer is run with the `-I env` option, the attachment points
// are isotopically labelled.

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwbits/fixed_bit_vector.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

//#include "Molecule_Tools/atom_separations.h"

namespace frag_to_map {

using std::cerr;
using fixed_bit_vector::FixedBitVector;

void
Usage(int rc) {

  cerr << " -c                remove chirality\n";
  cerr << " -l                strip to largest fragment\n";
  cerr << " -v                verbose output\n";

  ::exit(rc);
}

// To avoid passing around a lot of arguments, gather info about
// the molecule in one place.
class MolData {
  private:
    int _matoms;

    // The atom types assigned to the molecule.
    std::unique_ptr<uint32_t> _atom_type;
  
    // An array matoms*matoms set if the bond between those two atoms
    // is breakable.
    std::unique_ptr<int[]> _breakable_bond;

    // Set if the i'th atom has at least one breakble bond.
    std::unique_ptr<int[]> _has_breakable_bond;

    // a stack of atoms that during recursive function calls, holds
    // the id's of the atoms already included.
    // Note that this object does not retain the stack pointer.
    std::unique_ptr<int[]> _atom;

    // One of the operations we perform is to identify the inner atoms
    // defined by a set of atoms. that info will be stored in this array.
    std::unique_ptr<int[]> _subset;
    
    // When creating subsets, we need a cross reference.
    std::unique_ptr<int[]> _xref;

  // private functions
    // Allocate the arrays for a molecule with `matoms` atoms.
    int IdentifyAtomsBetweenStackMembers2(Molecule& m);
    int IdentifyAtomsBetweenStackMembers3(Molecule& m);

    atom_number_t AtomHeadingTowardsTarget(Molecule& m,
                        atom_number_t origin,
                        atom_number_t destination) const;
    atom_number_t AtomHeadingTowardsTarget(Molecule& m,
                        atom_number_t origin,
                        atom_number_t destination1, atom_number_t destination2) const;

    int IdentifyAtomsBetween(Molecule& m, atom_number_t zatom, atom_number_t destination);
    int IdentifyAtomsBetween(Molecule& m, atom_number_t zatom,
                              atom_number_t destination1, atom_number_t destination2);

  public:
    MolData(int max_connections, int matoms);

    int IdentifyBreakableBonds(Molecule& m);

    // Given `stack_size` atoms in the stack, populate
    // the _subset array with the atoms that are between
    // the stack atoms.
    int IdentifyAtomsBetweenStackMembers(Molecule& m, int stack_size);

    uint32_t* atom_type() {
      return _atom_type.get();
    }
    const uint32_t* atom_type() const {
      return _atom_type.get();
    }
    
    int* atom_stack() {
      return _atom.get();
    }
    const int* atom_stack() const {
      return _atom.get();
    }

    int HasBreakableBond(int a) {
      return _has_breakable_bond[a];
    }

    int IsBreakable(int a1, int a2) const {
      return _breakable_bond[a1 * _matoms + a2];
    }

    int* xref() {
      return _xref.get();
    }
    int* subset() {
      return _subset.get();
    }
};

MolData::MolData(int max_connections, int matoms) {
  _matoms = matoms;

  _atom.reset(new_int(max_connections + 1));

  _has_breakable_bond.reset(new_int(matoms));
  _breakable_bond.reset(new_int(matoms * matoms));
  _atom_type.reset(new_unsigned_int(matoms));
  _subset.reset(new_int(matoms));
  _xref.reset(new_int(matoms));
}

int
IsAmideLike(Molecule& m,
            const atom_number_t nitrogen,
            const atom_number_t carbon) {
  const Atom * nitrogen_atom = m.atomi(nitrogen);
  if (nitrogen_atom->atomic_number() != 7) {
    return 0;
  }

  const Atom * carbon_atom = m.atomi(carbon);
  if (carbon_atom->atomic_number() == 6) {
  } else if (carbon_atom->atomic_number() == 16) {
  } else {
    return 0;
  }

  if (carbon_atom->ncon() == carbon_atom->nbonds()) {
    return 0;
  }

  for (const Bond* b : *carbon_atom) {
    if (! b->is_double_bond()) {
      continue;
    }
    atom_number_t o = b->other(carbon);
    if (m.atomic_number(o) == 8) {
      return 1;
    }
  }

  return 0;
}

// Returns the number of atoms with breakable bonds.
int
MolData::IdentifyBreakableBonds(Molecule& m) {
  const int matoms = m.natoms();

  m.compute_aromaticity_if_needed();

  int rc = 0;
  for (const Bond* b : m.bond_list()) {
    if (b->nrings()) {
      continue;
    }

    if (! b->is_single_bond()) {
      continue;
    }

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();
    if (IsAmideLike(m, a1, a2) ||
        IsAmideLike(m, a2, a1)) {
      continue;
    }
    _breakable_bond[a1 * matoms + a2] = 1;
    _breakable_bond[a2 * matoms + a1] = 1;
    ++_has_breakable_bond[a1];
    ++_has_breakable_bond[a2];
    ++rc;
  }

  return rc;
}

int
MolData::IdentifyAtomsBetweenStackMembers(Molecule& m, int stack_size) {
  if (stack_size == 2) {
    return IdentifyAtomsBetweenStackMembers2(m);
  } else if (stack_size == 3) {
    return IdentifyAtomsBetweenStackMembers3(m);
  } else {
    cerr << "MolData::IdentifyAtomsBetweenStackMembers:how to process " << stack_size << '\n';
    return 0;
  }
}

// Return the one atom that is bonded to `origin` and
// is closer to `destination`. If there are multiple
// such atoms, we return INVALID_ATOM_NUMBER.
atom_number_t
MolData::AtomHeadingTowardsTarget(Molecule& m,
                        atom_number_t origin,
                        atom_number_t destination) const {
  const Atom* a = m.atomi(origin);
  int dmax = m.bonds_between(origin, destination);
  atom_number_t rc = INVALID_ATOM_NUMBER;
  for (const Bond* b : *a) {
    if (b->nrings()) {
      continue;
    }

    atom_number_t j = b->other(origin);
    if (m.bonds_between(j, destination) < dmax) {
      if (rc == INVALID_ATOM_NUMBER) {
        rc = j;
      } else {
        return INVALID_ATOM_NUMBER;
      }
    }
  }

  return rc;
}

// Return the atom number of the single atom bonded to `origin` that is
// along the way to both destinations.
atom_number_t
MolData::AtomHeadingTowardsTarget(Molecule& m,
                        atom_number_t origin,
                        atom_number_t destination1,
                        atom_number_t destination2) const {
  const Atom* a = m.atomi(origin);
  int dmax1 = m.bonds_between(origin, destination1);
  int dmax2 = m.bonds_between(origin, destination2);
  atom_number_t rc = INVALID_ATOM_NUMBER;
  for (const Bond* b : *a) {
    if (b->nrings()) {
      continue;
    }

    atom_number_t j = b->other(origin);
    if (_subset[j]) {
      continue;
    }
    if (m.bonds_between(j, destination1) >= dmax1 ||
        m.bonds_between(j, destination2) >= dmax2) {
      continue;
    }

    if (rc == INVALID_ATOM_NUMBER) {
      rc = j;
    } else {
      return INVALID_ATOM_NUMBER;
    }
  }

  return rc;
}

// First task is to determine whether or not there is a unique
// path from both atoms to the other. This avoids the case
// of a ring starting at either end.
int
MolData::IdentifyAtomsBetweenStackMembers2(Molecule& m) {
  atom_number_t a0 = _atom[0];
  atom_number_t a1 = _atom[1];
  cerr << "Atoms are " << a0 << " and " << a1 << '\n';
  if (m.bonds_between(a0, a1) < 3) {
    return 0;
  }

  atom_number_t towards_a1 = AtomHeadingTowardsTarget(m, a0, a1);
  if (towards_a1 == INVALID_ATOM_NUMBER) {
    return 0;
  }

  atom_number_t towards_a0 = AtomHeadingTowardsTarget(m, a1, a0);
  if (towards_a0 == INVALID_ATOM_NUMBER) {
    return 0;
  }
  cerr << " towards " << a0 << " " << towards_a1 << ' ' << a1 << '\n';
  cerr << " towards " << a1 << " " << towards_a0 << ' ' << a0 << '\n';
  if (towards_a0 == towards_a1) {
    return 0;
  }

  std::fill_n(_subset.get(), m.natoms(), 0);

  _subset[a0] = 2;
  _subset[a1] = 2;

  int rc = IdentifyAtomsBetween(m, towards_a1, a1);
  cerr << "From IdentifyAtomsBetween " << rc << '\n';

  _subset[a0] = 0;
  _subset[a1] = 0;

  _subset[towards_a0] = 2;
  _subset[towards_a1] = 2;

  return rc;
}

int
MolData::IdentifyAtomsBetween(Molecule& m,
                              atom_number_t zatom,
                              atom_number_t destination) {
  const int d = m.bonds_between(zatom, destination);
  _subset[zatom] = 1;

  const Atom* a = m.atomi(zatom);
  for (const Bond* b : *a) {
    atom_number_t j = b->other(zatom);
    if (_subset[j]) {
      continue;
    }

    if (j == destination) {
      return 1;
    }

    if (m.bonds_between(j, destination) >= d) {
      continue;
    }

    IdentifyAtomsBetween(m, j, destination);
  }

  return 1;
}

int
MolData::IdentifyAtomsBetween(Molecule& m,
                              atom_number_t zatom,
                              atom_number_t destination1,
                              atom_number_t destination2) {
  const int d1 = m.bonds_between(zatom, destination1);
  const int d2 = m.bonds_between(zatom, destination2);

  _subset[zatom] = 1;

  const Atom* a = m.atomi(zatom);
  for (const Bond* b : *a) {
    atom_number_t j = b->other(zatom);
    if (_subset[j]) {
      continue;
    }

    if (j == destination1 || j == destination2) {
      return 1;
    }

    if (m.bonds_between(j, destination1) >= d1 ||
        m.bonds_between(j, destination2) >= d2) {
      continue;
    }

    IdentifyAtomsBetween(m, j, destination1, destination2);
  }

  return 1;
}


int
MolData::IdentifyAtomsBetweenStackMembers3(Molecule& m) {
  atom_number_t a1 = _atom[0];
  atom_number_t a2 = _atom[1];
  atom_number_t a3 = _atom[2];

  // For each of a1, a2 and a3, identify an attached atom that heads
  // off towards the other two atoms.
  atom_number_t from_a1 = AtomHeadingTowardsTarget(m, a1, a2, a3);
  if (from_a1 == INVALID_ATOM_NUMBER) {
    return 0;
  }
  atom_number_t from_a2 = AtomHeadingTowardsTarget(m, a2, a1, a3);
  if (from_a2 == INVALID_ATOM_NUMBER) {
    return 0;
  }
  atom_number_t from_a3 = AtomHeadingTowardsTarget(m, a3, a2, a3);
  if (from_a3 == INVALID_ATOM_NUMBER) {
    return 0;
  }

  if (! IdentifyAtomsBetween(m, a1, a2, a3)) {
    return 0;
  }

  return 1;
}


class GenerateFragments {
  private:
    int _isotope;

    int _min_fragment_size;
    int _max_fragment_size;

    int _fragments_generated;

    // We keep track of which subsets of atoms have been processed before.
    resizable_array<FixedBitVector> _seen;

  // private functions
  int Process2(Molecule& m, MolData& mol_data, int natoms);
  int Process3(Molecule& m, MolData& mol_data, int natoms);

  public:
    GenerateFragments();

    int Process(Molecule& m, MolData& mol_data, int natoms);

    int Report(std::ostream& output) const;
};

GenerateFragments::GenerateFragments() {
  _isotope = 1;
  
  _min_fragment_size = 0;
  _max_fragment_size = std::numeric_limits<int>::max();

  _fragments_generated = 0;
}

// The `natoms` atoms in the atom_stack of `mol_data` form
// a set of OkDistance atoms.
int
GenerateFragments::Process(Molecule& m,
                           MolData& mol_data,
                           int natoms) {
  if (natoms == 2) {
    return Process2(m, mol_data, natoms);
  } else if (natoms == 3) {
    return Process3(m, mol_data, natoms);
  }

  return 1;
}

int
GenerateFragments::Process2(Molecule& m,
                            MolData& mol_data,
                            int natoms) {
  cerr << "Process2 callin IdentifyAtomsBetweenStackMembers\n";
  if (! mol_data.IdentifyAtomsBetweenStackMembers(m, natoms)) {
    return 0;
  }

  const int* atoms = mol_data.atom_stack();
  for (int i = 0; i < m.natoms(); ++i) {
    cerr << " atom " << i << " subset " << mol_data.subset()[i] << '\n';
  }

  const int d = m.bonds_between(atoms[0], atoms[1]);
  uint32_t at1 = mol_data.atom_type()[atoms[0]];
  uint32_t at2 = mol_data.atom_type()[atoms[1]];
  if (at1 > at2) {
    std::swap(at1, at2);
  }

  int* subset = mol_data.subset();
  int* xref = mol_data.subset();

  const int matoms = m.natoms();

  int atoms_in_subset = 0;

  // identify the join points, their value in _subset is 2
  atom_number_t attach1 = INVALID_ATOM_NUMBER;
  atom_number_t attach2 = INVALID_ATOM_NUMBER;
  for (int i = 0; i < matoms; ++i) {
    if (subset[i] == 0) {
      continue;
    }
    if (subset[i] == 1) {
      ++atoms_in_subset;
      continue;
    }
    if (subset[i] != 2) {
      continue;
    }
    ++atoms_in_subset;
    subset[i] = 1;
    if (attach1 == INVALID_ATOM_NUMBER) {
      attach1 = i;
    } else if (attach2 == INVALID_ATOM_NUMBER) {
      attach2 = i;
      break;
    }
  }

  if (atoms_in_subset < _min_fragment_size ||
      atoms_in_subset > _max_fragment_size) {
    return 0;
  }

  if (attach2 == INVALID_ATOM_NUMBER) {
    cerr << "GenerateFragments:Process2:no attachments\n";
    return 0;
  }

  Molecule frag;
  m.create_subset(frag, subset, 1, xref);

  frag.set_isotope(xref[attach1], _isotope);
  frag.set_isotope(xref[attach2], _isotope);

  ++_fragments_generated;

  cerr << frag.smiles() << '\n';

  return 1;
}

int
GenerateFragments::Process3(Molecule& m,
                            MolData& mol_data,
                            int natoms) {
  const int* atoms = mol_data.atom_stack();

  if (! mol_data.IdentifyAtomsBetweenStackMembers(m, 3)) {
    return 0;
  }
  atom_number_t a1 = mol_data.atom_stack()[0];
  atom_number_t a2 = mol_data.atom_stack()[1];
  atom_number_t a3 = mol_data.atom_stack()[2];

  int d12 = m.bonds_between(a1, a2);
  int d13 = m.bonds_between(a1, a3);
  int d23 = m.bonds_between(a2, a3);
  uint32_t at1 = mol_data.atom_type()[atoms[0]];
  uint32_t at2 = mol_data.atom_type()[atoms[1]];
  uint32_t at3 = mol_data.atom_type()[atoms[2]];
  // Manual sort.
  if (at1 > at2) {
    std::swap(at1, at2);
    std::swap(d13, d23);
  }
  if (at2 > at3) {
    std::swap(at2, at3);
    std::swap(d12, d13);
  }
  if (at1 > at2) {
    std::swap(at1, at2);
    std::swap(d13, d23);
  }

  std::fill_n(mol_data.subset(), m.natoms(), 0);

  return 1;
}

int
GenerateFragments::Report(std::ostream& output) const {
  output << "Generated " << _fragments_generated << " fragments\n";
  return 1;
}

class FragToMap {
  private:
    int _verbose;

    int _molecules_processed;

    int _reduce_to_largest_fragment;

    int _remove_chirality;

    FileType _input_type;

    Chemical_Standardisation _chemical_standardisation;

    Atom_Typing_Specification _atom_type;

    int _min_distance;
    int _max_distance;

    int _max_connections;

    // Some current molecule specific arrays. No thread safety here.

    std::unique_ptr<uint32_t> _atype;

    int _molecules_with_no_breakable_bonds;

  // private functions

    int OkDistance(const int d) const;
    template <typename P> int GenerateAtomPairs(Molecule& m, P& proc);
    template <typename P> int Recurse(Molecule& m,
                   MolData& mol_data, int nprev, P& proc);

  public:
    FragToMap();

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);

    int Report(std::ostream& output) const;

    FileType input_type() const {
      return _input_type;
    }

    // Generate atom pairs, triples, etc and pass each one
    // to `func`. During profiling a collection, func will be
    // something that accumulates data. During a lookup it will
    // be something that looks up the new data in a database.
    template <typename P> int Process(Molecule& m, P& proc);
};

FragToMap::FragToMap() {
  _verbose = 0;
  _molecules_processed = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _input_type = FILE_TYPE_INVALID;

  _min_distance = 0;
  _max_distance = std::numeric_limits<int>::max();

  _max_connections = 3;

  _molecules_with_no_breakable_bonds = 0;
}

int
FragToMap::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('c')) {
    _remove_chirality = 1;
    if (_verbose) {
      cerr << "Will remove chirality from input molecules\n";
    }
  }

  if (cl.option_present('g')) {
    if (! _chemical_standardisation.construct_from_command_line(cl, _verbose, 'g')) {
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

  if (cl.option_present('P')) {
    IWString p = cl.string_value('P');
    if (! _atom_type.build(p)) {
      cerr << "Cannot initialise atom typing (-P) '" << p << "'\n";
      return 0;
    }
  } else {
    IWString p = "UST:ACY";
    _atom_type.build(p);
  }

  if (cl.option_present('m')) {
    if (! cl.value('m', _min_distance) || _min_distance < 1) {
      cerr << "The min distance option (-m) must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will only process atoms at least " << _min_distance << " bonds apart\n";
    }
  }

  if (cl.option_present('M')) {
    if (! cl.value('m', _max_distance) || _max_distance < _min_distance) {
      cerr << "The max distance option (-M) must be a whole +ve number greater than " << _min_distance << '\n';
      return 0;
    }

    if (_verbose) {
      cerr << "Will only process atoms at most " << _max_distance << " bonds apart\n";
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
FragToMap::Preprocess(Molecule& m) {
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
FragToMap::OkDistance(const int d) const {
  if (d > _max_distance) {
    return 0;
  }

  if (d < _min_distance) {
    return 0;
  }

  return 1;
}

int
FragToMap::Report(std::ostream& output) const {
  cerr << "FragToMap:processed " << _molecules_processed << " molecules\n";
  if (_molecules_processed == 0) {
    return 1;
  }

  return 1;
}

template <typename P>
int
FragToMap::Process(Molecule& m,
                   P& proc) {
  ++_molecules_processed;

  return GenerateAtomPairs(m, proc);
}

template <typename P>
int
FragToMap::GenerateAtomPairs(Molecule& m,
                             P& proc) {
  const int matoms = m.natoms();

  MolData mol_data(_max_connections, matoms);

  if (! _atom_type.assign_atom_types(m, mol_data.atom_type())) {
    return 0;
  }

  if (! mol_data.IdentifyBreakableBonds(m)) {
    ++_molecules_with_no_breakable_bonds;
    return 1;
  }

  for (int i = 0; i < matoms; ++i) {
    if (! mol_data.HasBreakableBond(i)) {
      continue;
    }
    mol_data.atom_stack()[0] = i;
    for (int j = i + 1; j < matoms; ++j) {
      if (! mol_data.HasBreakableBond(j)) {
        continue;
      }
      const int d = m.bonds_between(i, j);
      if (d == 1) {
        continue;
      }

      if (! OkDistance(d)) {
        continue;
      }
      mol_data.atom_stack()[1] = j;
      cerr << " atoms are " << i << " and " << j << " dist " << d << '\n';
      proc.Process(m, mol_data, 2);

      if (_max_connections > 2) {
        Recurse(m, mol_data, 2, proc);
      }
    }
  }

  return 1;
}

int
FragToMap::BreakBonds(Molecule& m, int depth) {
  const int matoms = m.natoms();

  MolData mol_data(_max_connections, matoms);

  if (! _atom_type.assign_atom_types(m, mol_data.atom_type())) {
    return 0;
  }

  if (! mol_data.IdentifyBreakableBonds(m)) {
    ++_molecules_with_no_breakable_bonds;
    return 1;
  }

  for (const Bond* b : m.bond_list()) {
    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();
    if (! mol_data.IsBreakable(a1, a2)) {
      continue;
    }
    Molecule mcopy(m);
    mcopy.remove_bond_between_atoms(a1, a2);
    mcopy.set_isotope(a1, depth);
    mcopy.set_isotope(a2, depth);
  }

  return 1;
}

template <typename P>
int
FragToMap::Recurse(Molecule& m,
                   MolData& mol_data,
                   int nprev,
                   P& proc) {
  const int matoms = m.natoms();

  int* prev = mol_data.atom_stack();

  for (int i = prev[nprev - 1] + 1; i < matoms; ++i) {
    for (int j = 0; j < nprev; ++j) {
      if (! OkDistance(m.bonds_between(prev[j], i))) {
        continue;;
      }
      prev[nprev] = i;
    }
    proc.Process(m, mol_data, nprev + 1);
    if (nprev < _max_connections) {
      Recurse(m, mol_data, nprev + 1, proc);
    }
  }

  return 1;
}



// Load `breakable_bond` with 1 for each atom pair that bounds a
// brekable bond. For every atom that has at least one end of a 
// breakable bond, set `has_breakable_bond`.
// We assume the arrays are zero'd on entry.
int
ReplaceCore(FragToMap& frag_to_map,
            Molecule& m,
            GenerateFragments& proc) {
  return frag_to_map.Process(m, proc);
}

int
ReplaceCore(FragToMap& frag_to_map,
            data_source_and_type<Molecule>& input,
            GenerateFragments& proc) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! frag_to_map.Preprocess(*m)) {
      continue;
    }

    if (! ReplaceCore(frag_to_map, *m, proc)) {
      return 0;
    }
  }

  return 1;
}

int
ReplaceCore(FragToMap& frag_to_map,
            const char * fname,
            GenerateFragments& proc) {
  FileType input_type = frag_to_map.input_type();
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return ReplaceCore(frag_to_map, input, proc);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:i:g:lcP:");
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

  FragToMap frag_to_map;
  if (! frag_to_map.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  GenerateFragments proc;
  for (const char* fname : cl) {
    if (! ReplaceCore(frag_to_map, fname, proc)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    frag_to_map.Report(cerr);
  }

  return 0;
}

}  // namespace frag_to_map
int
main(int argc, char** argv) {
  int rc = frag_to_map::Main(argc, argv);

  return rc;
}
