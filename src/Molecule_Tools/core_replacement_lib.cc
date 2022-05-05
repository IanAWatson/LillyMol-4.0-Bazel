#include <iostream>

#include "Molecule_Tools/core_replacement_lib.h"

namespace separated_atoms {

// Returns only true or false.
int
IdentifyInterior(const Molecule& m, const atom_number_t zatom,
                 int * interior) {
#ifdef DEBUG_IDENTIFY_INTERIOR
  cerr << "IdentifyInterior continues with atom " << zatom << " d from start " << distance_from_start << '\n';
#endif

  const Atom& a = m.atom(zatom);
  for (const Bond* b : a) {
    const atom_number_t o = b->other(zatom);
    // If we have already identified this atom, continue.
    if (interior[o] == 2) {
      continue;
    }
    // We have found an outer atom. If adjacent to where we start, that
    // is fine. Otherwise we may have found a ring, and must fail.
    if (interior[o] == 1) {
      continue;
    }
    interior[o] = 2;
    if (! IdentifyInterior(m, o, interior)) {
      return 0;
    }
  }

#ifdef DEBUG_IDENTIFY_INTERIOR
  cerr << "From atom " << zatom << " returning\n";
#endif
  return 1;
}

// Identify the single atom, attached to `zatom` that is 'inside' the
// region defined by `matched_atoms`. `zatom` will be a member of `matched_atoms`.
// Works by looking for atoms that going from zatom gets closer to all the
// other members of matched_atoms.
atom_number_t
IdentifyFrontierAtom(Molecule& m,
              atom_number_t zatom,
              const Set_of_Atoms& matched_atoms,
              const int * to_be_removed) {
  const Atom& atom = m.atom(zatom);
#ifdef DEBUG_IDENTIFY_FRONTIER_ATOM
  cerr << "Looking for interior frontier from " << zatom << '\n';
#endif

  // The atom(s) that seem to be heading towards the other matched atoms from zatom.
  // There must be just one such atom.
  atom_number_t frontier = INVALID_ATOM_NUMBER;  // To be returned.
  for (const Bond* bond : atom) {
    atom_number_t f = bond->other(zatom);
#ifdef DEBUG_IDENTIFY_FRONTIER_ATOM
    cerr << "Check atom " << f << " tbr " << to_be_removed[f] << '\n';
#endif
    if (to_be_removed[f]) {  // Presumably another member of matched_atoms
      continue;
    }

    bool looks_interior = true;
    for (atom_number_t ma2: matched_atoms) {
      if (ma2 == zatom) {
        continue;
      }
#ifdef DEBUG_IDENTIFY_FRONTIER_ATOM
      cerr << "From " << zatom << " to " << ma2 << " dist " << m.bonds_between(f, ma2) << " and " << m.bonds_between(zatom, ma2) << '\n';
#endif
      if (m.bonds_between(f, ma2) > m.bonds_between(zatom, ma2)) {
        looks_interior = false;
        break;
      }
    }
    if (! looks_interior) {
      continue;
    }
    if (frontier == INVALID_ATOM_NUMBER) {
      frontier = f;
    } else {
      return INVALID_ATOM_NUMBER;
    }
  }

  if (frontier == INVALID_ATOM_NUMBER) {
    return INVALID_ATOM_NUMBER;
  }

  return frontier;
}

#define DEBUG_IDENTIFY_ATOMS_TO_REMOVE

// Identify the atoms that are 'inside' the region defined by
// the atoms in `matched_atoms`.
// The first task is, for each atom in `matched_atoms`, identify
// the first atom in the inner part of the molecule.
// There must be just one such atom for each item in `matched_atoms`.
int
IdentifyAtomsToRemove(Molecule& m,
                      const Set_of_Atoms& matched_atoms,
                      int * to_be_removed) {
  matched_atoms.set_vector(to_be_removed, 1);   // will later be reset.
  Set_of_Atoms frontier;
  for (atom_number_t ma1 : matched_atoms) {
    cerr << "Find frontier from atom " << ma1 << '\n';
    atom_number_t f = IdentifyFrontierAtom(m, ma1, matched_atoms, to_be_removed);
    if (f == INVALID_ATOM_NUMBER) {
      cerr << "No frontier atom from " << ma1 << '\n';
      return 0;
    }
    frontier.add_if_not_already_present(f);
    cerr << "Atom " << f << " is frontier from " << ma1 << '\n';
  }

#ifdef DEBUG_IDENTIFY_ATOMS_TO_REMOVE
  cerr << "After identifying frontier\n";
  for (int i = 0; i < m.natoms(); ++i) {
    cerr << "atom " << i << " value " << to_be_removed[i] << " type " << m.smarts_equivalent_for_atom(i) << '\n';
  }
#endif

  if (frontier.empty()) {  // Should not happen.
    return 0;
  }

  frontier.set_vector(to_be_removed, 2);

  // Expand all frontier atoms into the center.
  for (atom_number_t f : frontier) {
    if (IdentifyInterior(m, f, to_be_removed) < 0) {
      return 0;
    }
  }

#ifdef DEBUG_IDENTIFY_ATOMS_TO_REMOVE
  cerr << "After finding interior\n";
  for (int i = 0; i < m.natoms(); ++i) {
    cerr << "atom " << i << " value " << to_be_removed[i] << " type " << m.smarts_equivalent_for_atom(i) << '\n';
  }
#endif

  return 1;
}

}  // namespace separated_atoms
