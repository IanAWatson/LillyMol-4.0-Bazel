#include <algorithm>
#include <iostream>
#include <unordered_set>

#include "Foundational/iwmisc/misc.h"

#include "atom_separations.h"

namespace separated_atoms {

using std::cerr;

Hasher::Hasher() {
  std::fill_n(_atomic_number_to_index, HIGHEST_ATOMIC_NUMBER + 1, -1);
  _atomic_number_to_index[6] = 1;
  _atomic_number_to_index[7] = 2;
  _atomic_number_to_index[8] = 3;
  _atomic_number_to_index[15] = 4;
  _atomic_number_to_index[16] = 5;
}

// Return a hash function for atomic numbers `z1` and `z2` separated
// by `distance` bonds.
// In the case of a pair of atomic numbers, we form a 4 digit number
// where the first digit is the lowest atomic number, the middle 2
// are the distance, and the last one is the highest atomic number
uint32_t
Hasher::Value(int z1, int z2, int distance) const {
  if (_atomic_number_to_index[z1] < 0 ||
      _atomic_number_to_index[z2] < 0) {
    cerr << "Hasher::Value:invalid atomic number " << z1 << " " << z2 << '\n';
    return 0;
  }

  if (z1 > z2) {
    std::swap(z1, z2);
  }

  return _atomic_number_to_index[z1] * 1000 + distance * 10 + _atomic_number_to_index[z2];
}

// For a triple, we first sort the arguments to a canonical order.
uint32_t
Hasher::Value(int z1, int z2, int z3, int d12, int d13, int d23) const {
  z1 = _atomic_number_to_index[z1];
  z2 = _atomic_number_to_index[z2];
  z3 = _atomic_number_to_index[z3];
  if (z1 > z2) {
    std::swap(z1, z2);
    std::swap(d13, d23);
  }

  if (z2 > z3) {
    std::swap(z2, z3);
    std::swap(d12, d13);
  }

  if (z1 > z2) {
    std::swap(z1, z2);
    std::swap(d13, d23);
  }

  if (z1 == z2 && z2 == z3) {
    if (d12 > d13) {
      std::swap(d12, d13);
    }
    if (d13 > d23) {
      std::swap(d23, d12);
    }
    if (d12 > d23) {
      std::swap(d12, d23);
    }
  } else if (z1 == z2) {
    if (d13 > d23) {
      std::swap(d13, d23);
    }
  } else if (z2 == z3) {
    if (d12 > d13) {
      std::swap(d12, d13);
    }
  }

  // Use base 40 to encode the distances, 40^3 = 64000, so atoms start with 100k
  uint32_t rc = 100000 * 100 * z1 +
                100000 * 10 * z2 +
                100000 * z3 +
                40 * 40 * d12 +
                40 * d13 +
                d23;
  return rc;
}

std::unique_ptr<int[]>
BondSymmetry(Molecule& m) {
  const int nedges = m.nedges();
  std::unique_ptr<int[]> result(new_int(nedges, 1));

  const int * symmetry_class = m.symmetry_classes();
#ifdef DEBUG_INITIALISE_BONDS
  for (int i = 0; i < m.natoms(); ++i) {
    cerr << " atom " << i << " sym " << symmetry_class[i] << '\n';
  }
#endif

  // The number of times each symmetry class is found.
  const int matoms = m.natoms();
  std::unique_ptr<int[]> degeneracy(new_int(matoms + 1));
  int highest_value = 0;
  for (int i = 0; i < matoms; ++i) {
    ++degeneracy[symmetry_class[i]];
    if (symmetry_class[i] > highest_value) {
      highest_value = symmetry_class[i];
    }
  }

  // If no symmetry we are done.
  if (highest_value == matoms) {
    return result;
  }

  std::unordered_set<int> seen;
  for (int i = 0; i < nedges; ++i) {
    const Bond * b = m.bondi(i);
    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();
    if (degeneracy[symmetry_class[a1]] == 1 && degeneracy[symmetry_class[a2]] == 1) {
      continue;
    }
    int minsym, maxsym;
    if (symmetry_class[a1] <= symmetry_class[a2]) {
      minsym = symmetry_class[a1];
      maxsym = symmetry_class[a2];
    } else {
      minsym = symmetry_class[a2];
      maxsym = symmetry_class[a1];
    }
    int canonical = maxsym * matoms + minsym;
    auto [iter, inserted] = seen.emplace(canonical);
#ifdef DEBUG_INITIALISE_BONDS
    cerr << b->a1() << ',' << b->a2() << " minsym " << minsym << " maxsym " << maxsym << " canonical " << canonical << " inserted " << inserted << '\n';
#endif

    if (! inserted) {  // Seen before, is a symmetry duplicate.
      result[i] = 0;
    }
  }

  return result;
}

#ifdef THIS_VERSION_SLIGHTLY_SLOWER
int
DoublyBondedToOxygen(Molecule& m,
                     atom_number_t zatom) {
  const Atom * a = m.atomi(zatom);
  for (const Bond * b : *a) {
    if (b->is_single_bond()) {
      continue;
    }
    if (m.atomic_number(b->other(zatom)) == 8) {
      return 1;
    }
  }

  return 0;
}

int
TurnOffAmideBondsII(Molecule& m,
                  int * breakable) {
  const int nedges = m.nedges();
  int rc = 0;
  for (int i = 0; i < nedges; ++i) {
    const Bond * b = m.bondi(i);
    if (b->nrings()) {
      continue;
    }
    if (! b->is_single_bond()) {
      continue;
    }
    atom_number_t carbon = INVALID_ATOM_NUMBER;  // or sulphur
    if (m.atomic_number(b->a1()) == 7 && m.ncon(b->a1())) {
      carbon = b->a2();
    } else if (m.atomic_number(b->a2()) == 7 && m.ncon(b->a2())) {
      carbon = b->a1();
    } else {
      continue;
    }
    
    const Atom* c = m.atomi(carbon);
    if (c->fully_saturated()) {
      continue;
    }

    if (c->atomic_number() == 6 && c->ncon() == 3) {
    } else if (c->atomic_number() == 16 && c->ncon() == 4) {
    } else {
      continue;
    }

    if (DoublyBondedToOxygen(m, carbon)) {
      breakable[i] = 0;
      ++rc;
    }
  }

  return rc;
}
#endif

// Return true if atom `zatom` in `m` is singly bonded to
// a Nitrogen atom.
int
SingleBondToNitrogen(const Molecule& m,
                     atom_number_t zatom) {
  const Atom * atom = m.atomi(zatom);
  for (const Bond* b : *atom) {
    if (! b->is_single_bond()) {
      continue;
    }
    if (m.atomic_number(b->other(zatom)) == 7) {
      return 1;
    }
  }

  return 0;
}

// For each bond in `m` turn of the corresponding entry in 
// `breakable` if that bond is the C-N bond of an amide (or sulfonamide).
int
TurnOffAmideBonds(Molecule& m,
                  int * breakable) {
  const int matoms = m.natoms();
  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    const Atom * o = m.atomi(i);
    if (o->atomic_number() != 8) {
      continue;
    }
    if (o->ncon() != 1) {
      continue;
    }
    const Bond * b = o->item(0);
    if (! b->is_double_bond()) {
      continue;
    }
    const Atom* carbon = m.atomi(b->other(i));
    if (carbon->fully_saturated()) {
      continue;
    }
    if (carbon->atomic_number() == 6 && carbon->ncon() == 3) {
    } else if (carbon->ncon() == 4 && carbon->atomic_number() == 16) {
    } else {
      continue;
    }

    if (SingleBondToNitrogen(m, b->other(i))) {
      breakable[i] = 0;
      ++rc;
    }
  }

  return rc;
}

leveldb::Slice 
Key(const uint32_t& hash) {
  leveldb::Slice result(reinterpret_cast<const char*>(&hash), sizeof(hash));
  return result;
}

}  // namespace separated_atoms
