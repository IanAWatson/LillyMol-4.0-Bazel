#ifndef MOLECULE_TOOLS_ATOM_SEPARATION_H
#define MOLECULE_TOOLS_ATOM_SEPARATION_H

#include <cstdint>

#include "leveldb/db.h"

#include "Foundational/iwstring/iwstring.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/molecule.h"

namespace separated_atoms {

// We need to hash atom type - distance combinations. This
// class generates hash functions for those.
class Hasher {
  private:
    int _atomic_number_to_index[HIGHEST_ATOMIC_NUMBER + 1];

  public:
    Hasher();

    uint32_t Value(int z1, int z2, int distance) const;
    uint32_t Value(int z1, int z2, int z3, int d12, int d13, int d23) const;

    int OkAtomicNumber(int z) const {
      return _atomic_number_to_index[z] > 0;
    }

    // When writing files, it is convenient to have element
    // names, rather than numbers.
    // Given a 2 bond hash, convert to <element1><distance><element2>
    IWString TwoHashToString(uint32_t) const;
};

// Return a vector for each bond in `m` indicating whether or not
// that bond is unique - as derived from symmetry. For symmetric bonds
// the first one will be marked as 1, subsequent instances as 0
std::unique_ptr<int[]> BondSymmetry(Molecule& m);

int TurnOffAmideBonds(Molecule& m, int * breakable);

// Given a hash value, return a key for a leveldb database.
leveldb::Slice Key(const uint32_t& hash);

}  // namespace separated_atoms

#endif // MOLECULE_TOOLS_ATOM_SEPARATION_H
