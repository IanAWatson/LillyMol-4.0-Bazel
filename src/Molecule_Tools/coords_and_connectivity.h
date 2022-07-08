#ifndef MOLECULE_TOOLS_COORDS_AND_CONNECTIVITY_H
#define MOLECULE_TOOLS_COORDS_AND_CONNECTIVITY_H

#include <optional>
#include <string>
#include <vector>

#include "Molecule_Lib/molecule.h"

namespace cached_molecule {

// Class to do one-time ingestion of a molecule, which is then available for
// querying.
class CachedMolecule {
  private:
    // The molecules that is the source of data. Ingested via Build().
    Molecule _mol;

  public:
    // Build _mol from `smiles`. Returns the number of atoms in `smiles`.
    // Returns -1 if `smiles` is invalid.
    int Build(const std::string& smiles);
    // Try to build _mol from `smiles`. If successful, returns the number
    // of atoms.
    // Needs a different name than Build...
    std::optional<int> ParseSmiles(const std::string& smiles);

    // The number of atoms in _mol.
    int natoms() const {
      return _mol.natoms();
    }
    // Any Molecule member function could be implemented.

    // Fill `coords` with 3*natoms coordinates from _mol.
    // The vector is ordered with x,y,z of the first atom, then
    // xyz of the second atom...
    // Returns the number of atoms. coords.size() == 3 * natoms;
    int GetCoords(std::vector<float>& coords) const;

    // Fill `connections` with natoms*natoms data about connectivity.
    // If there is a single bond between two atoms,the entry on
    // `connections` will be 1, double bond-> 2, triple bond->3
    // No provision for aromatic - although that would be easy.
    int ConnectionMatrix(std::vector<int>& connections) const;
};

}  // namespace cached_molecule
#endif  // MOLECULE_TOOLS_COORDS_AND_CONNECTIVITY_H
