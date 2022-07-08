// Demo of app to demonstrate etracting information from a LillyMol molecule

#include <iostream>

#include "Molecule_Tools/coords_and_connectivity.h"

namespace cached_molecule {

using std::cerr;

int
CachedMolecule::Build(const std::string& smiles) {
  if (! _mol.build_from_smiles(smiles)) {
    // some kind of python appropriate error message
    cerr << "CachedMolecule::Build:invalid smiles '" << smiles << "'\n";
    return false;
  }

  return _mol.natoms();
}

std::optional<int>
CachedMolecule::ParseSmiles(const std::string& smiles) {
  if (! _mol.build_from_smiles(smiles)) {
    // some kind of python appropriate error message
    cerr << "CachedMolecule::Build:invalid smiles '" << smiles << "'\n";
    return std::nullopt;
  }

  return _mol.natoms();
}

int
CachedMolecule::GetCoords(std::vector<float>& coords) const {
  const int matoms = _mol.natoms();  // Number of atoms in _mol.

  coords.resize(matoms * 3);  // Each atom has 3 coordinates.

  if (matoms == 0) {
    return 0;
  }

  // Loop over all the atoms in _mol.
  for (int i = 0; i < matoms; ++i) {
    const Atom* a = _mol.atomi(i);
    coords[3 * i] = a->x();
    coords[3 * i + 1] = a->y();
    coords[3 * i + 2] = a->z();
  }

  return matoms;
}

// Convert `btype` to an integer. Single_bond->1...
int
BondTypeToInt(bond_type_t btype) {
  switch(btype) {
    case SINGLE_BOND:
      return 1;
    case DOUBLE_BOND:
      return 2;
    case TRIPLE_BOND:
      return 3;
    default:
      cerr << "BondTypeToInt:what is " << btype << '\n';
      return 0;
  }
}

int
CachedMolecule::ConnectionMatrix(std::vector<int>& connections) const {
  const int matoms = _mol.natoms();

  connections.resize(matoms * matoms, 0);

  if (matoms == 0) {
    return 0;
  }

  // Loop over all the bonds (edges) in the molecule graph.
  for (const Bond* b : _mol.bond_list()) {
    // Get the two atom numbers that define the bond.
    const int a1 = b->a1();
    const int a2 = b->a2();
    const int btype = BondTypeToInt(b->btype());
    connections[a1 * matoms + a2] = btype;
    connections[a2 * matoms + a1] = btype;
  }

  return matoms;
}

}  // namespace cached_molecule
