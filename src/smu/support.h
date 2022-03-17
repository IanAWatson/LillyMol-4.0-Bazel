#ifndef SMU_SUPPORT_H
#define SMU_SUPPORT_H

#include <optional>

#include "smu/dataset.pb.h"
#include "Molecule_Lib/molecule.h"

namespace smu {

std::optional<Molecule> MoleculeFromBondTopology(const BondTopology& bond_topology);

int AddGeometry(const Geometry& geometry, Molecule& m);

std::optional<int> AtomTypeToAtomicNumber(BondTopology::AtomType atype);

double DistanceBetweenAtoms(const Geometry& geom, int a1, int a2);

}  // namespace smu
#endif // SMU_SUPPORT_H
