#ifndef SMU_SUPPORT_H
#define SMU_SUPPORT_H

#include <optional>

#include "smu/dataset.pb.h"
#include "Molecule_Lib/molecule.h"

namespace smu {

std::optional<Molecule> MoleculeFromBondTopology(const GoogleSmu::BondTopology& bond_topology);

int AddGeometry(const GoogleSmu::Geometry& geometry, Molecule& m);

std::optional<int> AtomTypeToAtomicNumber(GoogleSmu::BondTopology::AtomType atype);

double DistanceBetweenAtoms(const GoogleSmu::Geometry& geom, int a1, int a2);

// Return the index of the BondTopology in `btids` that is the
// starting topology.
int IndexOfStartingBTid(const google::protobuf::RepeatedPtrField<GoogleSmu::BondTopology>& btids);

}  // namespace smu
#endif // SMU_SUPPORT_H
