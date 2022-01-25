#ifndef SMU_SUPPORT_H
#define SMU_SUPPORT_H

#include <optional>

#include "smu/dataset.pb.h"
#include "Molecule_Lib/molecule.h"

namespace smu {

std::optional<Molecule> MoleculeFromBondTopology(const BondTopology& bond_topology);

}  // namespace smu
#endif // SMU_SUPPORT_H
