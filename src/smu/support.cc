#include <cmath>

#include "smu/support.h"

namespace smu {

constexpr float kBohrToAngstrom = 0.529177f;

using std::cerr;

using GoogleSmu::BondTopology;
using GoogleSmu::Geometry;

std::optional<Molecule>
MoleculeFromBondTopology(const BondTopology& bond_topology) {
  Molecule result;

  for (int i = 0; i < bond_topology.atoms().size(); ++i) {
    const BondTopology::AtomType atype = bond_topology.atoms(i);
    switch (atype) {
      case BondTopology::ATOM_H: {
          const Element * e = get_element_from_atomic_number(1);
          result.add(e);
        }
        break;
      case BondTopology::ATOM_C: {
        const Element * e = get_element_from_atomic_number(6);
        result.add(e);
      }
        break;
      case BondTopology::ATOM_N: {
        const Element * e = get_element_from_atomic_number(7);
        result.add(e);
      }
        break;
      case BondTopology::ATOM_NPOS: {
        const Element * e = get_element_from_atomic_number(7);
        result.add(e);
        result.set_formal_charge(result.natoms() - 1, 1);
      }
        break;
      case BondTopology::ATOM_O: {
        const Element * e = get_element_from_atomic_number(8);
        result.add(e);
      }
        break;
      case BondTopology::ATOM_ONEG: {
        const Element * e = get_element_from_atomic_number(8);
        result.add(e);
        result.set_formal_charge(result.natoms() - 1, -1);
      }
        break;
      case BondTopology::ATOM_F: {
        const Element * e = get_element_from_atomic_number(9);
        result.add(e);
      }
        break;
      case BondTopology::ATOM_UNDEFINED: {
        cerr << "Undefined atom\n";
        return std::nullopt;
      }
        break;
      default:
        cerr << "Unrecognised atom type " << atype << '\n';
        return std::nullopt;
    }
  }

  for (int i = 0; i < bond_topology.bonds().size(); ++i) {
    const BondTopology::Bond& bond = bond_topology.bonds(i);
    const atom_number_t a1 = bond.atom_a();
    const atom_number_t a2 = bond.atom_b();
    switch (bond.bond_type()) {
      case BondTopology::BOND_SINGLE:
        result.add_bond(a1, a2, SINGLE_BOND, 1);
        break;
      case BondTopology::BOND_DOUBLE:
        result.add_bond(a1, a2, DOUBLE_BOND, 1);
        break;
      case BondTopology::BOND_TRIPLE:
        result.add_bond(a1, a2, TRIPLE_BOND, 1);
        break;
      default:
        cerr << "Unrecognised bond type " << bond.bond_type() << '\n';
        return  std::nullopt;
    }
  }

  return result;
}

int
AddGeometry(const Geometry& geometry,
            Molecule& m) {
  const int matoms = m.natoms();
  if (matoms != geometry.atom_positions().size()) {
    return 0;
  }

  for (int i = 0; i < matoms; ++i) {
    const Geometry::AtomPos& apos = geometry.atom_positions(i);
    coord_t x = apos.x() * kBohrToAngstrom;
    coord_t y = apos.y() * kBohrToAngstrom;
    coord_t z = apos.z() * kBohrToAngstrom;
    m.setxyz(i, x, y, z);
  }

  return 1;
}

std::optional<int>
AtomTypeToAtomicNumber(BondTopology::AtomType atype) {
  switch (atype) {
    case BondTopology::ATOM_H:
      return 1;
    case BondTopology::ATOM_C:
      return 6;
    case BondTopology::ATOM_N:
      return 7;
    case BondTopology::ATOM_NPOS:
      return 7;
    case BondTopology::ATOM_O:
      return 8;
    case BondTopology::ATOM_ONEG:
      return 8;
    case BondTopology::ATOM_F:
      return 9;
    case BondTopology::ATOM_UNDEFINED:
      cerr << "AtomTypeToAtomicNumber:unrecognised " << atype << '\n';
      return std::nullopt;
    default:
        cerr << "AtomTypeToAtomicNumber::invalid atom type " << atype << '\n';
        return std::nullopt;
  }
}

double
DistanceBetweenAtoms(const Geometry& geometry, int a1, int a2)
{
  const Geometry::AtomPos pos1 = geometry.atom_positions(a1);
  const Geometry::AtomPos pos2 = geometry.atom_positions(a2);
  return kBohrToAngstrom * sqrt(
        (pos1.x() - pos2.x()) * (pos1.x() - pos2.x()) +
        (pos1.y() - pos2.y()) * (pos1.y() - pos2.y()) +
        (pos1.z() - pos2.z()) * (pos1.z() - pos2.z())
  );
}

int
IndexOfStartingBTid(const google::protobuf::RepeatedPtrField<BondTopology>& btids) {
  for (int i = 0; i < btids.size(); ++i) {
    if (btids[i].is_starting_topology()) {
      return i;
    }
  }
  //cerr << "Warning starting bond topology not found\n";
  return -1;
}

}  // namespace smu
