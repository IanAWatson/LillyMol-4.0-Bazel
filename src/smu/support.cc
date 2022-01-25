#include "smu/support.h"

namespace smu {

using std::cerr;

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

}  // namespace smu
