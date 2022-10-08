// An introduction to LillyMol

#include <iostream>

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

namespace lillymol_introcution {

using std::cerr;

// An empty molecule has no atoms or bonds, but can have a name.
void
DemoEmptyMolecule() {
  Molecule mol;
  mol.set_name("foo");
  // Should have zero atoms and molecular weight.
  cerr << "Empty molecule '" << mol.name() << "' has " << mol.natoms()
       << " atoms, amw " << mol.molecular_weight() << '\n';
}

// When presented with an invalid smiles, build_from_smiles will fail.
void
DemoCannotParseBadSmiles() {
  Molecule mol;
  if (mol.build_from_smiles("invalid")) {
    cerr << "Building from bad smiles succeeded, this should not happen.\n";
  } else {
    cerr << "Bad smiles cannot be parsed, good outcome.\n";
  }
}

void
DemoAtomicNumbers() {
  Molecule mol;
  if (!mol.build_from_smiles("C")) {
    cerr << "Unable to build methane!\n";
    return;
  }
  mol.set_name("methane");

  cerr << mol.name() << " has " << mol.natoms() << " atoms\n";

  // There are three ways to find the atomic number of atom atom.
  // Ask the molecule for the atomic number of a given atom.
  cerr << "Molecule's atomic number 0 " << mol.atomic_number(0) << '\n';

  // We can fetch an Atom and ask the atom for its atomic number.
  // We only have one atom in the molecule, so it is atom number 0.
  const Atom& a = mol.atom(0);
  cerr << "Atom's atomic number " << a.atomic_number() << '\n';

  // We can ask the atom for its Element and query the atomic number
  const Element* e = a.element();
  cerr << "Element's atomic number " << e->atomic_number() << '\n';
}

void
DemoAtomIteration() {
  Molecule mol;
  if (! mol.build_from_smiles("CCC")) {
    cerr << "Cannot parse propane\n";
    return;
  }

  // Loop through the atoms in `mol`, writing the atomic number.
  // Note that the iterator returns a pointer to an Atom.
  for (const Atom* atom : mol) {
    // We can ask an atom for its atomic number of atomic symbol.
    cerr << " atomic_number " << atom->atomic_number() << '\n';
  }

  // Or we can iterate based on atom numbers. This time we
  // can get a reference to an atom.
  const int matoms = mol.natoms();
  for (int i = 0; i < matoms; ++i) {
    const Atom& atom = mol.atom(i);
    cerr << " atomic number " << atom.atomic_number() << '\n';
  }
}

void
DemoNcon() {
  Molecule mol;
  mol.build_from_smiles("CC(C)(C)C");

  // The molecule is
  //              [C:2]
  //                |
  //                |
  //     [C:0] -- [C:1] -- [C:4]
  //                |
  //                |
  //              [C:3]
  // Atoms have different connectivity.

  // We can get the degree of any atom by asking either the
  // molecule, or the atom.
  // In addition get a smarts equivalent for that atom - this
  // is mostly useful for debugging.

  const int matoms = mol.natoms();
  for (int i = 0; i < matoms; ++i) {
    const Atom& atom = mol.atom(i);
    cerr << " atom " << i << " has " << mol.ncon(i) << " connections, atom "
         << atom.ncon() << " type " << mol.smarts_equivalent_for_atom(i) << '\n';
  }
}

void
DemoBondIteration() {
  Molecule mol;
  mol.build_from_smiles("CC(C)(C)C");

  // For each atom, enumerate the connected atoms.
  // Iteration is over Bond objects associated with each atom.
  // The Bond contains two atom numbers, we use the bond->other
  // method to determine what the other atom is.

  const int matoms = mol.natoms();
  for (int i = 0; i < matoms; ++i) {
    const Atom& atom = mol.atom(i);
    for (const Bond* bond : atom) {
      cerr << "atom " << i << " bonded to " << bond->other(i) <<
              " type " << bond->btype() << '\n';
    }
  }
}

// Start with ethane, verify one fragment, then remove the bond
// to generate two methanes.
void
DemoFragments() {
  Molecule mol;
  mol.build_from_smiles("CC");

  cerr << "Methane has " << mol.number_fragments() << " fragments\n";
  cerr << "Ethane contains only one edge " << mol.nedges() << '\n';

  // there are only two atoms, remove the bond between them.
  mol.remove_bond_between_atoms(0, 1);

  // We should have two fragments now.
  cerr << "After bond removal " << mol.number_fragments() << " fragments\n";

  // Fragment membership is only known by the Molecule.
  // Fragments are assigned sequential numbers, starting with zero.
  // We can ask the molecule, in which fragment is the i'th atom?
  // Atom's do not know anything about fragment membership.

  const int matoms = mol.natoms();
  for (int i = 0; i < matoms; ++i) {
    cerr << " atom " << i << " in fragment " << mol.fragment_membership(i) << '\n';

    // const Atom& atom = mol.atom(i);
    // atom.fragment_membership() CANNOT WORK
  }

  // we can ask the molecule how many atoms are in each fragment.
  const int nf = mol.number_fragments();
  for (int i = 0; i < nf; ++i) {
    cerr << " Fragment " << i << " contains " << mol.atoms_in_fragment(i) << " atoms\n";
  }

  // To determine if two atoms are in the same fragment, compare 
  // mol.fragment_membership(i) == mol.fragment_membership(j);

  // Put the molecule back by adding a bond

  mol.add_bond(0, 1, SINGLE_BOND);
  cerr << "After recreating nfrag " << mol.number_fragments() << '\n';
  cerr << "Same frag " << (mol.fragment_membership(0) == mol.fragment_membership(1)) << '\n';
}

void
DemoRing() {
  Molecule mol;
  // Adjacent three and four membered rings.
  mol.build_from_smiles("C1CCC1C1CC1");

  // We should get two rings.
  cerr << "Molecule contains " << mol.nrings() << " rings\n";

  // As listed here, the ring will have no status for aromaticity,
  // because aromaticity has not been determined.
  const int nr = mol.nrings();
  for (int i = 0; i < nr; ++i) {
    const Ring* ring = mol.ringi(i);
    cerr << " ring " << i << " contains " << ring->size() << " atoms\n";
    cerr << " atoms " << *ring << '\n';
  }

  // We can check whether atoms are in the same ring
  cerr << "same ring 0,1 " << mol.in_same_ring(0, 1)  << '\n';
  cerr << "same ring 3,4 " << mol.in_same_ring(3, 4)  << '\n';
  cerr << "same ring 4,5 " << mol.in_same_ring(4, 5)  << '\n';
}

void
DemoRingSystem() {
  Molecule mol;
  // Fused 3 and 4 membered rings.
  //    C ---- C
  //    |      |  \
  //    |      |   C
  //    |      |  /
  //    C ---- C
  // 
  mol.build_from_smiles("C12CCC1C2");

  // Contains 2 rings. 
  cerr << "Molecule contains " << mol.nrings() << " rings\n";

  // The two rings must have the same fused system identifier
  // Generally rings are sorted by size so the smallest rings should be first.
  const int nr = mol.nrings();
  for (int i = 0; i < nr; ++i) {
    const Ring* r = mol.ringi(i);
    cerr << *r << '\n';
    cerr << "fused rings:fused_system_identifier " << r->fused_system_identifier() << '\n';
    cerr << "is_fused " << r->is_fused() << '\n';
  }

  // For each ring, we can fetch the Ring's attached.
  // In this case, where will be two atoms in common between the two rings.
  for (int i = 0; i < nr; ++i) {
    const Ring* ri = mol.ringi(i);
    const int nfused = ri->fused_ring_neighbours();
    for (int j = 0; j < nfused; ++j) {
      const Ring* nbr = ri->fused_neighbour(j);
      cerr << "ring " << *ri << " fused to " << *nbr << '\n';
    }
  }


  // Two separated rings will have different fused_system_identifier

  mol.build_from_smiles("C1CC1C1CC1");
  for (int i = 0; i < 2; ++i) {
    const Ring* r = mol.ringi(i);
    cerr << "separated rings: fused_system_identifier " <<r->fused_system_identifier() << '\n';
    cerr << "is_fused " << r->is_fused() << '\n';
  }
}

void
DemoSubstuctureSearch() {
  Molecule mol;
  mol.build_from_smiles("OCN");

  Substructure_Query qry;
  qry.create_from_smarts("OCN");

  Substructure_Results sresults;
  const int nhits = qry.substructure_search(mol, sresults);
  cerr << "Query makes " << nhits << " matches\n";


}

int
Main(int argc, char** argv) {
  DemoEmptyMolecule();
  cerr << '\n';
  DemoCannotParseBadSmiles();
  cerr << '\n';
  DemoAtomicNumbers();
  cerr << '\n';
  DemoNcon();
  cerr << '\n';
  DemoAtomIteration();
  cerr << '\n';
  DemoFragments();
  cerr << '\n';
  DemoRing();
  cerr << '\n';
  DemoRingSystem();
  cerr << '\n';
  DemoSubstuctureSearch();
  cerr << '\n';

  return 0;
}

}  // namespace lillymol_introcution

int
main(int argc, char** argv) {
  int rc = lillymol_introcution::Main(argc, argv);

  return rc;
}
