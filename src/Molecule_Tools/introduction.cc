// An introduction to LillyMol.

// A Molecule object contains a vector of Atoms, and these atoms
// can be iterated via a C++ range for loop. But almost all
// times an Atom is passed from the Molecule it will be const.
// The reason for this is that if certain atomic properties are
// altered, then certain molecular properties would need to be
// either changed or invalidated. For example, if an atom were
// to be assigned a new isotopic value, the smiles for the molecule
// must be recomputed.
// In reality, Molecule is lazy, and will only recompute the smiles
// if requested. Changing the isotope just invalidates the smiles.
// Similariy, removing a bond between atoms would need to invalidate
// fragment membership, ring membership, symmetry, aromaticity...

#include <iostream>

#include "google/protobuf/text_format.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/output.h"
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
  cerr << "Is the molecule empty " << mol.empty() << '\n';
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

// The each_atom method requires a struct that responds to operator()(const Atom&);
void
DemoEachAtom() {
  Molecule mol;
  mol.build_from_smiles("OCN");

  struct SumNcon {
    int sum = 0;
    void operator()(const Atom& a) {
      sum += a.atomic_number();
    }
  };

  SumNcon sum_ncon;
  mol.each_atom(sum_ncon);
  cerr << "Sum of ncon " << sum_ncon.sum << '\n';
}

void
DemoEachAtomLambda() {
  Molecule mol;
  mol.build_from_smiles("OCN");

  int max_atomic_number = 0;
  mol.each_atom_lambda([&max_atomic_number] (const Atom& atom) {
    max_atomic_number = std::max(max_atomic_number, atom.atomic_number());
  });
  cerr << "highest atomic number " << max_atomic_number << '\n';
}

void
DemoEachBond() {
  Molecule mol;
  mol.build_from_smiles("CC=CC");

  struct CountSingleBonds {
    int single_bonds = 0;
    void operator() (const Bond& b) {
      if (b.is_single_bond()) {
        ++single_bonds;
      }
    }
  };

  CountSingleBonds count_single_bonds;

  mol.each_bond(count_single_bonds);
  cerr << "Find " <<count_single_bonds.single_bonds << " single bonds\n";

  // This is the same as.
  int single_bonds = 0;
  for (const Bond* bond : mol.bond_list()) {
    if (bond->is_single_bond()) {
      ++single_bonds;
    }
  }
}

void
DemoEachRing() {
  Molecule mol;
  mol.build_from_smiles("C1CC1C1CCC1C1CCCC1");

  struct LargestRingSize {
    int largest_ring = 0;
    void operator() (const Ring& r) {
      largest_ring = std::max(r.number_elements(), largest_ring);
    }
  };

  LargestRingSize largest_ring;

  mol.each_ring(largest_ring);
  cerr << "Find " << largest_ring.largest_ring << " largest ring\n";

  // This is the same as.
  int big_ring = 0;
  const int nrings = mol.nrings();
  for (int i = 0; i < nrings; ++i) {
    const Ring* ri = mol.ringi(i);
    big_ring = std::max(ri->number_elements(), big_ring);
  }
}

// Apply a member function to each atom.
// This particular demo is quite unusual, normally some
// kind of numeric accumulation would be done. THis
// works because IWString has a += operator.
void
DemoEachIndex() {
  Molecule mol;
  mol.build_from_smiles("OCN");

  Molecule::const_member_fn<const IWString&> fn = &Molecule::atomic_symbol;
  IWString result = mol.each_index<IWString>(fn);
  cerr << "Concatenated symbols '" << result << "'\n";
}

void
Demoisotopes() {
  Molecule mol;
  mol.build_from_smiles("CC");
  mol.set_isotope(0, 1);
  mol.set_isotope(1, 2);

  // Isotopes can be retrieved by querying either the molecule or the atom.
  cerr << "Isotope on atom 0 " << mol.isotope(0) << '\n';
  for (const Atom* atom : mol) {
    cerr << " isotope on atom " << atom->isotope() << '\n';
  }
}

void
DemoFormalCharges() {
  Molecule mol;
  mol.build_from_smiles("CC[N+H3]");
  cerr << "Does the molecule have formal charges " << mol.has_formal_charges() <<
          " net " << mol.net_formal_charge() << '\n';
  cerr << "Charge " << mol.formal_charge(2) << '\n';
  cerr << "Charge on atom " << mol.atom(2).formal_charge() << '\n';

  mol.set_formal_charge(2, 0);
  cerr << "After reset " << mol.formal_charge(2) << '\n';
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

// It is possible to iterate through the connections to an atom
// by iterating through the Molecule. More efficient and elegant
// to iterate via the Atom.
void
DemoBondIterationMolecule() {
  Molecule mol;
  mol.build_from_smiles("CO");
  int ncon = mol.ncon(0);
  for (int i = 0; i < ncon; ++i) {
    atom_number_t j = mol.other(0, i);
    cerr << "Atom 0 is bonded to " << j << '\n';
  }
}

// While it is possible to get the connections to an
// atom, this is usually less efficient than iterating,
// since iterating does not require the creation of a
// Set_of_Atoms object.
void
DemoConnections() {
  Molecule mol;
  mol.build_from_smiles("C(C)O");
  const Set_of_Atoms connected = mol.connections(0);
  cerr  << " atoms connected to 0 " << connected << '\n';
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

  // We can ask the molecule about the ring membership of an atom.
  // Note that Atom objects do not know anything about their ring status.

  const int matoms = mol.natoms();
  for (int i = 0; i < matoms; ++i) {
    cerr << " Atom " << i << ' ' << mol.smarts_equivalent_for_atom(i) <<
         " in " << mol.nrings(i) << " rings, ring bond count " << mol.ring_bond_count(i) << '\n';
  }
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
DemoAreBonded() {
  Molecule mol;
  mol.build_from_smiles("CC.C");

  const int matoms = mol.natoms();
  for (int i = 0; i < matoms; ++i) {
    for (int j = i + 1; j < matoms; ++j) {
      if (mol.are_bonded(i, j)) {
        cerr << "Atoms " << i << " and " << j << " bonded\n";
      } else {
        cerr << "Atoms " << i << " and " << j << " not bonded\n";
      }
    }
  }
}

void
DemoCoordinatesInAtoms() {
  Molecule mol;
  mol.build_from_smiles("C{{-1,1,0}}C{{0,0,0}}C{{1,0,0}}C{{2,1,1}}");
  cerr << "Distance between 0 1 " << mol.distance_between_atoms(0, 1) << '\n';
  cerr << "Bond angle 0 1 2 " << (mol.bond_angle(0, 1, 2) * RAD2DEG) << '\n';
  cerr << "Dihedral angle 0 1 2 3 " << (mol.dihedral_angle(0, 1, 2, 3) * RAD2DEG) << '\n';

  // Note that the distances above are geometric distances.
  // Topological distances are frequently of interest.
  const int matoms = mol.natoms();
  for (int i = 0; i < matoms; ++i) {
    for (int j = i + 1; j < matoms; ++j) {
      cerr << "Bonds between " << i << " and " << j << ' ' << mol.bonds_between(i, j) << '\n';
    }
  }

  // We can ask the molecule or an individual atom for their coordinates.
  // Internally, the Molecule fetches the Atom.
  cerr << "Atom 0 at " << mol.x(0) << ',' << mol.y(0) << ',' << mol.z(0) << '\n';
  const Atom& atom = mol.atom(0);
  cerr << "Atom 0 at " << atom.x() << ',' << atom.y() << ',' << atom.z() << '\n';
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

  // Loop over the hits - one in this case
  for (int i = 0; i < nhits; ++i) {
    const Set_of_Atoms* e = sresults.embedding(i);
    cerr << " atoms " << *e << '\n';
  }
}

// Using a Molecule_to_Match as the target for the query
// can be more efficient because computed properties are cached.
void
DemoSubstuctureSearchWithTarget() {
  Molecule mol;
  mol.build_from_smiles("Oc1ccccc1");

  Molecule_to_Match target(&mol);

  Substructure_Query qry;
  qry.create_from_smarts("O-c:c");

  Substructure_Results sresults;

  int nhits = qry.substructure_search(target, sresults);
  cerr << "Got " << nhits << " hits to phenol (unconstrained)\n";

  // If only one embedding per start atom, just one match will be found.
  qry.set_find_one_embedding_per_atom(1);
  nhits = qry.substructure_search(target, sresults);
  cerr << "Got " << nhits << " hits to phenol (one embedding per start atom)\n";

  qry.set_find_one_embedding_per_atom(0);

  // Turn off perception of symmetry equivalent matches.
  qry.set_perceive_symmetry_equivalent_matches(0);
  nhits = qry.substructure_search(target, sresults);
  cerr << "Got " << nhits << " hits to phenol (no symmetry)\n";

  qry.set_perceive_symmetry_equivalent_matches(1);
  nhits = qry.substructure_search(target, sresults);
  cerr << "Got " << nhits << " hits to phenol (back to default)\n";
}

void
DemoQueryFromProto() {
  Molecule mol;
  mol.build_from_smiles("Oc1ccccc1");

  const std::string s = R"pb(
query {
    min_nrings: 1
    ring_specifier {
      aromatic: true
      fused: 0
      base {
        heteroatom_count: 0
      }
    }
    smarts: "O-c:c"
}
)pb";

  // Build a proto from the text.
  SubstructureSearch::SubstructureQuery proto;  
  if (! google::protobuf::TextFormat::ParseFromString(s, &proto)) {
    cerr << "Cannot parse proto\n";
    return;
  }

  // Now build a substructure query from the proto.

  Substructure_Query qry;
  qry.ConstructFromProto(proto);

  Substructure_Results sresults;
  qry.substructure_search(mol, sresults);
  for (const Set_of_Atoms* embedding : sresults.embeddings()) {
    cerr << "From proto matches atom " << *embedding << '\n';
  }
}

void
DemoOutput() {
  Molecule mol;
  mol.build_from_smiles("C");
  mol.set_name("methane");

  mol.write_molecule(std::cout, FILE_TYPE_SMI, "");
  // which in this case is the same as
  std::cout << mol.smiles() << ' ' << mol.name() << '\n';
}

int
Main(int argc, char** argv) {
  DemoEmptyMolecule();
  cerr << '\n';
  DemoCannotParseBadSmiles();
  cerr << '\n';
  DemoAtomicNumbers();
  cerr << '\n';
  Demoisotopes();
  cerr << '\n';
  DemoFormalCharges();
  cerr << '\n';
  DemoNcon();
  cerr << '\n';
  DemoBondIteration();
  cerr << '\n';
  DemoBondIterationMolecule();
  cerr << '\n';
  DemoConnections();
  cerr << '\n';
  DemoAtomIteration();
  cerr << '\n';
  DemoEachAtom();
  cerr << '\n';
  DemoEachAtomLambda();
  cerr << '\n';
  DemoEachBond();
  cerr << '\n';
  DemoEachRing();
  cerr << '\n';
  DemoEachIndex();
  cerr << '\n';
  DemoAreBonded();
  cerr << '\n';
  DemoCoordinatesInAtoms();
  cerr << '\n';
  DemoFragments();
  cerr << '\n';
  DemoRing();
  cerr << '\n';
  DemoRingSystem();
  cerr << '\n';
  DemoSubstuctureSearch();
  cerr << '\n';
  DemoSubstuctureSearchWithTarget();
  cerr << '\n';
  DemoQueryFromProto();
  cerr << '\n';
  DemoOutput();
  cerr << '\n';

  return 0;
}

}  // namespace lillymol_introcution

int
main(int argc, char** argv) {
  int rc = lillymol_introcution::Main(argc, argv);

  return rc;
}