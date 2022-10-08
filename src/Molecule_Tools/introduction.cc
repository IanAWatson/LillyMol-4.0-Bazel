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

// Return a newly created molecule build from `smiles`.
// The molecule name is not set.
// If smiles interpretation failes, an empty molecule is returned.
// A better approach would be to use std::optional.
Molecule
MolFromSmiles(const char* smiles) {
  Molecule result;
  if (!result.build_from_smiles(smiles)) {
    cerr << "Invalid smiles '" << smiles << "'\n";
    return Molecule();
  }

  return result;
}

// In LillyMol, there is a periodic table datastructure which contains
// an Element object for all elements known to the system. We can
// fetch pointers to those elements via various means.

// For historical reasons, to get an element from a two letter element
// we must use get_element_from_symbol_no_case_conversion. This should
// be fixed. We had some early mdl files where the symbols were all
// uppercase.
void
DemoElements() {
  const Element* carbon = get_element_from_atomic_number(6);
  const Element* uranium = get_element_from_symbol_no_case_conversion("U");
  const Element* chlorine = get_element_from_symbol_no_case_conversion("Cl");
  const Element* iron = get_element_from_atomic_number(26);

  cerr << "Carbon is " << carbon->symbol() << " atnum " << carbon->atomic_number() << '\n';
  cerr << "normal isotope " << carbon->normal_isotope() << '\n';
  cerr << "atomic_mass " << carbon->atomic_mass() << '\n';
  cerr << "exact_mass " << carbon->exact_mass() << '\n';
  cerr << "normal_valence " << carbon->normal_valence() << '\n';
  cerr << "number_alternate_valences " << carbon->number_alternate_valences() << '\n';
  cerr << "carbon organic " << carbon->organic() << '\n';
  cerr << "uranium organic " << uranium->organic() << '\n';
  cerr << "carbon needs_square_brackets " << carbon->needs_square_brackets() << '\n';
  cerr << "uranium needs_square_brackets " << uranium->needs_square_brackets() << '\n';
  cerr << "carbon is_halogen " << carbon->is_halogen() << '\n';
  cerr << "uranium is_halogen " << uranium->is_halogen() << '\n';
  cerr << "chlorine is_halogen " << chlorine->is_halogen() << '\n';
  cerr << "carbon is_in_periodic_table " << carbon->is_in_periodic_table() << '\n';
  cerr << "uranium is_in_periodic_table " << uranium->is_in_periodic_table() << '\n';
  cerr << "carbon is_metal " << carbon->is_metal() << '\n';
  cerr << "iron is_metal " << iron->is_metal() << '\n';

  // Computing the number of pi electrons and lone pairs from an element is
  // complex and you are better off asking a Molecule for the pi electron
  // or lone pair count on a given atom.
}

void
DemoAtomIteration() {
  Molecule mol = MolFromSmiles("CCC");

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
  Molecule mol = MolFromSmiles("OCN");

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

// Find the highest atomic number in the molecule.
void
DemoEachAtomLambda() {
  Molecule mol = MolFromSmiles("OCN");

  int max_atomic_number = 0;
  mol.each_atom_lambda([&max_atomic_number] (const Atom& atom) {
    max_atomic_number = std::max(max_atomic_number, atom.atomic_number());
  });
  cerr << "highest atomic number " << max_atomic_number << '\n';
}

// Iterate through all bonds(edges) in the molecule.
void
DemoEachBond() {
  Molecule mol = MolFromSmiles("CC=CC");

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

// Iterators through rings.
void
DemoEachRing() {
  Molecule mol = MolFromSmiles("C1CC1C1CCC1C1CCCC1");

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

// Apply a Molecule member function to each atom.
// This particular demo is quite unusual, normally some
// kind of numeric accumulation would be done. THis
// works because IWString has a += operator.
void
DemoEachIndex() {
  Molecule mol = MolFromSmiles("OCN");

  Molecule::const_member_fn<const IWString&> fn = &Molecule::atomic_symbol;
  IWString result = mol.each_index<IWString>(fn);
  cerr << "Concatenated symbols '" << result << "'\n";
}

void
DemoIsotopes() {
  Molecule mol = MolFromSmiles("CC");
  mol.set_isotope(0, 1);
  mol.set_isotope(1, 2);

  // Isotopes can be retrieved by querying either the molecule or the atom.
  cerr << "Isotope on atom 0 " << mol.isotope(0) << '\n';
  for (const Atom* atom : mol) {
    cerr << " isotope on atom " << atom->isotope() << '\n';
  }

  // Isotopes are arbitrary numbers.
  mol.set_isotope(0, 999);
  mol.set_isotope(1, std::numeric_limits<int>::max());
  cerr << "Ridiculous isotopes " << mol.smiles() << '\n';

  // We should be able to build a new molecule from `mol`
  // and their unique smiles should be the same.
  Molecule round_trip;
  round_trip.build_from_smiles(mol.smiles());
  cerr << "Same? " << (mol.unique_smiles() == round_trip.unique_smiles()) << '\n';
}

void
DemoIsotopesMolecularWeight() {
  Molecule mol = MolFromSmiles("OCN");
  float amw = mol.molecular_weight();
  cerr << "Non isotopic molecular weight " << amw << '\n';
  mol.set_isotope(0, 1);
  amw = mol.molecular_weight();
  cerr << "With isotopes, molecular weight is zero " << amw << '\n';

  amw = mol.molecular_weight_ignore_isotopes();
  cerr << "If isotopes are ignored, amw " << amw << '\n';

  // Isotope 1 gets counted as contributing 1.0 to the amw.
  amw = mol.molecular_weight_count_isotopes();
  cerr << "If isotopes included " << amw << '\n';
}

void
DemoImplicitAndExplicitHydrogens() {
  Molecule mol = MolFromSmiles("OCN");

  // Note that the ncon method, the number of connections to
  // an atom, does NOT reflect implicit Hydrogens.
  int matoms = mol.natoms();
  cerr << "Molecule begins with " << matoms << " atoms\n";
  for (int i = 0; i < matoms; ++i) {
    cerr << "Atom " << i << " ncon " << mol.ncon(i) << " connections has " << mol.hcount(i) << " Hydrogens\n";
  }

  // Convert implicit Hydrogens to explicit.
  // Note that the hcount for the existing atoms, natoms, should
  // be unchaged, because the hcount() method includes both
  // explicit and implicit Hydrogens.
  // but now that the formerly implicit Hydrogens are now
  // explicit members of the connection table, the ncon values
  // for the existing atoms will change.

  cerr << "After implicit Hydrogens become explicit\n";
  mol.make_implicit_hydrogens_explicit();
  for (int i = 0; i < matoms; ++i) {
    cerr << "Atom " << i << " ncon " << mol.ncon(i) << " connections has " << mol.hcount(i) << " Hydrogens\n";
  }

  cerr << "Molecule now contains " << mol.natoms() << " atoms\n";
}

// One of the most troublesome parts of LillyMol is the
// implicit hydrogens known attribute. If a smiles begins with
// C-[CH]-C where the centre Carbon is deficient one Hydrogen
// what happens when we add a Carbon atom to it, thereby
// satisfying the valence.
void 
DemoImplicitHydrogensKnown() {
  Molecule mol = MolFromSmiles("C-[CH]-C");
  cerr << "Valence? " << mol.valence_ok() << '\n';

  // Indeed the problem is with atom 1.
  for (int i = 0; i < mol.natoms(); ++i) {
    cerr << " atom " << i << " valence " << mol.valence_ok(i) << '\n';
  }

  // Add a carbon atom.
  // First fetch a pointer to the Element with atomic number 6.
  // Then call mol.add() with a pointer to the Element.
  // Internally, mol will instantiate an Atom which has
  // that element.
  const Element* carbon = get_element_from_atomic_number(6);
  mol.add(carbon);

  // The molecule now has two fragments.
  cerr << "Molecule has " << mol.number_fragments() << " fragments\n";

  // Add a bond fromthis last atom added to the middle atom to satisfy the valence.
  // The newly added atom will always have the highest atom number.
  mol.add_bond(1, 3, SINGLE_BOND);
  cerr << "Now " << mol.number_fragments() << " fragments\n";

  // Is the valence ok now?
  cerr << "After atom addition valence? " << mol.valence_ok() << '\n';

  // But the smiles shows that the centre atom still has its
  // implicit Hydrogens known attribute set. Even though it
  // is no longer necessary.
  cerr << "Smiles now " << mol.smiles() << '\n';
  cerr << "Is the value fixed " << mol.implicit_hydrogens_known(1) << '\n';

  // Remove the implicit Hydrogen known attribute.
  // Now the smiles does not have the square brackets.
  mol.set_implicit_hydrogens_known(1, 0);
  cerr << "Smiles now " << mol.smiles() << '\n';
  cerr << "Is the value fixed " << mol.implicit_hydrogens_known(1) << '\n';
}

// Random smiles are useful for testing. Algorithms should not
// depend on the order of the smiles.
void
DemoRandomSmiles() {
  // Caffeine
  Molecule mol = MolFromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C");

  // Note that we make a copy of the unique smiles.
  // Generally we get a reference to the unique smiles,
  //   IWString& starting_unique_smiles = mol.unique_smiles()
  // But in this case the smiles of `mol` will be destroyed when the
  // random smiles are formed below, and the reference would
  // be invalid.

  const IWString starting_unique_smiles = mol.unique_smiles();

  for (int i = 0; i < 10; ++i) {
    cerr << " random smiles " << mol.random_smiles() << '\n';
    Molecule hopefully_the_same;
    hopefully_the_same.build_from_smiles(mol.random_smiles());
    cerr << "Same? " << (starting_unique_smiles == hopefully_the_same.unique_smiles()) << '\n';
  }
}

void
DemoFormalCharges() {
  Molecule mol = MolFromSmiles("CC[N+H3]");

  cerr << "Does the molecule have formal charges " << mol.has_formal_charges() <<
          " net " << mol.net_formal_charge() << '\n';
  cerr << "Charge " << mol.formal_charge(2) << '\n';
  cerr << "Charge on atom " << mol.atom(2).formal_charge() << '\n';

  mol.set_formal_charge(2, 0);
  cerr << "After reset " << mol.formal_charge(2) << '\n';
}

void
DemoNcon() {
  Molecule mol = MolFromSmiles("CC(C)(C)C");

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
  Molecule mol = MolFromSmiles("CC(C)(C)C");

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
// by iterating through the Molecule. It is more efficient and elegant
// to iterate via the Atom rather than via the Molecule.
void
DemoBondIterationMolecule() {
  Molecule mol = MolFromSmiles("CO");

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
  Molecule mol = MolFromSmiles("C(C)O");

  const Set_of_Atoms connected = mol.connections(0);
  cerr  << " atoms connected to 0 " << connected << '\n';
}

// Start with ethane, verify one fragment, then remove the bond
// to generate two methanes.
void
DemoFragments() {
  Molecule mol = MolFromSmiles("CC");

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

// when we have a multi-fragment molecule, we can create
// molecules from each fragment.
void
DemoCreateComponents() {
  Molecule mol = MolFromSmiles("C.C.C.C.C.C.C.C.C");
  resizable_array_p<Molecule> fragments;
  mol.create_components(fragments);
  cerr << "Created " << fragments.size() << " fragments from " << mol.smiles() << '\n';

  // The fragments should all have the same smiles.
  // Note that we cannot use 'const Molecule* m' here because generating
  // a smiles is a potentially non-const operation - see lazy evaluation.
  for (Molecule* m : fragments) {
    if (m->unique_smiles() != "C") {
      cerr << "Invalid fragment unique smiles\n";
    }
  }
}

// Normally getting the largest fragment is straightforward and
// yields the expected result.
void
DemoReduceToLargestFragment() {
  Molecule mol = MolFromSmiles("CC.C");
  mol.reduce_to_largest_fragment();
  cerr << "after reduce_to_largest_fragment have " << mol.natoms() << " atoms " << mol.smiles() << '\n';
}

void
DemoRing() {
  // Adjacent three and four membered rings.
  Molecule mol = MolFromSmiles("C1CCC1C1CC1");

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
  // Fused 3 and 4 membered rings.
  //    C ---- C
  //    |      |  \
  //    |      |   C
  //    |      |  /
  //    C ---- C

  Molecule mol = MolFromSmiles("C12CCC1C2");

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

  // For each ring, we can fetch pointers to the Ring's attached.
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
  Molecule mol = MolFromSmiles("CC.C");

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

// A LillyMol smiles extension is to allow coordinates in smiles.
// {{x,y,z}} is appended after each atom.
// This is more compact than a .sdf representation, and allows
// molecules to be processed one per line.
void
DemoCoordinatesInAtoms() {
  Molecule mol = MolFromSmiles("C{{-1,1,0}}C{{0,0,0}}C{{1,0,0}}C{{2,1,1}}");

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

// Performing a substructure search requires instantiation of
// a Substructure_Query object. That is then constructed
// from smarts, a query file, a proto, or other means.
// Once the Substructure_Query is built, use it to perform
// substructure searches on Molecule's.
void
DemoSubstuctureSearch() {
  Molecule mol = MolFromSmiles("OCN");

  Substructure_Query qry;
  qry.create_from_smarts("OCN");

  // For historical reasons, this method takes a pointer.
  const int nhits = qry.substructure_search(&mol);
  cerr << "Query makes " << nhits << " matches\n";
}

// The example above did not save the matched atoms. If we need
// the matched atoms, pass a Substructure_Results to the search.
void
DemoSubstuctureSearchSaveResults() {
  Molecule mol = MolFromSmiles("OCN");

  Substructure_Query qry;
  qry.create_from_smarts("OCN");

  Substructure_Results sresults;
  const int nhits = qry.substructure_search(mol, sresults);

  // Loop over the hits - one in this case. There is an iterator
  // for shown below.
  for (int i = 0; i < nhits; ++i) {
    const Set_of_Atoms* e = sresults.embedding(i);
    cerr << " atoms " << *e << '\n';
  }
}

// If all we care about is whether or not a query matches
// and we do not care about how many matches there might 
// be, and we do not want to save the matched atoms, we
// can limit the scope of the search.
// The initial molecule/smarts generates 559872 query
// matches. Every 10,000 matches found, there will be
// a message to stderr informing you of a potentially
// catastrophic matching.
// 559872 is 2^8 * 3^7
void
DemoParsimoniousSearch() {
  const char* smiles = "C(C)(C)(C)c1c(C(C)(C)(C))c(C(C)(C)(C))c(C(C)(C)(C))c(C(C)(C)(C))c1C(C)(C)(C)";
  Molecule mol = MolFromSmiles(smiles);

  Substructure_Query qry;
  qry.create_from_smarts(smiles);

  Substructure_Results sresults;
  int nhits = qry.substructure_search(mol, sresults);
  cerr << "Without constraints, find " << nhits << " substructure matches\n";

  qry.set_max_matches_to_find(1);
  qry.set_save_matched_atoms(0);

  nhits = qry.substructure_search(mol, sresults);
  cerr << "With constraints find " << nhits << " substructure matches\n";
}

// Using a Molecule_to_Match as the target for the query
// can be more efficient because computed properties are cached.
// Definitely recommended if the same molecule is being searched
// by multiple queries.
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

// Set up a vector of Substructure_Query objects.
// A resizable_array_p is a vector that contains pointers
// to objects. The resizable_array_p manages ownership of
// the objects.
void
DemoMultipleQueries() {
  resizable_array_p<Substructure_Query> queries;
  constexpr int kNqueries = 10;

  // Build up a set of smarts that are C, CC, CCC ...
  IWString smarts = 'C';
  for (int i = 0; i < kNqueries; ++i) {
    std::unique_ptr<Substructure_Query> qry = std::make_unique<Substructure_Query>();
    qry->create_from_smarts(smarts);
    // The name of the query is the same as the smarts.
    qry->set_comment(smarts);
    queries << qry.release();

    smarts += 'C';
  }

  Molecule mol = MolFromSmiles("CCCCC");
  Molecule_to_Match target(&mol);

  // It is interesting to contemplate the number of matches found
  // to each of these queries.
  Substructure_Results sresults;
  for (Substructure_Query* qry : queries) {
    int nhits = qry->substructure_search(target, sresults);
    cerr << nhits << " matches to query " << qry->comment() << '\n';
  }
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

// If enabled, arbitrary two character combinations are valid elements.
void
DemoStrangeElements() {
  set_auto_create_new_elements(1);
  Molecule mol = MolFromSmiles("[Th][Eq][U]IC[K][Br]O[W]NFO[X][Ju][Mp]SO[Ve][Rt][La][Zy][D]O[G]");
  cerr << "Smiles " << mol.smiles() << '\n';

  cerr << "contains_non_periodic_table_elements " << mol.contains_non_periodic_table_elements() << '\n';
  cerr << "organic " << mol.organic_only() << '\n';

  Substructure_Query qry;
  qry.create_from_smarts("[Th][Eq][U]IC[K]");

  Substructure_Results sresults;
  int nhits = qry.substructure_search(mol, sresults);
  cerr << nhits << " hits to '[Th][Eq][U]'\n";

  // But if there are element names that might collide with smarts
  // directives, we can enclose that in #{}
  qry.create_from_smarts("[#{Rt}][La]");
  nhits = qry.substructure_search(mol, sresults);
  cerr << nhits << " hits to '[Rt][La]'\n";

  set_auto_create_new_elements(0);
}

// If enabled, any sequence of letters is an ok element
void
DemoAnyLengthElement() {
  set_auto_create_new_elements(1);
  set_atomic_symbols_can_have_arbitrary_length(1);

  Molecule mol = MolFromSmiles("[Ala]1[Arg][Asn][Asp]1");

  cerr << "Peptide " << mol.smiles() << '\n';

  // All regular smarts directives are available.
  Substructure_Query qry;
  qry.create_from_smarts("[#{Ala}D2R]1[#{Arg}D2x2][#{Asn}][#{Asp}D2]1");

  Substructure_Results sresults;
  const int nhits = qry.substructure_search(mol, sresults);
  cerr << nhits << " hits to Ala-Arg-Asn-Asp\n";

  set_atomic_symbols_can_have_arbitrary_length(0);
  set_auto_create_new_elements(0);
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
  DemoIsotopes();
  cerr << '\n';
  DemoIsotopesMolecularWeight();
  cerr << '\n';
  DemoImplicitAndExplicitHydrogens();
  cerr << '\n';
  DemoImplicitHydrogensKnown();
  cerr << '\n';
  DemoRandomSmiles();
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
  DemoElements();
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
  DemoCreateComponents();
  cerr << '\n';
  DemoReduceToLargestFragment();
  cerr << '\n';
  DemoRing();
  cerr << '\n';
  DemoRingSystem();
  cerr << '\n';
  DemoSubstuctureSearch();
  cerr << '\n';
  DemoSubstuctureSearchSaveResults();
  cerr << '\n';
  DemoParsimoniousSearch();
  cerr << '\n';
  DemoSubstuctureSearchWithTarget();
  cerr << '\n';
  DemoMultipleQueries();
  cerr << '\n';
  DemoQueryFromProto();
  cerr << '\n';
  DemoOutput();
  cerr << '\n';
  DemoStrangeElements();
  cerr << '\n';
  DemoAnyLengthElement();
  cerr << '\n';

  return 0;
}

}  // namespace lillymol_introcution

int
main(int argc, char** argv) {
  int rc = lillymol_introcution::Main(argc, argv);

  return rc;
}
