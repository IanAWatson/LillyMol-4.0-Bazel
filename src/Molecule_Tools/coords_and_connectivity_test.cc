// Tests for the CachedMolecule class

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Molecule_Tools/coords_and_connectivity.h"

namespace {

using cached_molecule::CachedMolecule;

TEST(TestCachedMolecule, BadSmiles) {
  CachedMolecule mol;
  EXPECT_EQ(mol.ParseSmiles("invalid smiles"), std::nullopt);
}

// Test that ParseSmiles returns the number of atoms in the smiles.
struct SmilesAtomCount {
  std::string smiles;
  int natoms;
};

class TestCachedMoleculeNatoms : public testing::TestWithParam<SmilesAtomCount> {
  protected:
    CachedMolecule _cached_molecule;
};

TEST_P(TestCachedMoleculeNatoms, TestCachedMoleculeNatoms) {
  const auto params = GetParam();
  EXPECT_THAT(_cached_molecule.ParseSmiles(params.smiles), testing::Optional(params.natoms));
}
INSTANTIATE_TEST_SUITE_P(TestCachedMoleculeNatoms, TestCachedMoleculeNatoms, testing::Values(
  SmilesAtomCount{"C", 1},
  SmilesAtomCount{"CC", 2},
  SmilesAtomCount{"CCC", 3},
  SmilesAtomCount{"CCCC", 4},
  SmilesAtomCount{"CC(C)(C)C", 5},
  SmilesAtomCount{"c1ccccc1", 6}
));

}  // namespace
