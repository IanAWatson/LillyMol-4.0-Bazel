// Tests for the CachedMolecule class
#include <cmath>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Molecule_Tools/coords_and_connectivity.h"

namespace {

using cached_molecule::CachedMolecule;

using testing::ElementsAreArray;
using testing::FloatNear;
using testing::Pointwise;

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

TEST(TestCachedMolecule, TestConnectivityEmptyMolecule) {
  CachedMolecule mol;
  ASSERT_EQ(mol.ParseSmiles(""), 0);
  std::vector<int> connections;
  EXPECT_EQ(mol.ConnectionMatrix(connections), 0);
}

TEST(TestCachedMolecule, TestConnectivitySingleAtom) {
  CachedMolecule mol;
  ASSERT_EQ(mol.ParseSmiles("C"), 1);
  std::vector<int> connections;
  EXPECT_EQ(mol.ConnectionMatrix(connections), 1);
  EXPECT_THAT(connections, ElementsAreArray({0}));
}

TEST(TestCachedMolecule, TestConnectivityTwoFragments) {
  CachedMolecule mol;
  constexpr int kNatoms = 2;
  ASSERT_EQ(mol.ParseSmiles("C.C"), kNatoms);
  std::vector<int> connections;
  EXPECT_EQ(mol.ConnectionMatrix(connections), kNatoms);
  EXPECT_EQ(connections.size(), kNatoms * kNatoms);
  EXPECT_THAT(connections, ElementsAreArray({0, 0, 0, 0}));
}

TEST(TestCachedMolecule, TestConnectivityCyclopropane) {
  CachedMolecule mol;
  constexpr int kNatoms = 3;
  ASSERT_EQ(mol.ParseSmiles("C1CC1"), kNatoms);
  std::vector<int> connections;
  EXPECT_EQ(mol.ConnectionMatrix(connections), kNatoms);
  EXPECT_EQ(connections.size(), kNatoms * kNatoms);
  EXPECT_THAT(connections, ElementsAreArray({0, 1, 1, 1, 0, 1, 1, 1, 0}));
}

TEST(TestCachedMolecule, TestConnectivityDoubleBond) {
  CachedMolecule mol;
  constexpr int kNatoms = 3;
  ASSERT_EQ(mol.ParseSmiles("C1C=C1"), kNatoms);
  std::vector<int> connections;
  EXPECT_EQ(mol.ConnectionMatrix(connections), kNatoms);
  EXPECT_EQ(connections.size(), kNatoms * kNatoms);
  EXPECT_THAT(connections, ElementsAreArray({0, 1, 1, 1, 0, 2, 1, 2, 0}));
}

TEST(TestCachedMolecule, TestConnectivityTripleBond) {
  CachedMolecule mol;
  constexpr int kNatoms = 3;
  ASSERT_EQ(mol.ParseSmiles("C1C#C1"), kNatoms);
  std::vector<int> connections;
  EXPECT_EQ(mol.ConnectionMatrix(connections), kNatoms);
  EXPECT_EQ(connections.size(), kNatoms * kNatoms);
  EXPECT_THAT(connections, ElementsAreArray({0, 1, 1, 1, 0, 3, 1, 3, 0}));
}

TEST(TestCachedMolecule, TestCoordsEmptyMolecule) {
  CachedMolecule mol;
  constexpr int kNatoms = 0;
  ASSERT_EQ(mol.ParseSmiles(""), kNatoms);
  std::vector<float> coords;
  EXPECT_EQ(mol.GetCoords(coords), kNatoms);
  EXPECT_EQ(coords.size(), kNatoms * kNatoms);
}

TEST(TestCachedMolecule, TestCoordsNoCoordinates) {
  CachedMolecule mol;
  constexpr int kNatoms = 1;
  ASSERT_EQ(mol.ParseSmiles("C"), kNatoms);
  std::vector<float> coords;
  EXPECT_EQ(mol.GetCoords(coords), kNatoms);
  EXPECT_EQ(coords.size(), kNatoms*3);
  EXPECT_THAT(coords, ElementsAreArray({0, 0, 0}));
}

TEST(TestCachedMolecule, TestCoordsSingleAtomWithCoords) {
  CachedMolecule mol;
  constexpr int kNatoms = 1;
  ASSERT_EQ(mol.ParseSmiles("C{{1,2,-3}}"), kNatoms);
  std::vector<float> coords;
  EXPECT_EQ(mol.GetCoords(coords), kNatoms);
  EXPECT_EQ(coords.size(), kNatoms*3);
  const std::vector<float> expected{1, 2, -3};
  EXPECT_THAT(coords, Pointwise(FloatNear(1.0e-05), expected));
}

// LillyMol smiles can have coordinates appended to each atom symbol
// surrounded by matching {{ and }}. No spaces are allowed.
TEST(TestCachedMolecule, TestCoordsTwoAtomsWithCoords) {
  CachedMolecule mol;
  constexpr int kNatoms = 2;
  ASSERT_EQ(mol.ParseSmiles("C{{1,2,-3}}C{{-4,-2,1}}"), kNatoms);
  std::vector<float> coords;
  EXPECT_EQ(mol.GetCoords(coords), kNatoms);
  EXPECT_EQ(coords.size(), kNatoms*3);
  const std::vector<float> expected{1, 2, -3, -4, -2, 1};
  EXPECT_THAT(coords, Pointwise(FloatNear(1.0e-05), expected));
}

TEST(TestCachedMolecule, TestCoordsComplexMoleculeRandomCoordinates) {
  CachedMolecule mol;
  constexpr int kNatoms = 6;
  ASSERT_EQ(mol.ParseSmiles("CCCCCC"), kNatoms);
  std::vector<float> expected;
  expected.reserve(kNatoms * 3);
  for (int i = 0; i < kNatoms; ++i) {
    float ifloat = static_cast<float>(i);
    float x = sin(ifloat);
    float y = tan(ifloat);
    float z = sqrt(ifloat);
    mol.molecule().setxyz(i, x, y, z);
    expected.push_back(x);
    expected.push_back(y);
    expected.push_back(z);
  }
  std::vector<float> coords;
  EXPECT_EQ(mol.GetCoords(coords), kNatoms);
  EXPECT_EQ(coords.size(), kNatoms*3);
  EXPECT_THAT(coords, Pointwise(FloatNear(1.0e-05), expected));
}

TEST(TestCachedMolecule, TestSubstructureQuery) {
  CachedMolecule mol;
  constexpr int kNatoms = 4;
  ASSERT_EQ(mol.ParseSmiles("NCOC"), kNatoms);
}


}  // namespace
