// Tests for generate_atom_separation

#include <iostream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "Molecule_Tools/atom_separations.h"

namespace {
using separated_atoms::Hasher;

TEST(TestHasher, TestTwo) {
  Hasher hasher;
  EXPECT_EQ(hasher.Value(6, 6, 1), 1011);
  EXPECT_EQ(hasher.Value(7, 7, 1), 2012);
  EXPECT_EQ(hasher.Value(8, 8, 1), 3013);
  EXPECT_EQ(hasher.Value(15, 15, 1), 4014);
  EXPECT_EQ(hasher.Value(16, 16, 1), 5015);
}

TEST(TestBondSymmetry, Test1) {
  const IWString smiles("CCc1cc(CC)cc(CC)c1");
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(smiles));
  std::unique_ptr<int[]> bond_sym = separated_atoms::BondSymmetry(m);
  for (int i = 0; i < m.natoms(); ++i) {
    if (i < 3) {
      EXPECT_EQ(bond_sym[i], 1);
    } else {
      EXPECT_EQ(bond_sym[i], 0);
    }
  }
}

}  // namespace

