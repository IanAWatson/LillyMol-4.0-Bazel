// Tests for ring_closure_enumeration

#include <array>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "Molecule_Tools/ring_closure_enumeration_lib.h"

namespace {
using ring_closure_enumeration::BondFormation;
using ring_closure_enumeration::Reagents;

TEST(TestBondFormation, NoComponent) {
  BondFormation bond_formation;
  EXPECT_FALSE(bond_formation.Build(""));
}

TEST(TestBondFormation, ComponentNothingElse) {
  BondFormation bond_formation;
  EXPECT_FALSE(bond_formation.Build("3"));
}

TEST(TestBondFormation, ComponentColonOnly) {
  BondFormation bond_formation;
  EXPECT_FALSE(bond_formation.Build("3:"));
}

TEST(TestBondFormation, MissingIsotope) {
  BondFormation bond_formation;
  EXPECT_FALSE(bond_formation.Build("1:-3"));
}

TEST(TestBondFormation, BadIsotope1) {
  BondFormation bond_formation;
  EXPECT_FALSE(bond_formation.Build("1:foo-3"));
}

TEST(TestBondFormation, ZeroIsotope) {
  BondFormation bond_formation;
  EXPECT_FALSE(bond_formation.Build("1:0-3"));
}

TEST(TestBondFormation, MissingRingClosure) {
  BondFormation bond_formation;
  EXPECT_FALSE(bond_formation.Build("1:1-"));
}

TEST(TestBondFormation, BadRingClosure) {
  BondFormation bond_formation;
  EXPECT_FALSE(bond_formation.Build("1:1-foo"));
}

TEST(TestBondFormation, ZeroRingClosure) {
  BondFormation bond_formation;
  EXPECT_FALSE(bond_formation.Build("1:1-0"));
}

TEST(TestBondFormation, BuildSingle) {
  BondFormation bond_formation;
  ASSERT_TRUE(bond_formation.Build("1:91-92"));
  EXPECT_EQ(bond_formation.isotope(), 91);
  EXPECT_EQ(bond_formation.btype(), '-');
  EXPECT_EQ(bond_formation.ring_closure(), 92);
}

TEST(TestBondFormation, BuildDouble) {
  BondFormation bond_formation;
  ASSERT_TRUE(bond_formation.Build("1:91=92"));
  EXPECT_EQ(bond_formation.isotope(), 91);
  EXPECT_EQ(bond_formation.btype(), '=');
  EXPECT_EQ(bond_formation.ring_closure(), 92);
}

TEST(TestBondFormation, BuildTriple) {
  BondFormation bond_formation;
  ASSERT_TRUE(bond_formation.Build("1:91#92"));
  EXPECT_EQ(bond_formation.isotope(), 91);
  EXPECT_EQ(bond_formation.btype(), '#');
  EXPECT_EQ(bond_formation.ring_closure(), 92);
}

TEST(TestReagents, TestSingle) {
  Reagents reagents;
  Molecule* m = new Molecule();
  ASSERT_TRUE(m->build_from_smiles("[3CH3]-N"));
  reagents.Add(m);
  BondFormation bond_formation;
  ASSERT_TRUE(bond_formation.Build("1:3-91"));
  EXPECT_TRUE(reagents.CreateFragments(bond_formation));
  const auto& smiles = reagents.smiles();
  EXPECT_EQ(smiles.size(), 1);
  EXPECT_STREQ(*smiles[0], "[3CH3]-%91N");
}

TEST(TestReagents, TestDouble) {
  Reagents reagents;
  Molecule* m = new Molecule();
  ASSERT_TRUE(m->build_from_smiles("[3CH3]-N"));
  reagents.Add(m);
  BondFormation bond_formation;
  ASSERT_TRUE(bond_formation.Build("1:3=91"));
  EXPECT_TRUE(reagents.CreateFragments(bond_formation));
  const auto& smiles = reagents.smiles();
  EXPECT_EQ(smiles.size(), 1);
  EXPECT_STREQ(*smiles[0], "[3CH3]=%91N");
}

TEST(TestReagents, TestTriple) {
  Reagents reagents;
  Molecule* m = new Molecule();
  ASSERT_TRUE(m->build_from_smiles("[3CH3]-N"));
  reagents.Add(m);
  BondFormation bond_formation;
  ASSERT_TRUE(bond_formation.Build("1:3#91"));
  EXPECT_TRUE(reagents.CreateFragments(bond_formation));
  const auto& smiles = reagents.smiles();
  EXPECT_EQ(smiles.size(), 1);
  EXPECT_STREQ(*smiles[0], "[3CH3]#%91N");
}

TEST(TestReagents, TestMultiple) {
  Reagents reagents;
  Molecule* m = new Molecule();
  ASSERT_TRUE(m->build_from_smiles("[3CH3]-NCC[4CH2]O"));
  reagents.Add(m);

  std::array<BondFormation, 2> bond_formations;

  ASSERT_TRUE(bond_formations[0].Build("1:3#91"));
  ASSERT_TRUE(bond_formations[1].Build("1:4=92"));

  resizable_array<BondFormation*> ptrs;
  ptrs << bond_formations.data();
  ptrs << bond_formations.data() + 1;

  EXPECT_TRUE(reagents.CreateFragments(ptrs));
  const auto& smiles = reagents.smiles();
  EXPECT_EQ(smiles.size(), 1);
  EXPECT_STREQ(*smiles[0], "[3CH3]#%91NCC[4CH2]=%92O");
}

}  // namespace
