#include <algorithm>

#include "googletest/include/gtest/gtest.h"

#include "aromatic.h"
#include "iwstandard.h"
#include "molecule.h"

namespace {


class TestStandardisation : public testing::Test
{
  protected:
    void SetUp() override;

  protected:
    Chemical_Standardisation _chemical_standardisation;

    IWString _smiles;
    Molecule _m1;
    Molecule _m2;
};

void
TestStandardisation::SetUp()
{
  set_global_aromaticity_type(Daylight);
}

TEST_F(TestStandardisation, EmptyMolecule)
{
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.activate_all();
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

TEST_F(TestStandardisation, TestAcidYes)
{
  _smiles = "CC(=O)[O-]";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  EXPECT_EQ(_m1.molecular_formula(), "C2O2H3");
  EXPECT_EQ(_m1.smiles(), "CC(=O)[O-]");

  _chemical_standardisation.Activate(CS_ACID, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);
  EXPECT_EQ(_m1.molecular_formula(), "C2O2H4");
  EXPECT_EQ(_m1.smiles(), "CC(=O)O");
  EXPECT_EQ(_m1.unique_smiles(), "OC(=O)C");
}

TEST_F(TestStandardisation, TestChargedImidazole)
{
  _smiles = "CN1C=C[N+](CC)=C1";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  EXPECT_EQ(_m1.molecular_formula(), "C6N2H11");
  EXPECT_EQ(_m1.smiles(), "CN1C=C[N+](=C1)CC");

  _chemical_standardisation.Activate(CS_CHARGED_IMIDAZOLE, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);

  _smiles = "CCN1C=C[N+](C)=C1";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_m1.molecular_formula(), "C6N2H11");
  EXPECT_EQ(_m1.unique_smiles(), "C[n+]1c[n](CC)cc1");

  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);
  EXPECT_EQ(_m1.unique_smiles(), "C[n]1c[n+](CC)cc1");

  _m1.invalidate_smiles();

  // Some random runs to make sure things do not crash.

  for (int i = 0; i < 5; ++i) {
    const IWString & smiles = _m1.random_smiles();
    Molecule m;
    ASSERT_TRUE(m.build_from_smiles(smiles));
    _chemical_standardisation.process(m);
    EXPECT_EQ(m.unique_smiles(), "C[n]1c[n+](CC)cc1");
  }
}

TEST_F(TestStandardisation, TestMisdrawnSulfonamideNoChange)
{
  _smiles = "CCS(=O)(=O)NC";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_MSDSA, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
}

TEST_F(TestStandardisation, TestMisdrawnSulfonamideChanges)
{
  _smiles = "CCS(=O)(O)=NC";
  ASSERT_TRUE(_m1.build_from_smiles(_smiles));
  EXPECT_EQ(_chemical_standardisation.process(_m1), 0);
  _chemical_standardisation.Activate(CS_MSDSA, /*verbose*/ false);
  EXPECT_EQ(_chemical_standardisation.process(_m1), 1);
  EXPECT_EQ(_m1.unique_smiles(), "CNS(=O)(=O)CC");
}

}  // namespace
