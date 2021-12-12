// Tests for aromaticity settings

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "aromatic.h"
#include "molecule.h"
#include "smiles.h"

namespace {

struct SmilesSettingsExpected {
  IWString smiles;
  int aromaticity;
  int min_aromatic_ring_size;
  int max_aromatic_ring_size;
  int two_electrons_aromatic;

  IWString expected;
};

class TestAromaticity: public testing::TestWithParam<SmilesSettingsExpected> {
  protected:
    Molecule _m;

    void SetUp();
};

void
TestAromaticity::SetUp() {
  reset_aromatic_file_scope_variables();
}

TEST_P(TestAromaticity, Various) {
  const auto params = GetParam();
  set_global_aromaticity_type(params.aromaticity);
  set_default_unique_smiles_aromaticity(params.aromaticity);
  set_min_aromatic_ring_size(params.min_aromatic_ring_size);
  set_max_aromatic_ring_size(params.max_aromatic_ring_size);
  set_allow_two_electron_systems_to_be_aromatic(params.two_electrons_aromatic);
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  EXPECT_EQ(_m.unique_smiles(), params.expected);
}
INSTANTIATE_TEST_SUITE_P(TestSettings, TestAromaticity, testing::Values(
  SmilesSettingsExpected{"O=C1C=CC1=O", Daylight, 4, 6, 0, "O=C1C(=O)C=C1"},
  SmilesSettingsExpected{"O=C1C=CC1=O", Daylight, 4, 6, 1, "O=c1c(=O)cc1"},
  SmilesSettingsExpected{"O=C1C=CC1=O", Daylight, 5, 6, 1, "O=C1C(=O)C=C1"},
  SmilesSettingsExpected{"N1=CC=CC=CC=C1", Daylight, 5, 6, 1, "N1=CC=CC=CC=C1"},
  SmilesSettingsExpected{"N1=CC=CC=CC=C1", Daylight, 5, 8, 1, "N1=CC=CC=CC=C1"},
  SmilesSettingsExpected{"N1=CC=CC=CC=C1", ANY_EVEN_NUMBER_OF_PI_ELECTRONS, 5, 9, 1, "[n]1ccccccc1"},
  SmilesSettingsExpected{"N1=CC=CC=CC=C1", EVERYTHING_HAS_A_PI_ELECTRON, 5, 9, 1, "[n]1ccccccc1"}
));

}  //namespace
