// Tests for smiles reading and generation

#include "googlemock/include/gmock/gmock.h"
#include "googletest/include/gtest/gtest.h"

#include "smiles.h"

namespace {

// For tests that need an input smiles and an expected output.
struct SmilesAndExpected {
  IWString input_smiles;
  IWString expected_smiles;
};

// For tests that need an input, a parameter to set, and an expected output.
struct SmilesAndParam {
  IWString input_smiles;
  int int_param;
  IWString expected_smiles;
};

class TestSmilesCharges : public testing::TestWithParam<SmilesAndParam> {
};

TEST_P(TestSmilesCharges, MultiplePlusCharges) {
  const auto params = GetParam();
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(params.input_smiles));
  set_write_formal_charge_as_consecutive_signs(params.int_param);
  EXPECT_EQ(m.smiles(), params.expected_smiles);
}

INSTANTIATE_TEST_SUITE_P(TestMultipleCharges, TestSmilesCharges, testing::Values(
   SmilesAndParam{"C", 0, "C"},
   SmilesAndParam{"C", 1, "C"},
   SmilesAndParam{"[Fe+++]", 1, "[Fe+++]"},
   SmilesAndParam{"[Fe+++]", 0, "[Fe+3]"},
   SmilesAndParam{"[P---]", 1, "[P---]"},
   SmilesAndParam{"[P---]", 0, "[P-3]"},
   SmilesAndParam{"[Zn--]", 1, "[Zn--]"},
   SmilesAndParam{"[Zn--]", 0, "[Zn-2]"}
));

class TestAromaticSmiles : public testing::TestWithParam<SmilesAndExpected> {
  public:
    void SetUp() {
      set_include_aromaticity_in_smiles(1);
    }
};

TEST_P(TestAromaticSmiles, AromaticSmiles) {
  const auto params = GetParam();
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(params.input_smiles));
  EXPECT_EQ(m.smiles(), params.expected_smiles);
}

INSTANTIATE_TEST_SUITE_P(TestAromaticSmiles, TestAromaticSmiles, testing::Values(
  SmilesAndExpected{"C", "C"},
  SmilesAndExpected{"C1=CC=CC=C1", "c1ccccc1"}
));

class TestHOnAromaticNP: public testing::TestWithParam<SmilesAndParam> {
};

TEST_P(TestHOnAromaticNP, TestHOnAromaticNP) {
  const auto params = GetParam();
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(params.input_smiles));
  set_include_implicit_hydrogens_on_aromatic_n_and_p(params.int_param);
  EXPECT_EQ(m.unique_smiles(), params.expected_smiles);
}

INSTANTIATE_TEST_SUITE_P(TestHOnAromaticNP, TestHOnAromaticNP, testing::Values(
  SmilesAndParam{"C", 1, "C"},
  SmilesAndParam{"C1=CC=CN1", 1, "[nH]1cccc1"},
  SmilesAndParam{"C1=CC=CN1", 0, "n1cccc1"}
));

class TestAddHToIsotopes : public testing::TestWithParam<SmilesAndParam> {
};

TEST_P(TestAddHToIsotopes, TestAddHToIsotopes) {
  const auto params = GetParam();
  set_unset_implicit_hydrogens_known_if_possible(params.int_param);
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(params.input_smiles));
  EXPECT_EQ(m.unique_smiles(), params.expected_smiles);
}

INSTANTIATE_TEST_SUITE_P(TestAddHToIsotopes, TestAddHToIsotopes, testing::Values(
  SmilesAndParam{"C", 1, "C"},
  SmilesAndParam{"[2C]", 1, "[2C]"},
  SmilesAndParam{"[2C]C", 1, "[2C]C"},
  SmilesAndParam{"C[2C]C", 1, "C[2C]C"}
//SmilesAndParam{"[9C]1=CC=CC=C1", 1, "[9cH]1ccccc1"}  not sure why this does not work, investigate.
));


}  // namespace
