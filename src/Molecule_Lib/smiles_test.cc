// Tests for smiles reading and generation

#include "googlemock/include/gmock/gmock.h"
#include "googletest/include/gtest/gtest.h"

#include "smiles.h"

namespace {

struct TestCharges {
  IWString input_smiles;
  int write_formal_charge_as_consecutive_signs;
  IWString expected_smiles;
};

class TestSmilesCharges : public testing::TestWithParam<TestCharges> {
};

TEST_P(TestSmilesCharges, MultiplePlusCharges) {
  const auto params = GetParam();
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(params.input_smiles));
  set_write_formal_charge_as_consecutive_signs(params.write_formal_charge_as_consecutive_signs);
  EXPECT_EQ(m.smiles(), params.expected_smiles);
}

INSTANTIATE_TEST_SUITE_P(TestMultipleCharges, TestSmilesCharges, testing::Values(
   TestCharges{"C", 0, "C"},
   TestCharges{"C", 1, "C"},
   TestCharges{"[Fe+++]", 1, "[Fe+++]"},
   TestCharges{"[Fe+++]", 0, "[Fe+3]"},
   TestCharges{"[P---]", 1, "[P---]"},
   TestCharges{"[P---]", 0, "[P-3]"},
   TestCharges{"[Zn--]", 1, "[Zn--]"},
   TestCharges{"[Zn--]", 0, "[Zn-2]"}
));
}  // namespace
