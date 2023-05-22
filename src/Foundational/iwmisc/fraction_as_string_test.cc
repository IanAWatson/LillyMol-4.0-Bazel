// Tests for fraction as string

#include "googletest/include/gtest/gtest.h"

#include "iwdigits.h"

namespace {
struct FloatAndString {
  // The parameters of the Fraction_As_String
  float minval;
  float maxval;
  int digits;

  // An input value and what should return.
  float value;
  IWString expected;
};

class TestFractionAsString: public testing::TestWithParam<FloatAndString> {
  protected:
    Fraction_as_String _fraction_as_string;

    IWString _as_string;
};

TEST_P(TestFractionAsString, Lots) {
  const auto params = GetParam();
  ASSERT_TRUE(_fraction_as_string.initialise(params.minval, params.maxval, params.digits));
  _fraction_as_string.append_number(_as_string, params.value);
  EXPECT_EQ(_as_string, params.expected) <<
     "mismatch from " << params.value << " encoded " << _as_string << " expected " << params.expected;
}
INSTANTIATE_TEST_SUITE_P(TestFractionAsString, TestFractionAsString, testing::Values(
   FloatAndString{0.0, 1.0, 1, 0.5, "0.5"},
   FloatAndString{0.0, 1.0, 2, 0.5, "0.50"},
   // The ends are treated specially.
   FloatAndString{0.0, 1.0, 2, 0.0, "0"},
   FloatAndString{0.0, 1.0, 2, 1.0, "1"},
   // Out of range
   FloatAndString{0.0, 1.0, 2, 1.1, "1.1"},
   // Works for ranges other than fractions
   FloatAndString{0.0, 2.0, 2, 1.103, "1.10"}
));

}  // namespace
