// Tests for tanimoto_float

#include <iostream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "tanimoto_float.h"

namespace {

using tanimoto_float::AvxTanimoto;

// Compute the tanimoto coefficient scalar.
float
Expected(const float* v1, const float* v2, int n) {
  float result = 0.0f;
  for (int i = 0; i < n; ++i) {
    if (v1[i] < v2[i]) {
      result += v1[i] / v2[i];
    } else if (v1[i] > v2[i]) {
      result += v2[i] / v1[i];
    } else {
      result += 1.0f;
    }
  }

  return result / static_cast<float>(n);
}

TEST(TestTanimoto, All1) {
  alignas(32) float v1[]{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  alignas(32) float v2[]{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  EXPECT_FLOAT_EQ(AvxTanimoto(v1, v2, 8), 1.0);
}

TEST(TestTanimoto, OneandTwo) {
  alignas(32) float v1[]{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  alignas(32) float v2[]{2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
  EXPECT_FLOAT_EQ(AvxTanimoto(v1, v2, 8), 0.5);
}

TEST(TestTanimoto, OneandFour) {
  alignas(32) float v1[]{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  alignas(32) float v2[]{4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0};
  EXPECT_FLOAT_EQ(AvxTanimoto(v1, v2, 8), 0.25);
  EXPECT_FLOAT_EQ(AvxTanimoto(v2, v1, 8), 0.25);
}

TEST(TestTanimoto, OneandVarious1) {
  alignas(32) float v1[]{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  alignas(32) float v2[]{2.0, 4.0, 2.0, 4.0, 2.0, 4.0, 2.0, 4.0};
  EXPECT_FLOAT_EQ(AvxTanimoto(v1, v2, 8), 0.375);
  EXPECT_FLOAT_EQ(AvxTanimoto(v2, v1, 8), 0.375);
}

TEST(TestTanimoto, OneandVarious2) {
  alignas(32) float v1[]{5.0, 4.0, 2.0, 4.0, 2.0, 4.0, 2.0, 4.0};
  alignas(32) float v2[]{1.0, 3.0, 2.0, 3.0, 3.0, 4.0, 3.0, 8.0};
  const float expected = Expected(v1, v2, 8);
  EXPECT_FLOAT_EQ(AvxTanimoto(v1, v2, 8), expected);
  EXPECT_FLOAT_EQ(AvxTanimoto(v2, v1, 8), expected);
}

}  // namespace
