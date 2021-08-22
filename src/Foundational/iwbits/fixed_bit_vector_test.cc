// Tests for FixedbitVector

#include "fixed_bit_vector.h"

#include "googlemock/include/gmock/gmock.h"
#include "googletest/include/gtest/gtest.h"
#include "google/protobuf/text_format.h"

namespace {
using fixed_bit_vector::FixedBitVector;

// Really just a test that nothing crashes.
TEST(TestFixedBitVector, Empty) {
  FixedBitVector foo;
  EXPECT_EQ(1, 1);
}

TEST(TestFixedBitVector, ResizeFromNothing1) {
  FixedBitVector foo;
  foo.resize(3);
  EXPECT_EQ(foo.nbits(), 64);
  for (int b = 0; b < 64; ++b) {
    EXPECT_FALSE(foo.is_set(b));
  }
}

TEST(TestFixedBitVector, ResizeFromNothing2) {
  FixedBitVector foo;
  foo.resize(125);
  EXPECT_EQ(foo.nbits(), 64 * 2);
  for (int b = 0; b < 64 * 2; ++b) {
    EXPECT_FALSE(foo.is_set(b));
  }
}

TEST(TestFixedBitVector, ResizeFromNothingExact) {
  FixedBitVector foo;
  foo.resize(128);
  EXPECT_EQ(foo.nbits(), 64 * 2);
  for (int b = 0; b < 64 * 2; ++b) {
    EXPECT_FALSE(foo.is_set(b));
  }
}

TEST(TestFixedBitVector, BitsSpecifiedInConstructor) {
  FixedBitVector foo(231);
  EXPECT_EQ(foo.nbits(), 64 * 4);
  for (int b = 0; b < 64 * 4; ++b) {
    EXPECT_FALSE(foo.is_set(b));
  }
}

TEST(TestFixedBitVector, BitsSpecifiedInConstructorExact) {
  FixedBitVector foo(256);
  EXPECT_EQ(foo.nbits(), 64 * 4);
  for (int b = 0; b < 64 * 4; ++b) {
    EXPECT_FALSE(foo.is_set(b));
  }
}

TEST(TestFixedBitVector, TestSet) {
  constexpr int nbits = 256;
  FixedBitVector foo(nbits);
  for (int i = 0; i < nbits; ++i) {
    foo.set_bit(i);
    for (int j = 0; j < nbits; ++j) {
      if (i == j)
        EXPECT_TRUE(foo.is_set(j));
      else
        EXPECT_FALSE(foo.is_set(j));
    }
    foo.unset_bit(i);
  }
}

TEST(TestFixedBitVector, TestBitsInCommon) {
  constexpr int nbits = 128;
  FixedBitVector f1(nbits), f2(nbits);
  for (int i = 0; i < nbits; ++i) {
    f1.set_bit(i);
    for (int j = 0; j < nbits; ++j) {
      f2.set_bit(j);
      if (i == j) {
        EXPECT_EQ(f1.BitsInCommon(f2), 1);
      } else {
        EXPECT_EQ(f1.BitsInCommon(f2), 0);
      }
      f2.unset_bit(j);
    }
    f1.unset_bit(i);
  }
}

TEST(TestFixedBitVector, TestNsetSingleBot) {
  constexpr int nbits = 512;
  FixedBitVector foo(nbits);
  for (int i = 0; i < nbits; ++i) {
    foo.set_bit(i);
    EXPECT_EQ(foo.nset(), 1);
    foo.unset_bit(i);
  }
}

TEST(TestFixedBitVector, TestNsetAll) {
  constexpr int nbits = 2048;
  FixedBitVector foo(nbits);
  for (int i = 0; i < nbits; ++i) {
    foo.set_bit(i);
    EXPECT_EQ(foo.nset(), i + 1);
  }
}

TEST(TestLog2, TestSingleBitSet) {
  constexpr uint64_t one = 1;
  for (int i = 1; i < 64; ++i) {
    uint64_t v = one << i;
    EXPECT_EQ(fixed_bit_vector::first_bit_set(v), i);
  }
}

TEST(TestLog2, TestOtherValues) {
  constexpr uint64_t one = 1;
  for (int i = 1; i < 63; ++i) {
    uint64_t v = one << i;
    EXPECT_EQ(fixed_bit_vector::first_bit_set(v), i);
    v |= one << (i+1);
    EXPECT_EQ(fixed_bit_vector::first_bit_set(v), i);
    v |= one << (i-1);
    EXPECT_EQ(fixed_bit_vector::first_bit_set(v), i - 1);
  }
}

TEST(TestFixedBitVector, TestFirstBitSetSingleValue) {
  constexpr int nbits = 2048;
  FixedBitVector foo(nbits);
  for (int i = 0; i < nbits; ++i) {
    foo.set_bit(i);
    EXPECT_EQ(foo.FirstBitSet(), i);
    foo.reset();
  }
}

TEST(TestLog2, TestUnsetSingleValues) {
  constexpr uint64_t one = 1;
  EXPECT_EQ(fixed_bit_vector::first_unset_bit(0), 0);
  EXPECT_LT(fixed_bit_vector::first_unset_bit(std::numeric_limits<uint64_t>::max()), 0);

  for (int i = 1; i < 64; ++i) {
    uint64_t v = one << i;
    EXPECT_EQ(fixed_bit_vector::first_unset_bit(v), 0);
  }
}

TEST(TestLog2, TestUnsetSingleUnset) {
  constexpr uint64_t one = 1;
  for (int i = 0; i < 64; ++i) {
    uint64_t v = std::numeric_limits<uint64_t>::max();
    v ^= (one << i);
    EXPECT_EQ(fixed_bit_vector::first_unset_bit(v), i);
  }
}

TEST(TestLog2, TestUnsetManyUnset) {
  constexpr uint64_t one = 1;
  for (int i = 0; i < 64; ++i) {
    uint64_t v = std::numeric_limits<uint64_t>::max();
    v ^= (one << i);
    for (int j = i + 1; j < 64; j += 3) {
      v ^= (one << j);
    }
    EXPECT_EQ(fixed_bit_vector::first_unset_bit(v), i);
  }
}

TEST(TestFixedBitVector, TestFirstBitMultipleValues) {
  constexpr int nbits = 1024;
  for (int i = 0; i < nbits - 64; ++i) {
    FixedBitVector foo(nbits);
    foo.set_bit(i);
    foo.set_bit(i + 32);
    foo.set_bit(i + 33);
    EXPECT_EQ(foo.FirstBitSet(), i);
  }
}

TEST(TestFixedBitVector, TestOperatorLtLtSingleBit) {
  constexpr int nbits = 64;
  constexpr uint64_t one = 1;
  for (int i = 0; i < nbits; ++i) {
    FixedBitVector foo(nbits);
    foo.set_bit(i);
    std::stringstream ss;
    ss << foo;
    std::stringstream expected;
    expected << "FixedBitVector 64 " << std::hex << (one << i);
    EXPECT_EQ(ss.str(), expected.str());
  }
}
}  // namespace
