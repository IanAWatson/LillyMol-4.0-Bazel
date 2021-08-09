// Test the coordinate box idea.

#include "coordinate_box.h"

#include "googlemock/include/gmock/gmock.h"
#include "googletest/include/gtest/gtest.h"
#include "google/protobuf/text_format.h"


namespace {

using coordinate_box::CoordinateBox;

TEST(TestCoordinateBox, BadBoxSpec1) {
  const const_IWSubstring s;
  CoordinateBox box;
  EXPECT_EQ(box.BuildFromSmilesToken(s), 0);
}
TEST(TestCoordinateBox, BadBoxSpecTooManyTokens) {
  const const_IWSubstring s = "1,2,3,4";
  CoordinateBox box;
  EXPECT_EQ(box.BuildFromSmilesToken(s), 0);
}
TEST(TestCoordinateBox, BadBoxSpecBadCellSize1) {
  const const_IWSubstring s = "0.0,12.0,3.0";
  CoordinateBox box;
  EXPECT_EQ(box.BuildFromSmilesToken(s), 0);
}
TEST(TestCoordinateBox, BadBoxSpecBadCellSize2) {
  const const_IWSubstring s = "-0.2,12.0,3.0";
  CoordinateBox box;
  EXPECT_EQ(box.BuildFromSmilesToken(s), 0);
}
TEST(TestCoordinateBox, BadBoxSpecBadX) {
  const const_IWSubstring s = "0.1,0.0,3.0";
  CoordinateBox box;
  EXPECT_EQ(box.BuildFromSmilesToken(s), 0);
}
TEST(TestCoordinateBox, BadBoxSpecBadY) {
  const const_IWSubstring s = "0.1,1.0,0.0";
  CoordinateBox box;
  EXPECT_EQ(box.BuildFromSmilesToken(s), 0);
}
TEST(TestCoordinateBox, BadBoxSpecNotAMultipleX) {
  const const_IWSubstring s = "0.1,1.05,1.0";
  CoordinateBox box;
  EXPECT_EQ(box.BuildFromSmilesToken(s), 0);
}
TEST(TestCoordinateBox, BadBoxSpecNotAMultipleY) {
  const const_IWSubstring s = "0.1,1.00,0.93";
  CoordinateBox box;
  EXPECT_EQ(box.BuildFromSmilesToken(s), 0);
}
TEST(TestCoordinateBox, OriginIsZero) {
  const const_IWSubstring s = "0.1,1.00,1.00";
  CoordinateBox box;
  ASSERT_TRUE(box.BuildFromSmilesToken(s));
  EXPECT_EQ(box.CellNumber(0.0, 0.0, 0.0), 0);
}
#ifdef CELL_NUMBER_CHECK_IN_BOUNDS
TEST(TestCoordinateBox, OutOfRangeXneg) {
  const const_IWSubstring s = "0.1,1.00,1.00";
  CoordinateBox box;
  ASSERT_TRUE(box.BuildFromSmilesToken(s));
  EXPECT_LT(box.CellNumber(-0.1, 0.0, 0.0), 0);
}
TEST(TestCoordinateBox, OutOfRangeX) {
  const const_IWSubstring s = "0.1,1.00,1.00";
  CoordinateBox box;
  ASSERT_TRUE(box.BuildFromSmilesToken(s));
  EXPECT_LT(box.CellNumber(1.1, 0.0, 0.0), 0);
}
TEST(TestCoordinateBox, OutOfRangeYneg) {
  const const_IWSubstring s = "0.1,1.00,1.00";
  CoordinateBox box;
  ASSERT_TRUE(box.BuildFromSmilesToken(s));
  EXPECT_LT(box.CellNumber(0.1, -0.1, 0.0), 0);
}
TEST(TestCoordinateBox, OutOfRangeY) {
  const const_IWSubstring s = "0.1,1.00,1.00";
  CoordinateBox box;
  ASSERT_TRUE(box.BuildFromSmilesToken(s));
  EXPECT_LT(box.CellNumber(0.5, 1.01, 0.0), 0);
}
TEST(TestCoordinateBox, OutOfRangeZneg) {
  const const_IWSubstring s = "0.1,1.00,1.00";
  CoordinateBox box;
  ASSERT_TRUE(box.BuildFromSmilesToken(s));
  EXPECT_LT(box.CellNumber(0.1, 0.7, -1.001), 0);
}
#endif
TEST(TestCoordinateBox, AllCloseX) {
  const const_IWSubstring s = "0.001,1.00,1.00";
  CoordinateBox box;
  ASSERT_TRUE(box.BuildFromSmilesToken(s));
  for (int i = 0; i < 1000; ++i) {
    const float x = i * 0.001;
    const float y = 0.0;
    const float z = 0.0;
    Space_Vector<float> initial(x, y, z);
    const int cell_number = box.CellNumber(x, y, z);
    ASSERT_GE(cell_number, 0);
    const Space_Vector<float> back = box.CoordinatesAsVector<float>(cell_number);
    EXPECT_LT(back.distance(initial), 0.001);
  }
}
TEST(TestCoordinateBox, AllCloseY) {
  const const_IWSubstring s = "0.001,1.00,1.00";
  CoordinateBox box;
  ASSERT_TRUE(box.BuildFromSmilesToken(s));
  for (int i = 0; i < 1000; ++i) {
    const float x = 0.0;
    const float y = i * 0.001;
    const float z = 0.0;
    Space_Vector<float> initial(x, y, z);
    const int cell_number = box.CellNumber(x, y, z);
    ASSERT_GE(cell_number, 0);
    const Space_Vector<float> back = box.CoordinatesAsVector<float>(cell_number);
    EXPECT_LT(back.distance(initial), 0.001);
  }
}
TEST(TestCoordinateBox, AllCloseZ) {
  const const_IWSubstring s = "0.001,1.00,1.00";
  CoordinateBox box;
  ASSERT_TRUE(box.BuildFromSmilesToken(s));
  for (int i = 0; i < 1000; ++i) {
    const float x = 0.0;
    const float y = 0.0;
    const float z = i * 0.001;
    Space_Vector<float> initial(x, y, z);
    const int cell_number = box.CellNumber(x, y, z);
    ASSERT_GE(cell_number, 0);
    const Space_Vector<float> back = box.CoordinatesAsVector<float>(cell_number);
    EXPECT_LT(back.distance(initial), 0.001);
  }
}
TEST(TestCoordinateBox, ArbitraryPositions) {
  const const_IWSubstring s = "0.001,1.00,2.00";
  CoordinateBox box;
  ASSERT_TRUE(box.BuildFromSmilesToken(s));
  for (int i = 0; i < 37; ++i) {
    const float x = static_cast<double>(i) / 37.0;
    for (int j = 0; j < 57; ++j) {
      const float y = static_cast<double>(j) / 57.0;
      for (int k = 0; k < 87; ++k) {
        const float z = static_cast<double>(k) / 87.0;
        Space_Vector<float> initial(x, y, z);
        const int cell_number = box.CellNumber(x, y, z);
        ASSERT_GE(cell_number, 0);
        const Space_Vector<float> back = box.CoordinatesAsVector<float>(cell_number);
        EXPECT_LT(back.distance(initial), 0.001);
      }
    }
  }
}

}  // namespace
