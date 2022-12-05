// Tests for IW_Logical_Expression

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Foundational/iwmisc/logical_expression.h"
#include "Foundational/iwstring/iwstring.h"

namespace {

TEST(TestLogicalExpression, EmptyCompare) {
  IW_Logical_Expression e1, e2;
  EXPECT_EQ(e1, e2);
}

TEST(TestLogicalExpression, EmptyCompareDuringBUild) {
  IW_Logical_Expression e1, e2;
  EXPECT_EQ(e1, e2);
  e1.add_operator(IW_LOGEXP_OR);
  EXPECT_NE(e1, e2);
  e2.add_operator(IW_LOGEXP_OR);
  EXPECT_EQ(e1, e2);

  // The default unary operator is 1, so they are still equal
  e1.set_unary_operator(0, 1);
  EXPECT_EQ(e1, e2);
  e2.set_unary_operator(0, 1);
  EXPECT_EQ(e1, e2);

  e1.set_unary_operator(0, 0);
  EXPECT_NE(e1, e2);
  e2.set_unary_operator(0, 0);
  EXPECT_EQ(e1, e2);

  EXPECT_TRUE(e1.ok());
  EXPECT_TRUE(e2.ok());
}

TEST(TestLogicalExpression, TestEmptyString) {
  IW_Logical_Expression e1;
  EXPECT_EQ(e1.ToString(), "");
}

TEST(TestLogicalExpression, TestOr) {
  IW_Logical_Expression e1;
  e1.add_operator(IW_LOGEXP_OR);
  EXPECT_EQ(e1.ToString(), ".,.");
  e1.set_unary_operator(0, 0);
  EXPECT_EQ(e1.ToString(), "!.,.");
  e1.set_unary_operator(1, 0);
  EXPECT_EQ(e1.ToString(), "!.,!.");
  e1.set_unary_operator(0, 1);
  EXPECT_EQ(e1.ToString(), ".,!.");
}

TEST(TestLogicalExpression, TestAnd) {
  IW_Logical_Expression e1;
  e1.add_operator(IW_LOGEXP_AND);
  EXPECT_EQ(e1.ToString(), ".&.");
  e1.set_unary_operator(0, 0);
  EXPECT_EQ(e1.ToString(), "!.&.");
  e1.set_unary_operator(1, 0);
  EXPECT_EQ(e1.ToString(), "!.&!.");
  e1.set_unary_operator(0, 1);
  EXPECT_EQ(e1.ToString(), ".&!.");
}

TEST(TestLogicalExpression, TestXor) {
  IW_Logical_Expression e1;
  e1.add_operator(IW_LOGEXP_XOR);
  EXPECT_EQ(e1.ToString(), ".^.");
  e1.set_unary_operator(0, 0);
  EXPECT_EQ(e1.ToString(), "!.^.");
  e1.set_unary_operator(1, 0);
  EXPECT_EQ(e1.ToString(), "!.^!.");
  e1.set_unary_operator(0, 1);
  EXPECT_EQ(e1.ToString(), ".^!.");
}

TEST(TestLogicalExpression, TestLPAnd) {
  IW_Logical_Expression e1;
  e1.add_operator(IW_LOGEXP_LOW_PRIORITY_AND);
  EXPECT_EQ(e1.ToString(), ".;.");
  e1.set_unary_operator(0, 0);
  EXPECT_EQ(e1.ToString(), "!.;.");
  e1.set_unary_operator(1, 0);
  EXPECT_EQ(e1.ToString(), "!.;!.");
  e1.set_unary_operator(0, 1);
  EXPECT_EQ(e1.ToString(), ".;!.");
}

TEST(TestLogicalExpression, TestBuildAnd) {
  IW_Logical_Expression e1;
  e1.add_operator(IW_LOGEXP_AND);
  EXPECT_EQ(e1.ToString(), ".&.");
  IW_Logical_Expression e2;
  ASSERT_TRUE(e2.BuildFromString(e1.ToString()));
  EXPECT_EQ(e1, e2);
  e1.set_unary_operator(0, 0);
  EXPECT_EQ(e1.ToString(), "!.&.");
  ASSERT_TRUE(e2.BuildFromString(e1.ToString()));
  EXPECT_EQ(e1, e2);
  e1.set_unary_operator(1, 0);
  EXPECT_EQ(e1.ToString(), "!.&!.");
  ASSERT_TRUE(e2.BuildFromString(e1.ToString()));
  EXPECT_EQ(e1, e2);
  e1.set_unary_operator(0, 1);
  EXPECT_EQ(e1.ToString(), ".&!.");
  ASSERT_TRUE(e2.BuildFromString(e1.ToString()));
  EXPECT_EQ(e1, e2);
}

TEST(TestLogicalExpression, TestBuildOr) {
  IW_Logical_Expression e1;
  e1.add_operator(IW_LOGEXP_OR);
  EXPECT_EQ(e1.ToString(), ".,.");
  IW_Logical_Expression e2;
  ASSERT_TRUE(e2.BuildFromString(e1.ToString()));
  EXPECT_EQ(e1, e2);
  e1.set_unary_operator(0, 0);
  EXPECT_EQ(e1.ToString(), "!.,.");
  ASSERT_TRUE(e2.BuildFromString(e1.ToString()));
  EXPECT_EQ(e1, e2);
  e1.set_unary_operator(1, 0);
  EXPECT_EQ(e1.ToString(), "!.,!.");
  ASSERT_TRUE(e2.BuildFromString(e1.ToString()));
  EXPECT_EQ(e1, e2);
  e1.set_unary_operator(0, 1);
  EXPECT_EQ(e1.ToString(), ".,!.");
  ASSERT_TRUE(e2.BuildFromString(e1.ToString()));
  EXPECT_EQ(e1, e2);
}

TEST(TestLogicalExpression, TestBuildXOr) {
  IW_Logical_Expression e1;
  e1.add_operator(IW_LOGEXP_XOR);
  EXPECT_EQ(e1.ToString(), ".^.");
  IW_Logical_Expression e2;
  ASSERT_TRUE(e2.BuildFromString(e1.ToString()));
  EXPECT_EQ(e1, e2);
  e1.set_unary_operator(0, 0);
  EXPECT_EQ(e1.ToString(), "!.^.");
  ASSERT_TRUE(e2.BuildFromString(e1.ToString()));
  EXPECT_EQ(e1, e2);
  e1.set_unary_operator(1, 0);
  EXPECT_EQ(e1.ToString(), "!.^!.");
  ASSERT_TRUE(e2.BuildFromString(e1.ToString()));
  EXPECT_EQ(e1, e2);
  e1.set_unary_operator(0, 1);
  EXPECT_EQ(e1.ToString(), ".^!.");
  ASSERT_TRUE(e2.BuildFromString(e1.ToString()));
  EXPECT_EQ(e1, e2);
}

TEST(TestLogicalExpression, TestBuildLPAnd) {
  IW_Logical_Expression e1;
  e1.add_operator(IW_LOGEXP_LOW_PRIORITY_AND);
  EXPECT_EQ(e1.ToString(), ".;.");
  IW_Logical_Expression e2;
  ASSERT_TRUE(e2.BuildFromString(e1.ToString()));
  EXPECT_EQ(e1, e2);
  e1.set_unary_operator(0, 0);
  EXPECT_EQ(e1.ToString(), "!.;.");
  ASSERT_TRUE(e2.BuildFromString(e1.ToString()));
  EXPECT_EQ(e1, e2);
  e1.set_unary_operator(1, 0);
  EXPECT_EQ(e1.ToString(), "!.;!.");
  ASSERT_TRUE(e2.BuildFromString(e1.ToString()));
  EXPECT_EQ(e1, e2);
  e1.set_unary_operator(0, 1);
  EXPECT_EQ(e1.ToString(), ".;!.");
  ASSERT_TRUE(e2.BuildFromString(e1.ToString()));
  EXPECT_EQ(e1, e2);
}

class TestFailParse : public testing::TestWithParam<IWString> {
  protected:
    IW_Logical_Expression _logexp;
};

TEST_P(TestFailParse, TestFailParse) {
  const auto& params = GetParam();
  EXPECT_FALSE(_logexp.BuildFromString(params));
}
INSTANTIATE_TEST_SUITE_P(TestFailParse, TestFailParse, testing::Values(
  "&",
  ".&",
  "!.&",
  "!",
  "..",
  ".!",
  ".!."
));

// A string to construct an IW_Logical_Expression, a set of results
// and an expected result.
struct LogexpValuesResult {
  IWString _as_string;
  resizable_array<int> _values;
  int _result;
};

class TestLogexpValues : public testing::TestWithParam<LogexpValuesResult> {
  protected:
    IW_Logical_Expression _logexp;
};

TEST_P(TestLogexpValues, TestLogexpValues) {
  const auto& params = GetParam();
  ASSERT_TRUE(_logexp.BuildFromString(params._as_string));
  for (int i = 0; i < params._values.number_elements(); ++i) {
    _logexp.set_result(i, params._values[i]);
  }
  int res;
  EXPECT_TRUE(_logexp.evaluate(res));
  EXPECT_EQ(res, params._result);
}

INSTANTIATE_TEST_SUITE_P(TestLogexpValues, TestLogexpValues, testing::Values(
  LogexpValuesResult{ {".,."}, {1,1}, 1},
  LogexpValuesResult{ {".,."}, {0,0}, 0},
  LogexpValuesResult{ {".,."}, {1,0}, 1},
  LogexpValuesResult{ {".,."}, {0,1}, 1},

  LogexpValuesResult{ {".&."}, {0,1}, 0},
  LogexpValuesResult{ {".&."}, {1,0}, 0},
  LogexpValuesResult{ {".&."}, {0,0}, 0},
  LogexpValuesResult{ {".&."}, {1,1}, 1},

  LogexpValuesResult{ {".^."}, {0,1}, 1},
  LogexpValuesResult{ {".^."}, {1,0}, 1},
  LogexpValuesResult{ {".^."}, {0,0}, 0},
  LogexpValuesResult{ {".^."}, {1,1}, 0},

  LogexpValuesResult{ {".;."}, {0,1}, 0},
  LogexpValuesResult{ {".;."}, {1,0}, 0},
  LogexpValuesResult{ {".;."}, {0,0}, 0},
  LogexpValuesResult{ {".;."}, {1,1}, 1},

  LogexpValuesResult{ ".,.,.", {0, 0, 0}, 0},
  LogexpValuesResult{ ".,.,.", {0, 0, 1}, 1},
  LogexpValuesResult{ ".,.,.", {0, 1, 0}, 1},
  LogexpValuesResult{ ".,.,.", {1, 0, 0}, 1},
  LogexpValuesResult{ ".,.,.", {1, 0, 1}, 1},
  LogexpValuesResult{ ".,.,.", {1, 1, 1}, 1},

  LogexpValuesResult{ ".;.,.", {1, 1, 1}, 1},
  LogexpValuesResult{ ".;.,.", {1, 0, 1}, 1},
  LogexpValuesResult{ ".;.,.", {1, 1, 0}, 1},
  LogexpValuesResult{ ".;.,.", {1, 0, 0}, 0},
  LogexpValuesResult{ ".;.,.", {0, 1, 1}, 0},
  LogexpValuesResult{ ".;.,.", {0, 0, 1}, 0},
  LogexpValuesResult{ ".;.,.", {0, 1, 0}, 0},
  LogexpValuesResult{ ".;.,.", {0, 0, 0}, 0},

  LogexpValuesResult{ ".,.;.,.", {0, 0, 0, 0}, 0},
  LogexpValuesResult{ ".,.;.,.", {0, 1, 0, 0}, 0},
  LogexpValuesResult{ ".,.;.,.", {0, 1, 1, 0}, 1},

  LogexpValuesResult{ ".&.;.,.", {0, 1, 1, 0}, 0},
  LogexpValuesResult{ ".&.;.,.", {1, 1, 0, 0}, 0},
  LogexpValuesResult{ ".&.;.,.", {1, 1, 1, 0}, 1},

  LogexpValuesResult{ "!.,.", {1,0}, 0},
  LogexpValuesResult{ "!.,.", {0,0}, 1},
  LogexpValuesResult{ "!.;.", {0,1}, 1},

  LogexpValuesResult{ ".;!.", {0,1}, 0},
  LogexpValuesResult{ ".;!.", {0,0}, 0},
  LogexpValuesResult{ ".;!.", {1,0}, 1},
  LogexpValuesResult{ ".;!.", {1,1}, 0}
));

}  // namespace
