// Tester for the string classes

//#include "googlemock/include/gmock/gmock.h"
//#include "googletest/include/gtest/gtest.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "iwstring.h"

namespace {
TEST(TestIWString, TestAsString) {
  IWString s("hello");
  const std::string as_string = s.AsString();
  EXPECT_EQ(as_string, "hello");
}

TEST(TestConstIWSubstring, TestAsString) {
  const_IWSubstring s("hello");
  const std::string as_string = s.AsString();
  EXPECT_EQ(as_string, "hello");
}

TEST(TestInitializerList, TestInitializerList) {
  resizable_array<IWString> foo({"abc", "def"});
  EXPECT_EQ(foo.size(), 2);
  EXPECT_THAT(foo, testing::ElementsAre("abc", "def"));
}

}  // namespace
