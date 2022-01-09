// Tester for fileconv

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "fileconv_opts.h"

namespace {

using std::vector;

class FileconvTest : public testing::Test {
  protected:
    Molecule _m;
    fileconv::FileconvConfig _config;
    void SetUp();
    void TearDown();
};

void
FileconvTest::SetUp() {
}

void
FileconvTest::TearDown() {
}

TEST_F(FileconvTest, TestLongestPathPasses) {
  int argc = 3;
  const char * argv[] = {"fileconv", "-Y", "maxplen=4"};
  Command_Line cl(argc, const_cast<char**>(argv), "Y:");
  ASSERT_TRUE(_config.Build(cl));

  ASSERT_TRUE(_m.build_from_smiles("CCCCC"));
  fileconv::FileconvResult result = _config.Process(_m);
  EXPECT_FALSE(result.rejected);
}

TEST_F(FileconvTest, TestLongestPathRejected) {
  int argc = 3;
  const char * argv[] = {"fileconv", "-Y", "maxplen=4"};
  Command_Line cl(argc, const_cast<char**>(argv), "Y:");
  ASSERT_TRUE(_config.Build(cl));

  ASSERT_TRUE(_m.build_from_smiles("CCCCCC"));
  fileconv::FileconvResult result = _config.Process(_m);
  EXPECT_TRUE(result.rejected);
}

// Pattern is established, now write parameterized tests.

struct OptionsSmilesResult {
  // First inputs.
  std::vector<const char*> argv;
  const char * options;
  const char * smiles;
  const char * name;

  // And then the results.
  int expected_rejected;
  // For actions that change the name of the molecule.
  const char * expected_smiles;
  const char * expected_name;
};

class FileConvTestP : public testing::TestWithParam<OptionsSmilesResult> {
  protected:
    Molecule _m;
    fileconv::FileconvConfig _config;
};

TEST_P(FileConvTestP, FileConvTestP) {
  const auto& params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  _m.set_name(params.name);
  const int argc = params.argv.size();
  ASSERT_GT(argc, 0);
  char ** argv = new char*[argc];
  std::unique_ptr<char *[]> free_argv(argv);
  int ndx = 0;
  for (const auto& c : params.argv) {
    argv[ndx] = const_cast<char*>(c);
    ++ndx;
  }
  Command_Line cl(argc, argv, params.options);
  ASSERT_TRUE(_config.Build(cl));
  fileconv::FileconvResult result = _config.Process(_m);
  EXPECT_EQ(result.rejected, params.expected_rejected);
  EXPECT_EQ(_m.smiles(), params.expected_smiles);
  EXPECT_EQ(_m.name(), params.expected_name);
}
INSTANTIATE_TEST_SUITE_P(FileConvTestP, FileConvTestP, testing::Values(
  OptionsSmilesResult{vector<const char*>{"_", "-c", "4"}, "c:", 
                      "CC", "name", 1, "CC", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-c", "4"}, "c:",
                      "CCCC", "name", 0, "CCCC", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-C", "4"}, "C:",
                      "CCCCCCC", "name", 1, "CCCCCCC", "name"},
  OptionsSmilesResult{vector<const char*>{"_", "-C", "4"}, "C:",
                      "CCCC", "name", 0, "CCCC", "name"},

  OptionsSmilesResult{vector<const char*>{"_", "-p", "AMW"}, "p:", 
                      "OC(F)(Cl)C(P)(S)NI", "name", 0, "OC(F)(Cl)C(P)(S)NI", "name AMW = 303.4616"}
));

}  // namespace
