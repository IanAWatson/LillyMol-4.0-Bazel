// Tester for scaffold enumeration

#include "googletest/include/gtest/gtest.h"
#include "googlemock/include/gmock/gmock.h"

#include "scaffolds.h"

namespace {

using testing::UnorderedElementsAreArray;

struct MoleculeResult {
  IWString smiles;
  resizable_array<IWString> results;
};

class TestScaffolds : public testing::TestWithParam<MoleculeResult> {
  protected:
    Molecule _m;

    scaffolds::ScaffoldFinder _scaffold_finder;

    Scaffolds::ScaffoldData _scaffolds;

    // To do the final test, we convert the scaffold molecules to string form.
    resizable_array<IWString> _scaffold_smiles;
};

TEST_P(TestScaffolds, TestScaffolds) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));

#ifdef DEBUG_TEST_ADASD
  std::cerr << "Processing " << params.smiles << " results " << params.results.size() << " items\n";
  for (const IWString& s : params.results) {
    std::cerr << "Expecting '" << s << "'\n";
  }
#endif
  _scaffold_finder.MakeScaffolds(_m, _scaffolds);

#ifdef DEBUG_TEST_ADASD
  std::cerr << "Got " << _scaffolds.size() << " values back\n";
  std::cerr << "params.results.size() " << params.results.size() << '\n';
#endif

  EXPECT_EQ(_scaffolds.subset_size(), params.results.number_elements());

  if (_scaffolds.subset().empty()) {
    return;
  }

  for (const auto& p : _scaffolds.subset()) {
    Molecule m;
    const IWString smiles = p.smi();
    m.build_from_smiles(smiles);
    // std::cerr << "Generated " << m->smiles() << '\n';
    _scaffold_smiles << m.smiles();
  }

#ifdef DEBUG_TEST_ADASD
  for (const auto& s : params.results) {
    std::cerr << "Expecting " << s << '\n';
  }
#endif

  EXPECT_THAT(_scaffold_smiles, UnorderedElementsAreArray(params.results));
}

INSTANTIATE_TEST_SUITE_P(TestScaffolds, TestScaffolds, testing::Values(
  MoleculeResult{"C", {}},
  MoleculeResult{"C1CC1", {{"C1CC1"}}},
  MoleculeResult{"C1CC1C", {{"C1CC1"}}},
  MoleculeResult{"C12CC1C2CC1CC1", {{"C1CC1", "C12CC1C2CC1CC1", "C12CC1C2"}}},
  MoleculeResult{"C1CC1CC1CC1", {{"C1CC1", "C1CC1CC1CC1"}}},
  MoleculeResult{"C1CC1CC1CC1CC1CC1", {{"C1CC1", "C1CC1CC1CC1", "C1CC1CC1CC1CC1CC1"}}},
  MoleculeResult{"C1CC1C(C1CC1)C1CC1", {{"C1CC1", "C(C1CC1)C1CC1", "C1CC1C(C1CC1)C1CC1"}}}
));

}  // namespace
