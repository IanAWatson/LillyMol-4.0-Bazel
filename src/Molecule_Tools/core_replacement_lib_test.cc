// Tests for the core replacement lib

#include <memory>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "Foundational/iwmisc/misc.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/substructure.h"

#include "Molecule_Tools/core_replacement_lib.h"

namespace {
TEST(TestIdentifyAtomsToRemove, rm1) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CO[1C](=O)C1=C(C=CS1)[1NH]C(=N)N"));
  const int matoms = m.natoms();
  Substructure_Query query;
  query.create_from_smarts("[1].[1]");
  query.set_find_unique_embeddings_only(1);
  Substructure_Results sresults;

  int nhits = query.substructure_search(m, sresults);
  ASSERT_EQ(nhits, 1);
  std::unique_ptr<int[]> interior(new_int(matoms));
  sresults.each_embedding_set_vector(interior.get(), 1);
  EXPECT_EQ(separated_atoms::IdentifyAtomsToRemove(m, *sresults.embedding(0), interior.get()), 1);
  std::vector<int> got(interior.get(), interior.get() + matoms);
  EXPECT_THAT(got, testing::ElementsAreArray({0, 0, 1, 0, 2, 2, 2, 2, 2, 1, 0, 0, 0}));
}

class TestIdentifyCoreAtoms : public testing::Test {
  protected:
    Molecule _m;
    Substructure_Query _query;
    Substructure_Results _sresults;
    std::unique_ptr<int[]> _interior;
    std::vector<int> _got;
};

TEST_F(TestIdentifyCoreAtoms, Test1) {
  ASSERT_TRUE(_m.build_from_smiles("C[1CH2]c1cc([1CH2]O)cc([1CH2]N)c1"));
  ASSERT_TRUE(_query.create_from_smarts("[1].[1].[1]"));
  _query.set_find_unique_embeddings_only(1);
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const int matoms = _m.natoms();
  _interior.reset(new_int(matoms));
  _sresults.each_embedding_set_vector(_interior.get(), 1);
  EXPECT_EQ(separated_atoms::IdentifyAtomsToRemove(_m, *_sresults.embedding(0), _interior.get()), 1);
  _got.assign(_interior.get(), _interior.get() + matoms);
  EXPECT_THAT(_got, testing::ElementsAreArray({0, 1, 2, 2, 2, 1, 0, 2, 2, 1, 0, 2}));
}

TEST_F(TestIdentifyCoreAtoms, Bad1) {
  ASSERT_TRUE(_m.build_from_smiles("CC[1c]1c[1c](CO)c[1c](CN)c1"));
  ASSERT_TRUE(_query.create_from_smarts("[1].[1].[1]"));
  _query.set_find_unique_embeddings_only(1);
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const int matoms = _m.natoms();
  _interior.reset(new_int(matoms));
  _sresults.each_embedding_set_vector(_interior.get(), 1);
  EXPECT_EQ(separated_atoms::IdentifyAtomsToRemove(_m, *_sresults.embedding(0), _interior.get()), 0);
}

}  // namespace
