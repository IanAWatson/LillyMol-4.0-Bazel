// Tests for 3D reactions

#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "Molecule_Lib/iwreaction.h"
#include "Molecule_Lib/smiles.h"

namespace {

TEST(TestBumpCheck, NoProblem) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("N{{0,1,0}}C{{0,0,0}}C{{1,0,0}}C{{1,1,0}}"));
  const std::string string_proto = R"pb(
scaffold {
  id: 0
  smarts: "NCCC"
  dihedral_angle {
    a1: 0
    a2: 1
    a3: 2
    a4: 3
    angle: 90.0
  }
}
  )pb";

  ReactionProto::Reaction proto;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(string_proto, &proto));
  IWReaction rxn;
  ASSERT_TRUE(rxn.ConstructFromProto(proto));
  EXPECT_EQ(rxn.in_place_transformations(m), 1);
}

struct SmilesAngleDistance {
  IWString smiles;
  float angle;
  int expected;
};

class TestDihedral : public testing::TestWithParam<SmilesAngleDistance> {
  protected:
    Molecule _m;
    ReactionProto::Reaction _proto;
    IWReaction _rxn;
};

TEST_P(TestDihedral, TestVarious) {
  const auto& params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  const std::string string_proto = R"pb(
scaffold {
  id: 0
  smarts: "NCCC"
  dihedral_angle {
    a1: 0
    a2: 1
    a3: 2
    a4: 3
    angle: 90.0
  }
}
  )pb";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(string_proto, &_proto));
  _proto.mutable_scaffold()->mutable_dihedral_angle(0)->set_angle(params.angle);
  EXPECT_TRUE(_rxn.ConstructFromProto(_proto));
  EXPECT_EQ(_rxn.in_place_transformations(_m), params.expected);
}
INSTANTIATE_TEST_SUITE_P(TestDihedral, TestDihedral, testing::Values(
  SmilesAngleDistance{"N{{0,1,0}}C{{0,0,0}}C{{1,0,0}}C{{1,1,0}}", 0.0, 1},
  SmilesAngleDistance{"N{{0,1,0}}C{{0,0,0}}C{{1,0,0}}C{{1,1,0}}", 5.0, 1},
  SmilesAngleDistance{"N{{0,1,0}}C{{0,0,0}}C{{1,0,0}}C{{1,1,0}}", -5.0, 1},
  SmilesAngleDistance{"N{{0,1,0}}C{{0,0,0}}C{{1,0,0}}C{{1,1,0}}", 5.0, 1}
));

}  // namespace
