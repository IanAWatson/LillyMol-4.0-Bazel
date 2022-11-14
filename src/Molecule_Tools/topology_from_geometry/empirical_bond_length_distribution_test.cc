// Tester for the Empirical Bond Length Distribution class

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "Molecule_Tools/topology_from_geometry/empirical_bond_length_distribution.h"
#include "Molecule_Tools/topology_from_geometry/bond_length_distribution.pb.h"

namespace {

using topology_from_geometry::EmpiricalBondLengthDistribution;
using topology_from_geometry::EmpiricalBondLengthDistributions;

TEST(TestPdf, TestFlat) {
  const std::string string_proto = R"pb(
    atomic_number1: 6
    atomic_number2: 6
    btype: SS_SINGLE_BOND

    min_distance: 1.4
    dx: 0.01
    count: [
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1
    ]
)pb";

  BondLengthDistribution::Distribution proto;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(string_proto, &proto));
  EXPECT_EQ(proto.count_size(), 10);

  EmpiricalBondLengthDistribution bld;
  ASSERT_TRUE(bld.Build(proto));
  float expected_max_distance = proto.min_distance() + proto.count_size() * proto.dx();
  EXPECT_NEAR(bld.max_distance(), expected_max_distance, 0.001);

  EXPECT_FLOAT_EQ(bld.Score(proto.min_distance() - 0.001), 0.0f);
  EXPECT_FLOAT_EQ(bld.Score(expected_max_distance + 0.001), 0.0f);

  EXPECT_NEAR(bld.Score(proto.min_distance()), 1.0/30.0f, 0.00001);
  EXPECT_FLOAT_EQ(bld.Score(expected_max_distance), 0.0f);
}
}  // namespace
