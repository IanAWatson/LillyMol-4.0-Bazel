#include <string>

#include "aromatic.h"
#include "substructure.h"

#include "googlemock/include/gmock/gmock.h"
#include "googletest/include/gtest/gtest.h"
#include "google/protobuf/text_format.h"

namespace {

using testing::UnorderedElementsAre;
using testing::Contains;

class TestSubstructureRingSystem : public testing::Test
{
  protected:
    void SetUp() override;

  protected:
    std::string _string_proto;

    IWString _smiles;

    SubstructureSearch::SubstructureQuery _proto;

    Substructure_Query _query;

    Substructure_Results _sresults;

    Molecule _m;

    static constexpr int _ntest = 10;

  protected:
    void _WriteQuery(const char * fname) {
      IWString tmp(fname);
      _query.write_msi(tmp);
    }

    const Set_of_Atoms _FirstAtomEachEmbedding()
    {
      Set_of_Atoms to_be_returned;
      for (const auto * e : _sresults.embeddings()) {
        to_be_returned.add(e->item(0));
      }

      return to_be_returned;
    }

    bool _DoPermutationTests(const int expected);
};

void
TestSubstructureRingSystem::SetUp()
{
  set_global_aromaticity_type(Daylight);
}

bool TestSubstructureRingSystem::_DoPermutationTests(const int expected)
{
  const IWString initial_smiles = _m.smiles();

  for (int i = 0; i < _ntest; ++i) {
    const IWString smiles = _m.random_smiles();
    if (! _m.build_from_smiles(smiles)) 
    {
      cerr << "_DoPermutationTests:invalid smiles '" << smiles << "'\n";
      return false;
    }

    const int hits = _query.substructure_search(_m, _sresults);

    if (expected == hits)
      continue;

    cerr << "_DoPermutationTests:permutation test failure, expected " << expected << " got " << hits << " hits\n";
    cerr << "Initial smiles " << initial_smiles << endl;
    cerr << "Random  smiles " << smiles << endl;

    return false;
  }

  return true;
}

TEST_F(TestSubstructureRingSystem, TestNoRingsNoMatch)
{
  _string_proto = R"(
    query {
      smarts: "C"
      ring_system_specifier {
        rings_in_system: 1
        base {
          hits_needed: 0
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CCCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPermutationTests(0));
}

TEST_F(TestSubstructureRingSystem, TestRingIsRingSystem)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        rings_in_system: 1
        base {
          hits_needed: 1
        }
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 3);

  EXPECT_TRUE(_DoPermutationTests(3));
}

TEST_F(TestSubstructureRingSystem, TestRingMatches1)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        rings_in_system: 2
        base {
          hits_needed: 1
        }
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1C2CC12";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 4);

  EXPECT_TRUE(_DoPermutationTests(4));
}

TEST_F(TestSubstructureRingSystem, TestRingMatches2)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        rings_in_system: 2
        base {
          hits_needed: 2
        }
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1C2CC12.C1C2CC12.C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 9);

  EXPECT_TRUE(_DoPermutationTests(9));
}

TEST_F(TestSubstructureRingSystem, TestRingHitsNeededNoMatch)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        rings_in_system: 2
        base {
          hits_needed: 1
        }
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1C2CC12.C1C2CC12.C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPermutationTests(0));
}

TEST_F(TestSubstructureRingSystem, TestCubane)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        rings_in_system: 5  # SSSR rings
        base {
          hits_needed: 1
        }
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C12C3C4C1C5C2C3C45";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 8);

  EXPECT_TRUE(_DoPermutationTests(8));
}

TEST_F(TestSubstructureRingSystem, TestRingSizes)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        rings_in_system: 2
        ring_sizes: [3, 4]
        base {
          hits_needed: 1
        }
      }
      smarts: "[0C]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1C2CCC12";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 5);

  EXPECT_TRUE(_DoPermutationTests(5));
}

TEST_F(TestSubstructureRingSystem, TestRingSizesExtra)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        ring_sizes: [3, 4]
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1C2C3CC3C12";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 6);

  EXPECT_TRUE(_DoPermutationTests(6));
}

TEST_F(TestSubstructureRingSystem, TestRingSizeCountMatches)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        ring_size_count {
          ring_size: 3
          count: 2
        }
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1C2C3CC3C12";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 6);

  EXPECT_TRUE(_DoPermutationTests(6));
}

TEST_F(TestSubstructureRingSystem, TestRingSizeCountNoMatches)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        ring_size_count {
          ring_size: 3
          count: 2
        }
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1C2C3CC34CC421";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPermutationTests(0));
}

TEST_F(TestSubstructureRingSystem, TestRingSizeCountMultipleMatches)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        ring_size_count {
          ring_size: 3
          max_count: 3
        }
        ring_size_count {
          ring_size: 4
          min_count: 1
        }
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1C2C3CC34CC421";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 7);

  EXPECT_TRUE(_DoPermutationTests(7));
}

TEST_F(TestSubstructureRingSystem, TestAromaticRingCountNone)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        aromatic_ring_count: 1
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1C2C3CC34CC421";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPermutationTests(0));
}

TEST_F(TestSubstructureRingSystem, TestAromaticRingCountOne)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        aromatic_ring_count: 1
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "o1cc(C)cc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPermutationTests(1));
}

TEST_F(TestSubstructureRingSystem, TestAromaticRingCountTwo)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        min_aromatic_ring_count: 2
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "o1ccc2ccc(C)cc12";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_TRUE(_DoPermutationTests(1));
}

TEST_F(TestSubstructureRingSystem, TestNonAromaticRingNone)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        min_non_aromatic_ring_count: 1
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPermutationTests(0));
}

TEST_F(TestSubstructureRingSystem, TestNonAromaticRingMatchesOne)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        min_non_aromatic_ring_count: 1
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), _m.natoms());

  EXPECT_TRUE(_DoPermutationTests(_m.natoms()));
}

TEST_F(TestSubstructureRingSystem, TestNonAromaticRingMatchesFour)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        max_non_aromatic_ring_count: 5
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1C2C3CC34CC421";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), _m.natoms());

  EXPECT_TRUE(_DoPermutationTests(_m.natoms()));
}

TEST_F(TestSubstructureRingSystem, TestNonAromaticRingMatchesNoSpiro)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        non_aromatic_ring_count: 5
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1C2C3CC34C7(CC7)C421";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 0);

  EXPECT_TRUE(_DoPermutationTests(0));
}

TEST_F(TestSubstructureRingSystem, TestDegreeOfFusionNoneFused)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        degree_of_fusion: 0
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_nhits = 3;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_nhits);

  EXPECT_TRUE(_DoPermutationTests(expected_nhits));
}

TEST_F(TestSubstructureRingSystem, TestDegreeOfFusionOneFused)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        degree_of_fusion: 0
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1C2CC12";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_nhits = 0;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_nhits);

  EXPECT_TRUE(_DoPermutationTests(expected_nhits));
}

TEST_F(TestSubstructureRingSystem, TestDegreeOfFusionOneFusedMatches)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        degree_of_fusion: 1
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1C2CC12";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_nhits = 4;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_nhits);

  EXPECT_TRUE(_DoPermutationTests(expected_nhits));
}

TEST_F(TestSubstructureRingSystem, TestDegreeOfFusionThreeFusedMatches)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        degree_of_fusion: 3
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1C2C3CC34C7(CC7)C421";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_nhits = 9;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_nhits);

  EXPECT_TRUE(_DoPermutationTests(expected_nhits));
}

TEST_F(TestSubstructureRingSystem, TestAtomsInSystemNoSpiro)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        atoms_in_system: 9
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1C2C3CC34C7(CC7)C421";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_nhits = 0;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_nhits);

  EXPECT_TRUE(_DoPermutationTests(expected_nhits));
}

TEST_F(TestSubstructureRingSystem, TestAtomsInSystemXSpiro)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        atoms_in_system: 7
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1C2C3CC34C7(CC7)C421";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_nhits = 9;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_nhits);

  EXPECT_TRUE(_DoPermutationTests(expected_nhits));
}

TEST_F(TestSubstructureRingSystem, TestNumberSpinachGroupsNoSPiro)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        min_number_spinach_groups: 1
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  // Even though spiro fusions do not extend the ring system, they are not spinach!
  _smiles = "C1C2C3CC34C7(CC7)C421";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_nhits = 0;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_nhits);

  EXPECT_TRUE(_DoPermutationTests(expected_nhits));
}

TEST_F(TestSubstructureRingSystem, TestNumberSpinachGroupsNoMatch)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        number_spinach_groups: 4
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "c1c(F)cc(F)cc1C(=O)Nc1c(C)cc(O)cc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_nhits = 0;  // Not all in the same ring system.

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_nhits);

  EXPECT_TRUE(_DoPermutationTests(expected_nhits));
}

TEST_F(TestSubstructureRingSystem, TestNumberSpinachGroupsMatch)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        number_spinach_groups: 2
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "c1c(F)cc(F)cc1C(=O)Nc1c(C)cc(O)cc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_nhits = 2;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_nhits);

  EXPECT_TRUE(_DoPermutationTests(expected_nhits));
}

TEST_F(TestSubstructureRingSystem, TestNumberNonSpinachGroupsMatch)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        number_non_spinach_groups: 1
      }
      smarts: "C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "c1c(F)cc(F)cc1C(=O)Nc1c(C)cc(O)cc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_nhits = 2;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_nhits);

  EXPECT_TRUE(_DoPermutationTests(expected_nhits));
}

TEST_F(TestSubstructureRingSystem, TestAtomsInSpinachNone)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        atoms_in_spinach_group: 1
      }
      smarts: "c"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "c1cnccc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_nhits = 0;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_nhits);

  EXPECT_TRUE(_DoPermutationTests(expected_nhits));
}

TEST_F(TestSubstructureRingSystem, TestAtomsInSpinachMatch)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        atoms_in_spinach_group: 4
      }
      smarts: "A=A"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "c1cncc(CC(=O)N)c1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_nhits = 2;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_nhits);

  EXPECT_TRUE(_DoPermutationTests(expected_nhits));
}

TEST_F(TestSubstructureRingSystem, TestLengthOfSpinachGroup)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        length_of_spinach_group: 3
      }
      smarts: "A=A"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "c1cncc(CC(=O)N)c1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_nhits = 2;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_nhits);

  EXPECT_TRUE(_DoPermutationTests(expected_nhits));
}

TEST_F(TestSubstructureRingSystem, TestLengthOfSpinachGroupNoMatch)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        length_of_spinach_group: 2
      }
      smarts: "A=A"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "c1cncc(CC(=O)N)c1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_nhits = 0;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_nhits);

  EXPECT_TRUE(_DoPermutationTests(expected_nhits));
}

TEST_F(TestSubstructureRingSystem, TestDistanceToAnotherRingNone)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        distance_to_another_ring: 1
      }
      smarts: "A=A"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC1.C1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_nhits = 0;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_nhits);

  EXPECT_TRUE(_DoPermutationTests(expected_nhits));
}

// Hmmm, this does not work as expected. A reasonable expectation would be
// that you could find a spiro fused ring by specifying a zero value here.
// But that does not work, and it would be messy to fix it.
TEST_F(TestSubstructureRingSystem, TestDistanceToAnotherRingSpiro)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        distance_to_another_ring: 0
      }
      smarts: "A-A"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC12CC2";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_nhits = 0; //  should be 12;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_nhits);

  EXPECT_TRUE(_DoPermutationTests(expected_nhits));
}

TEST_F(TestSubstructureRingSystem, TestDistanceToAnotherRingOne)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        distance_to_another_ring: 1
      }
      smarts: "A-A"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC1C2CC2";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_nhits = 14;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_nhits);

  EXPECT_TRUE(_DoPermutationTests(expected_nhits));
}

TEST_F(TestSubstructureRingSystem, TestDistanceToAnotherRingTwo)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        distance_to_another_ring: 2
      }
      smarts: "A-A"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC1CC2CC2";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_nhits = 16;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_nhits);

  EXPECT_TRUE(_DoPermutationTests(expected_nhits));
}

TEST_F(TestSubstructureRingSystem, TestStronglyFusedNone)
{
  _string_proto = R"(
    query {
      ring_system_specifier {
        strongly_fused_ring_count: 4
      }
      smarts: "A-A"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1C3CC2CC(CC1C2)C3";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  constexpr int expected_nhits = 24;

  ASSERT_EQ(_query.substructure_search(_m, _sresults), expected_nhits);

  EXPECT_TRUE(_DoPermutationTests(expected_nhits));
}

}  // namespace
