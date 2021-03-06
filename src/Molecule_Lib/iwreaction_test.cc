#include <cmath>

#include <string>

#include "googlemock/include/gmock/gmock.h"
#include "googletest/include/gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "Molecule_Lib/aromatic.h"

#include "iwreaction.h"

namespace {

class TestReaction : public testing::Test
{
  protected:
    void SetUp() override;

  protected:
    std::string _string_proto;

    ReactionProto::Reaction _proto;

    IWReaction _rxn;

    Molecule _m;

    bool _MoleculeFromSmiles(const char * smiles);

    Substructure_Results _sresults;

  protected:
    void _WriteReaction(const char * fname) {
      IWString tmp(fname);
      _rxn.write_msi(tmp);
    }
};

void
TestReaction::SetUp()
{
  set_global_aromaticity_type(Daylight);
  interpret_d_as_deuterium();
}

bool
TestReaction::_MoleculeFromSmiles(const char * smiles)
{
  return _m.build_from_smiles(smiles);
}

TEST_F(TestReaction, TestIsotope)
{
  _string_proto = R"(
    scaffold: {
      id: 0
      smarts: "[CD1]-c"
      isotope {
        atom: 0
        isotope: 4
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("Cc1ccccc1"));

  EXPECT_EQ(_rxn.in_place_transformations(_m), 1);

  EXPECT_EQ(_m.unique_smiles(), "[4CH3]c1ccccc1");
}

TEST_F(TestReaction, TestIncrementIsotope)
{
  _string_proto = R"(
    scaffold: {
      id: 0
      smarts: "[CD1]-c"
      change_isotope {
        atom: 0
        delta: 4
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("[6CH3]c1ccccc1"));

  EXPECT_EQ(_rxn.in_place_transformations(_m), 1);

  EXPECT_EQ(_m.unique_smiles(), "[10CH3]c1ccccc1");
}

TEST_F(TestReaction, TestInvertIsotope)
{
  _string_proto = R"(
    scaffold: {
      id: 0
      smarts: "FCN"
      invert_isotope {
        atom: 0
        isotope: 4
      }
      invert_isotope {
        atom: 2
        isotope: 3
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("[1F]-[2CH2]-N"));

  EXPECT_EQ(_rxn.in_place_transformations(_m), 1);

  EXPECT_EQ(_m.unique_smiles(), "F[2CH2][3NH2]");
}

TEST_F(TestReaction, TestChangeElement)
{
  _string_proto = R"(
    scaffold: {
      id: 0
      smarts: "[Pb]"
      change_element {
        atom: 0
        element: "Au"
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("[Pb].[Pb]"));

  EXPECT_EQ(_rxn.in_place_transformations(_m), 2);

  EXPECT_EQ(_m.unique_smiles(), "[Au].[Au]");
}

TEST_F(TestReaction, TestFormalCharge)
{
  _string_proto = R"(
    scaffold: {
      id: 0
      smarts: "O=[C,S]-[OD1]"
      formal_charge {
        atom: 2
        formal_charge: -1
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("c1ccccc1C(=O)O"));

  EXPECT_EQ(_rxn.in_place_transformations(_m), 1);

  EXPECT_EQ(_m.unique_smiles(), "O=C([O-])c1ccccc1");
}

TEST_F(TestReaction, TestChangeFormalCharge)
{
  _string_proto = R"(
    scaffold: {
      id: 0
      smarts: "O=[C,S]-[OD1-]"
      change_formal_charge {
        atom: 2
        delta: 1
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("c1ccccc1C(=O)[O-]"));

  EXPECT_EQ(_rxn.in_place_transformations(_m), 1);

  EXPECT_EQ(_m.unique_smiles(), "OC(=O)c1ccccc1");
}

TEST_F(TestReaction, TestBreakBond)
{
  _string_proto = R"(
    scaffold: {
      id: 0
      smarts: "O=[C,S]-[ND2]"
      break_bond {
        a1: 1
        a2: 2
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("c1ccccc1CC(=O)NC"));

  EXPECT_EQ(_rxn.in_place_transformations(_m), 1);

  EXPECT_EQ(_m.unique_smiles(), "O=CCc1ccccc1.NC");
}

TEST_F(TestReaction, TestMakeBond)
{
  _string_proto = R"(
    scaffold: {
      id: 0
      smarts: "O=[C,S].N"
      make_bond {
        a1: 1
        a2: 2
        btype: SS_SINGLE_BOND
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("c1ccccc1CC(=O).NC"));

  EXPECT_EQ(_rxn.in_place_transformations(_m), 1);

  EXPECT_EQ(_m.unique_smiles(), "O=C(NC)Cc1ccccc1");
}

TEST_F(TestReaction, TestRemoveAtom)
{
  _string_proto = R"(
    scaffold: {
      id: 0
      smarts: "N.C"
      remove_atom: [0, 1]
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("Nc1c(C)cccc1"));

  EXPECT_EQ(_rxn.in_place_transformations(_m), 1);

  EXPECT_EQ(_m.unique_smiles(), "c1ccccc1");
}

TEST_F(TestReaction, TestRemoveFragment)
{
  _string_proto = R"(
    scaffold: {
      id: 0
      smarts: "c!@-*"
      break_bond {
        a1: 0
        a2: 1
      }
      remove_fragment: 1
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("CNc1c(CC)cccc1"));

  EXPECT_EQ(_rxn.in_place_transformations(_m), 2);

  EXPECT_EQ(_m.unique_smiles(), "c1ccccc1");
}

// This one is counterintuitive. The intent is to remove the aromatic ring, and retain
// the substituents. BUT, the substructure hits get processed independently, and
// the first match will eliminate the ring and everything attached to it.
TEST_F(TestReaction, TestKeepFragment)
{
  _string_proto = R"(
    scaffold: {
      id: 0
      smarts: "c!@-*"
      break_bond {
        a1: 0
        a2: 1
      }
      keep_fragment: 1
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("CNc1c(CC)cccc1"));

  EXPECT_EQ(_rxn.in_place_transformations(_m), 2);

  EXPECT_EQ(_m.unique_smiles(), "CC");
}

TEST_F(TestReaction, TestMakeSingleBond)
{
  _string_proto = R"(
    scaffold: {
      id: 0
      smarts: "[CD1].[CD1]"
      make_bond {
        a1: 0
        a2: 1
        btype: SS_SINGLE_BOND
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("CCCC"));

  EXPECT_EQ(_rxn.in_place_transformations(_m), 2);

  EXPECT_EQ(_m.unique_smiles(), "C1CCC1");
}

TEST_F(TestReaction, TestMakedoubleBond)
{
  _string_proto = R"(
    scaffold: {
      id: 0
      smarts: "[CD1].[CD1]"
      make_bond {
        a1: 0
        a2: 1
        btype: SS_DOUBLE_BOND
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("CCCC"));

  EXPECT_EQ(_rxn.in_place_transformations(_m), 2);

  EXPECT_EQ(_m.unique_smiles(), "C1CC=C1");
}

TEST_F(TestReaction, TestMakeTripleBond)
{
  _string_proto = R"(
    scaffold: {
      id: 0
      smarts: "[CD1].[ND0]"
      make_bond {
        a1: 0
        a2: 1
        btype: SS_TRIPLE_BOND
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("N.CC(=O)O"));

  EXPECT_EQ(_rxn.in_place_transformations(_m), 1);

  EXPECT_EQ(_m.unique_smiles(), "OC(=O)C#N");
}

TEST_F(TestReaction, TestBondLength)
{
  _string_proto = R"(
    scaffold: {
      id: 0
      smarts: "[CD1]-[ND1]"
      bond_length {
        a1: 0
        a2: 1
        distance: 1.4
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("C{{0,0,0}}N{{1,0,0}}"));

  EXPECT_EQ(_rxn.in_place_transformations(_m), 1);

  EXPECT_EQ(_m.unique_smiles(), "NC");

  EXPECT_FLOAT_EQ(_m.bond_length(0, 1), 1.4);
}

TEST_F(TestReaction, TestBondAngle)
{
  _string_proto = R"(
    scaffold: {
      id: 0
      smarts: "FCN"
      bond_angle {
        a1: 0
        a2: 1
        a3: 2
        angle: 109
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  const float angle = 109.0 / M_PI/ 180.0;

  _proto.mutable_scaffold()->mutable_bond_angle(0)->set_angle(angle);

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("F{{0,0,0}}C{{1,0,0}}N{{2,0,0}}"));

  EXPECT_EQ(_rxn.in_place_transformations(_m), 1);

  EXPECT_EQ(_m.unique_smiles(), "FCN");

  EXPECT_NEAR(_m.bond_angle(0, 1, 2), angle, 1.0e-06);
}

TEST_F(TestReaction, TestDihedralAngle)
{
  _string_proto = R"(
    scaffold: {
      id: 0
      smarts: "CCCN"
      dihedral_angle {
        a1: 0
        a2: 1
        a3: 2
        a4: 3
        angle: 0.5
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  const float angle = 0.5;

  _proto.mutable_scaffold()->mutable_dihedral_angle(0)->set_angle(angle);

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("C{{-1,-1,0}}C{{-0.5,0,0}}C{{0.5,0,0}}N{{1,1,0}}"));

  EXPECT_EQ(_rxn.in_place_transformations(_m), 1);

  EXPECT_EQ(_m.unique_smiles(), "NCCC");

  EXPECT_NEAR(fabs(_m.dihedral_angle(0, 1, 2, 3)), angle, 1.0e-06);
}

// Use of an Substructure_Query is complex. It is
// embedded inside the ReactionSite, and then a Substructure_Query
// consists of multiple Single_Substructure_Query objects.
// That is why we get to so many levels here.
TEST_F(TestReaction, TestSubstructureQuery)
{
  _string_proto = R"(
    scaffold: {
      id: 0
      query {
        query {
          smarts: "[CD1]-C"
          max_matches_to_find: 1
        }
      }
      break_bond {
        a1: 0
        a2: 1
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("CCNCCC"));

  ASSERT_EQ(_rxn.in_place_transformations(_m), 1);

  EXPECT_EQ(_m.unique_smiles(), "C.CNCCC");
}

TEST_F(TestReaction, TestReplaceAtom)
{
  _string_proto = R"(
    scaffold: {
      id: 0
      smarts: "NC(C)F.Cl"
      replace_atom {
        a1: 3
        a2: 4
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("N[C@H](C)F.Cl"));

  ASSERT_EQ(_rxn.in_place_transformations(_m), 1);

  EXPECT_EQ(_m.unique_smiles(), "F.Cl[C@H](N)C");
}

TEST_F(TestReaction, TestInvertChirality)
{
  _string_proto = R"(
    scaffold: {
      id: 0
      smarts: "NC(C)F"
      invert_chirality: 1
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("F[C@H](N)C"));

  ASSERT_EQ(_rxn.in_place_transformations(_m), 1);

  EXPECT_EQ(_m.unique_smiles(), "F[C@@H](N)C");
}

TEST_F(TestReaction, TestRemoveChirality)
{
  _string_proto = R"(
    scaffold: {
      id: 0
      smarts: "NC(C)F"
      remove_chirality: 1
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("F[C@H](N)C"));

  ASSERT_EQ(_rxn.in_place_transformations(_m), 1);

  EXPECT_EQ(_m.unique_smiles(), "FC(N)C");
}

TEST_F(TestReaction, TestFixedReagent)
{
  _string_proto = R"(
    scaffold: {
      id: 0
      smarts: "[OD1]-C=O"
      remove_atom: 0
    }
    sidechain {
      id: 1
      smarts: "[ND1]"
      reagent: "NCC"
      join {
        a1: 1
        a2: 0
        btype: SS_SINGLE_BOND
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("CC(=O)O"));

  EXPECT_EQ(_rxn.substructure_search(&_m, _sresults), 1);

  Reaction_Iterator iterator;
  for (iterator.initialise(_rxn); iterator.active(); iterator++)
  {
    Molecule result;
    ASSERT_EQ(_rxn.perform_reaction(&_m, _sresults, iterator, result), 1);
    ASSERT_EQ(result.smiles(), "CC(=O)NCC");
  }
}

TEST_F(TestReaction, TestMultipleReagents)
{
  _string_proto = R"(
    scaffold: {
      id: 0
      smarts: "[OD1]-C=O"
      remove_atom: 0
    }
    sidechain {
      id: 1
      smarts: "[ND1]"
      reagent: "NC"
      reagent: "NCC"
      reagent: "NCCC"
      join {
        a1: 1
        a2: 0
        btype: SS_SINGLE_BOND
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("CC(=O)O"));

  EXPECT_EQ(_rxn.substructure_search(&_m, _sresults), 1);

  const std::vector<IWString> expected = {"CC(=O)NC", "CC(=O)NCC", "CC(=O)NCCC"};

  Reaction_Iterator iterator;
  int ndx = 0;
  for (iterator.initialise(_rxn); iterator.active(); iterator++, ++ndx)
  {
    Molecule result;
    ASSERT_EQ(_rxn.perform_reaction(&_m, _sresults, iterator, result), 1);
    ASSERT_EQ(result.smiles(), expected[ndx]);
  }
}

TEST_F(TestReaction, TestAppendReagentName) {
  _string_proto = R"(
    append_reagent_name: true
    scaffold: {
      id: 0
      smarts: "[CD1]"
    }
    sidechain {
      id: 1
      smarts: "C"
      reagent: "C carbon"
      join {
        a1: 0
        a2: 0
        btype: SS_SINGLE_BOND
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("CC ethane"));

  EXPECT_EQ(_rxn.substructure_search(&_m, _sresults), 2);

  Reaction_Iterator iterator;
  for (iterator.initialise(_rxn); iterator.active(); iterator++)
  {
    Molecule result;
    ASSERT_EQ(_rxn.perform_reaction(&_m, _sresults, iterator, result), 1);
    EXPECT_EQ(result.smiles(), "C(C)CC");
    EXPECT_EQ(result.name(), "ethane + carbon + carbon");
  }
}

TEST_F(TestReaction, TestAppendReagentNameEnumerate) {
  _string_proto = R"(
    append_reagent_name: true
    scaffold: {
      id: 0
      smarts: "[CD1]"
    }
    sidechain {
      id: 1
      smarts: "C"
      reagent: "C carbon"
      join {
        a1: 0
        a2: 0
        btype: SS_SINGLE_BOND
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("CC ethane"));

  EXPECT_EQ(_rxn.substructure_search(&_m, _sresults), 2);

  Reaction_Iterator iterator;
  int ndx = 0;
  for (iterator.initialise(_rxn); iterator.active(); iterator++, ++ndx)
  {
    Molecule result;
    ASSERT_EQ(_rxn.perform_reaction(&_m, _sresults.embedding(ndx), iterator, result), 1);
    EXPECT_EQ(result.smiles(), "C(C)C");
    EXPECT_EQ(result.name(), "ethane + carbon");
  }
}

TEST_F(TestReaction, TestAppendToName) {
  _string_proto = R"(
    append_to_name: "hello"
    scaffold: {
      id: 0
      smarts: "[CD1]"
      change_element {
        atom: 0
        element: "N"
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  ASSERT_TRUE(_rxn.ConstructFromProto(_proto));

  ASSERT_TRUE(_MoleculeFromSmiles("CC1CC1 name"));

  EXPECT_EQ(_rxn.substructure_search(&_m, _sresults), 1);

  EXPECT_EQ(_rxn.in_place_transformations(_m, _sresults), 1);

  EXPECT_EQ(_m.smiles(), "NC1CC1");
  EXPECT_EQ(_m.name(), "namehello");
}

}  // namespace
