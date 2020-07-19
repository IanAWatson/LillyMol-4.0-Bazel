#include <string>
#include <vector>

// Tests for Substructre_Atom_Specifier

#include "aromatic.h"
#include "substructure.h"

#include "googletest/include/gtest/gtest.h"
#include "google/protobuf/text_format.h"
#include "googlemock/include/gmock/gmock.h"

namespace {

using testing::UnorderedElementsAre;

//using google::protobuf::TextFormat::ParseFromString;

class TestSubstructureSpec : public testing::Test
{
  protected:
    void SetUp() override;

  protected:
    std::string _string_proto;

    IWString _smiles;

    Substructure_Query _query;

    Substructure_Results _sresults;

    Molecule _m;

    SubstructureSearch::SubstructureQuery _proto;

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
};

void
TestSubstructureSpec::SetUp()
{
  set_global_aromaticity_type(Daylight);
}

TEST_F(TestSubstructureSpec, AtomicNumber)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  ASSERT_TRUE(_m.build_from_smiles("C"));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);
}

TEST_F(TestSubstructureSpec, AtomicSymbol)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_symbol: "C"
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CN";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1);

  const Set_of_Atoms * e = _sresults.embedding(0);

  EXPECT_EQ(e->number_elements(), 1);
  EXPECT_EQ(_m.atomic_number(e->item(0)), 6);
}

TEST_F(TestSubstructureSpec, Ncon1)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          ncon: 1
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C.NCN.COC.FC(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  for (int i = 0; i < 2; ++i) {
    const Set_of_Atoms * e = _sresults.embedding(i);
    EXPECT_EQ(e->number_elements(), 1);
    EXPECT_EQ(_m.atomic_number(e->item(0)), 6);
    EXPECT_EQ(_m.ncon(e->item(0)), 1);
    EXPECT_EQ(_m.attached_heteroatom_count(e->item(0)), 1);
  }

  const std::vector<int> matched_atoms = {_sresults.embedding(0)->item(0),
                                          _sresults.embedding(1)->item(0)};
  EXPECT_THAT(matched_atoms, UnorderedElementsAre(4, 6));
}

TEST_F(TestSubstructureSpec, MinNcon1)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          min_ncon: 3
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C.NCN.COC.NC(N)N.C(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2);

  const std::vector<int> matched_atoms = {_sresults.embedding(0)->item(0),
                                          _sresults.embedding(1)->item(0)};
  EXPECT_THAT(matched_atoms, UnorderedElementsAre(8, 11));
}

TEST_F(TestSubstructureSpec, Ncon2)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          ncon2: 3
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CC(C(F)(F)F)C(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms* e = _sresults.embedding(0);

  EXPECT_EQ(e->item(0), 0);
}

TEST_F(TestSubstructureSpec, Nbonds)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          nbonds: 3
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CC=N";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms* e = _sresults.embedding(0);

  EXPECT_EQ(e->item(0), 1);
}

TEST_F(TestSubstructureSpec, FormalChargePos)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
          formal_charge: 1
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "[N+]CC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms* e = _sresults.embedding(0);

  EXPECT_EQ(e->item(0), 0);
}

TEST_F(TestSubstructureSpec, FormalChargeNeg)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 8
          formal_charge: -1
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "[O-]C(=O)c1ccccc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms* e = _sresults.embedding(0);

  EXPECT_EQ(e->item(0), 0);
}

TEST_F(TestSubstructureSpec, FormalChargeZero)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
          formal_charge: 0
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "[N+]CN";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms* e = _sresults.embedding(0);

  EXPECT_EQ(e->item(0), 2);
}

TEST_F(TestSubstructureSpec, NringsPositive)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 8
          nrings: 1
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "OC1OC1.O1C2CCC1CC2";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms* e = _sresults.embedding(0);

  EXPECT_EQ(e->item(0), 2);
}

TEST_F(TestSubstructureSpec, RingBondCount)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 8
          ring_bond_count: 2
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "OC1OC1.O1C2CCC1CC2";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 2);
  Set_of_Atoms matched;
  for (const auto* e : _sresults.embeddings()) {
    matched.add(e->item(0));
  }

  EXPECT_THAT(matched, UnorderedElementsAre(2, 4));
}

TEST_F(TestSubstructureSpec, RingSize)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
          ring_size: 5
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "N1CC1.N1CCC1.[1N]1CCCC1.N1CCCCC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms * e = _sresults.embedding(0);

  EXPECT_EQ(e->item(0), 7);
  EXPECT_EQ(_m.isotope(e->item(0)), 1);
}

TEST_F(TestSubstructureSpec, Hcount)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          hcount: 2
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  const Set_of_Atoms * e = _sresults.embedding(0);

  EXPECT_EQ(e->item(0), 1);
}

TEST_F(TestSubstructureSpec, Aromatic)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          aromatic: true
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "c1ncccc1C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 5);
}

TEST_F(TestSubstructureSpec, Chirality)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          chirality: true
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C[C@H](N)CC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  const Set_of_Atoms * e = _sresults.embedding(0);
  EXPECT_EQ(e->item(0), 1);
}

TEST_F(TestSubstructureSpec, AromaticRingSize)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          aromatic_ring_size: 6
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "n1ccc2cnccc12";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 5);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();

  EXPECT_THAT(e, UnorderedElementsAre(3, 4, 6, 7, 8));
}

TEST_F(TestSubstructureSpec, AliphaticRingSize)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          aliphatic_ring_size: 4
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1OC1.C1CCC1.C1CCCC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 4);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();

  EXPECT_THAT(e, UnorderedElementsAre(3, 4, 5, 6));
}

TEST_F(TestSubstructureSpec, AttachedHeteroatomCount)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          attached_heteroatom_count: 3
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C.CC.C(C)C.C(C)(C)C.N[1CH](N)N";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();

  EXPECT_EQ(_m.isotope(e[0]), 1);
}

TEST_F(TestSubstructureSpec, LonePairCount1)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
          lone_pair_count: 1
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NC.CNC.CN(C)C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 3);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(0, 3, 6));
}

TEST_F(TestSubstructureSpec, LonePairCount2)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 8
          lone_pair_count: 2
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "O.OC.COC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 3);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(0, 1, 4));
}

TEST_F(TestSubstructureSpec, Unsaturation1)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 8
          unsaturation: 1
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "O.OC.COC.CC(=O)O";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_EQ(e[0], 8);
}

TEST_F(TestSubstructureSpec, Unsaturation2)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          unsaturation: 2
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CCN.C=C.C#N.C=C=N";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 2);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(5, 8));
}

TEST_F(TestSubstructureSpec, DaylightX1)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          daylight_x: 1
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "[CH].CN.[C]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_EQ(e[0], 0);
}

TEST_F(TestSubstructureSpec, DaylightX2)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          daylight_x: 2
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C[C]C.[CH].CN.[C]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_EQ(e[0], 1);
}

TEST_F(TestSubstructureSpec, DaylightX3)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
          daylight_x: 3
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "N.CN.[N+]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 2);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(0, 2));
}

TEST_F(TestSubstructureSpec, Isotope)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 92
          isotope: [3, 4]
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "[3U].[4U].[U]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 2);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(0, 1));
}

TEST_F(TestSubstructureSpec, Aryl)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
          aromatic: false
          aryl: [1, 2]
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCNc1ccccc1.c1cccnc1Nc1ncccc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 2);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(2, 15));
}


TEST_F(TestSubstructureSpec, Vinyl)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          vinyl: 1
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CNC.COC.CC(=O)C.Cc1ccccc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 2);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(6, 9));
}

TEST_F(TestSubstructureSpec, FusedSystemSize)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
          fused_system_size: 2
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "N1CCC1.N1C2CCC1CC2";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_EQ(e[0], 4);
}

TEST_F(TestSubstructureSpec, HeteroatomsInRing)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
          heteroatoms_in_ring: 2
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "N1CCC1.N1C2CCC1CC2.C1NCNC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 2);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(12, 14));
}

TEST_F(TestSubstructureSpec, SpinachOnlyNoRings1)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
          match_spinach_only: 1
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "N.CN";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 0);
}

TEST_F(TestSubstructureSpec, SpinachOnlyNoRings0)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
          match_spinach_only: 0
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "N.CN.C1CC1NC1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 3);
}

TEST_F(TestSubstructureSpec, SpinachOnly)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 7
          match_spinach_only: 1
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NC1CCN1.N1CCC1.N1C2CCC1CC2.C1NCNC1.C1CC1NC1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_EQ(e[0], 0);
}

TEST_F(TestSubstructureSpec, ScaffoldBondsAttachedToRing1)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          scaffold_bonds_attached_to_ring: 1
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC1CC1CC1.C1CC1CCC1CC1.C1CC1C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 12);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(0, 1, 2, 4, 5, 6, 7, 8, 9, 12, 13, 14));
}

TEST_F(TestSubstructureSpec, ScaffoldBondsAttachedToRing2)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 6
          scaffold_bonds_attached_to_ring: 2
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CC1CC1CC1CCC1CC1F.C1CC1CC1CC1.C1CC1CCC1CC1.C1CC1C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 3);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(4, 5, 6));
}

TEST_F(TestSubstructureSpec, SymmetryDegree1)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 9
          symmetry_degree: 1
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "FCC(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_EQ(e[0], 0);
}

TEST_F(TestSubstructureSpec, SymmetryDegree3)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 9
          symmetry_degree: 3
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "FCC(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 3);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(3,4,5));
}

TEST_F(TestSubstructureSpec, SymmetryDegree6)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 9
          symmetry_degree: 6
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "FC(F)(F)c1ccc(C(F)(F)F)cc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 6);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(0, 2, 3, 9, 10, 11));
}

TEST_F(TestSubstructureSpec, SymmetryDegree3b)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 9
          symmetry_degree: 3
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "FC(F)(F)c1ccc(C(F)(F)F)cc1F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 6);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(0, 2, 3, 9, 10, 11));
}

TEST_F(TestSubstructureSpec, SymmetryGroup1)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number : 9
          symmetry_degree: 3
          symmetry_group: 75
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number : 9
          symmetry_degree: 3
          symmetry_group: 75
        }
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CC(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 6);

  const Set_of_Atoms e = _FirstAtomEachEmbedding();
  EXPECT_THAT(e, UnorderedElementsAre(2, 2, 3, 3, 4, 4));
}

TEST_F(TestSubstructureSpec, TestXor)
{
  _string_proto = R"(query {
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
        }
        atom_properties {
          ncon: 2
          nrings: 0
          logical_operator: SS_XOR
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number : 6
          aromatic: true
        }
        single_bond: 0
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "c1ccc(N)cc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAre(3, 4));

  _smiles = "c1ccc(CC)cc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAre(3, 4));

  // Both parts of the XOR are true.
  _smiles = "c1ccc(NC)cc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 0);
}

TEST_F(TestSubstructureSpec, TestAtomType)
{
  _string_proto = R"(query {
      atom_type: "UST:Y"
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
          atom_type: 6007
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number : 6
          aromatic: true
          atom_type: 3001
        }
        single_bond: 0
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "c1ccc(N)cc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAre(3, 4));
}

TEST_F(TestSubstructureSpec, TestAtomTypeGroupMatches)
{
  _string_proto = R"(query {
      atom_type: "UST:C"
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
          aryl: 1
        }
        atom_type_group: 3
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number : 6
          aryl: 1
        }
        atom_type_group: 3
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "Cc1ccc(N)cc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);

  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAre(0, 5));
}

TEST_F(TestSubstructureSpec, TestAtomTypeGroupNoMatch)
{
  _string_proto = R"(query {
      atom_type: "UST:C"
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
          aromatic: false
          aryl: 1
        }
        atom_type_group: 3
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number : 6
          aromatic: false
          aryl: 1
        }
        atom_type_group: 2
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "Cc1ccc(N)cc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 0);
}

TEST_F(TestSubstructureSpec, TestAtomTypeGroupOrProblem)
{
  _string_proto = R"(query {
      atom_type: "UST:Y"
      query_atom {
        id: 0
        atom_properties {
          atomic_number: [6, 7]
          aromatic: false
        }
        atom_type_group: 3
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number : [6, 7]
          aromatic: false
          aryl: 1
        }
        atom_type_group: 2
        single_bond: 0
      }
    }
  )";

  cerr << "Building query from proto\n";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CNc1ccccc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  if (0 == _query.substructure_search(_m, _sresults))
    cerr << "No matches, hit " << _query.max_query_atoms_matched_in_search() << " atoms\n";
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAre(0, 1));

  _smiles = "NNc1ccccc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 0);

  _smiles = "CCc1ccccc1";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  ASSERT_EQ(_query.substructure_search(_m, _sresults), 0);
}

}  // namespace

