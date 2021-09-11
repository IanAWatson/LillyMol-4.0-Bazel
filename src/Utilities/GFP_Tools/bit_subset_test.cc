// Tests for the bit_subset namespace
// Testing is complicated by the fact that gfp depends on
// some globally hidden state - number of fixed/sparse
// fingerprints and tags. For that reason, we use a single
// gfp here.

#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

#include "googlemock/include/gmock/gmock.h"
#include "googletest/include/gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "Utilities/GFP_Tools/bit_subset.h"
#include "Utilities/GFP_Tools/gfp_to_svm_lite.pb.h"

namespace {

IWString fixed_tag = "FPD<";
IWString sparse_tag = "NCS<";
IWString smiles_tag = "$SMI<";
IWString identifier_tag = "PCN<";

class TestBitSubset : public testing::Test
{
  protected:
    void SetUp() override;

  protected:
    IW_General_Fingerprint _gfp;

    std::optional<IW_General_Fingerprint> BuildGfp(const std::vector<uint32_t>& fixed,
                  const std::unordered_map<uint32_t, uint32_t>& sparse);
};

void
TestBitSubset::SetUp() {
}

std::optional<IW_General_Fingerprint>
TestBitSubset::BuildGfp(const std::vector<uint32_t>& fixed,
         const std::unordered_map<uint32_t, uint32_t>& sparse) {
  const int nbits = fixed.size();
  IW_Bits_Base bits(nbits);
  for (int i = 0; i < nbits; ++i) {
    if (fixed[i]) {
      bits.set(i);
    }
  }
  IWString ascii_fixed;
  bits.daylight_ascii_representation_including_nset_info(ascii_fixed);

  Sparse_Fingerprint_Creator sfc;
  for (const auto [key, value] : sparse) {
    sfc.hit_bit(key, value);
  }

  Sparse_Fingerprint sfp;
  if (! sfp.build_from_sparse_fingerprint_creator(sfc)) {
    cerr << "BuildGfp:build_from_sparse_fingerprint_creator failed\n";
    return std::nullopt;
  }

  IWString ascii_sparse;
  sfp.append_daylight_ascii_form_with_counts_encoded(ascii_sparse);

  IWString tdt;
  tdt << smiles_tag << "C>\n";
  tdt << identifier_tag << "C>\n";
  tdt << fixed_tag << ascii_fixed << ">\n";
  tdt << sparse_tag << ascii_sparse << ">\n";
  tdt << "|\n";

  IW_TDT iw_tdt;
  if (! iw_tdt.Build(tdt)) {
    cerr << "Cannot build tdt from '" << tdt << "'\n";
    return std::nullopt;
  }

  IW_General_Fingerprint gfp;
  int fatal = 0;
  if (! gfp.construct_from_tdt(iw_tdt, fatal)) {
    cerr << "BuildGfp:cannot build from tdt\n";
    return std::nullopt;
  }

  return gfp;
}

TEST_F(TestBitSubset, TestGfpFormationEmpty) {
  constexpr int nbits = 64;
  std::vector<uint32_t> bits(nbits);
  bits.assign(nbits, 0);
  std::unordered_map<uint32_t, uint32_t> sparse;
  std::optional<IW_General_Fingerprint> gfp = BuildGfp(bits, sparse);
  ASSERT_TRUE(gfp);
  IWDYFP& fp = (*gfp)[0];
  EXPECT_EQ(fp.nbits(), nbits);
  EXPECT_EQ(fp.nset(), 0);
  const Sparse_Fingerprint& sfp = gfp->sparse_fingerprint(0);
  EXPECT_EQ(sfp.nbits(), 0);
}

TEST_F(TestBitSubset, TestGfpFormation1Bit) {
  constexpr int nbits = 64;
  std::vector<uint32_t> bits(nbits);
  bits.assign(nbits, 0);
  constexpr int set_bit = 53;  // Will be set in both the dense and sparse fp.
  bits[set_bit] = 1;
  std::unordered_map<uint32_t, uint32_t> sparse;
  constexpr int count = 3;
  sparse[set_bit] = count;
  std::optional<IW_General_Fingerprint> gfp = BuildGfp(bits, sparse);
  ASSERT_TRUE(gfp);
  IWDYFP& fp = (*gfp)[0];
  EXPECT_EQ(fp.nbits(), nbits);
  EXPECT_EQ(fp.nset(), 1);
  EXPECT_TRUE(fp.is_set(set_bit));
  const Sparse_Fingerprint& sfp = gfp->sparse_fingerprint(0);
  EXPECT_EQ(sfp.nbits(), 1);
  EXPECT_EQ(sfp.nset(), count);
  EXPECT_EQ(sfp.count_for_bit(set_bit), count);
}

void
SetBits(const std::vector<uint32_t>& on_bits,
        std::vector<uint32_t>& destination) {
  for (auto b : on_bits) {
    destination[b] = 1;
  }
}

TEST_F(TestBitSubset, TestGfpFormationSeveralBits) {
  constexpr int nbits = 64;
  std::vector<uint32_t> bits(nbits);
  bits.assign(nbits, 0);
  const std::vector<uint32_t> set_bits{0, 5, 9, 19, 21, 22, 23, 41, 58, 59, 61, 63};
  const std::vector<uint32_t> counts  {1, 1, 2,  3, 21, 11, 12,  1,  8,  3, 61, 63};

  SetBits(set_bits, bits);
  for (uint32_t b : set_bits) {
    bits[b] = 1;
  }

  std::unordered_map<uint32_t, uint32_t> sparse;
  for (uint32_t i = 0; i < set_bits.size(); ++i) {
    sparse[set_bits[i]] = counts[i];
  }
  std::optional<IW_General_Fingerprint> gfp = BuildGfp(bits, sparse);
  ASSERT_TRUE(gfp);
  IWDYFP& fp = (*gfp)[0];
  EXPECT_EQ(fp.nbits(), nbits);
  EXPECT_EQ(fp.nset(), set_bits.size());
  for (const auto b : set_bits) {
    EXPECT_TRUE(fp.is_set(b));
  }
  const Sparse_Fingerprint& sfp = gfp->sparse_fingerprint(0);
  EXPECT_EQ(sfp.nbits(), set_bits.size());
  EXPECT_EQ(sfp.nset(), std::accumulate(counts.cbegin(), counts.cend(), 0));
  for (uint32_t i = 0; i < set_bits.size(); ++i) {
    EXPECT_EQ(sfp.count_for_bit(set_bits[i]), counts[i]);
  }
}

// Generate an arbitrary sparse fingerprint.
Sparse_Fingerprint
ASparseFingerprint(const int nbits,
                   uint32_t bit,
                   const uint32_t bit_delta,
                   int count,
                   const int count_delta) {
  Sparse_Fingerprint_Creator sfc;
  for (int i = 0; i < nbits; ++i) {
    sfc.hit_bit(bit, count);
    bit += bit_delta;
    count += count_delta;
    if (count > std::numeric_limits<uint8_t>::max()) {
      count = 10;   // An arbitrary choice.
    }
  }

  Sparse_Fingerprint result;
  result.build_from_sparse_fingerprint_creator(sfc);
  return result;
}

TEST_F(TestBitSubset, NothingSpecifiedFailsToBuild) {
  constexpr int nbits = 77;
  const Sparse_Fingerprint sfp = ASparseFingerprint(nbits, 1099, 1056, 5, 3);
  ASSERT_EQ(sfp.nbits(), nbits);

  const std::string proto_text = R"pb(
    params {
    }
  )pb";
  GfpBitSubset::GfpBitSubset proto;
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(proto_text, &proto));
  bit_subset::BitSubset subset;
  EXPECT_FALSE(subset.Build(proto));
}

// No overlap between what is retained and what is in the starting fingerprint.
TEST_F(TestBitSubset, AllBitsRemoved) {
  std::vector<uint32_t> fixed(64);
  const std::vector<uint32_t> keep_fixed {0, 10, 20, 30, 40, 50, 60};
  SetBits(keep_fixed, fixed);
  std::unordered_map<uint32_t, uint32_t> sparse { {1, 1}, {2, 22}, {30, 1}, {55, 55}, {1024, 1}};

  std::optional<IW_General_Fingerprint> gfp = BuildGfp(fixed, sparse);
  ASSERT_TRUE(gfp);

  // No overlap between bits set and bits to be retained.
  const std::string proto_text = R"pb(
params {
}
xref {
  key: "FPD"
  value {
    bits: 55
    bits: 61
  }
}
xref {
  key: "NCS"
  value {
  }
}
)pb";

  GfpBitSubset::GfpBitSubset proto;

#ifdef FOR_DEBUGGING
  std::unordered_map<uint32_t, uint32_t> xref;
  xref[3] = 1;
  xref[5] = 1;

  std::string tag("FPD<");
  GfpBitSubset::BitSubset my_subset;
  my_subset.mutable_bits()->Add(55);
  my_subset.mutable_bits()->Add(58);
  google::protobuf::MapPair<std::string, GfpBitSubset::BitSubset> to_insert(tag, my_subset);
  proto.mutable_xref()->insert(to_insert);
  tag = "NCS<";
  google::protobuf::MapPair<std::string, GfpBitSubset::BitSubset> to_insert2(tag, my_subset);
  proto.mutable_xref()->insert(to_insert2);
  IWString fname("/tmp/a.dat");
  iwmisc::WriteBinaryProto(proto, fname);
#endif

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(proto_text, &proto));

  bit_subset::BitSubset bit_subset;
  ASSERT_TRUE(bit_subset.Build(proto));
  ASSERT_TRUE(bit_subset.InitialiseGfpKnown(*gfp));

  // ALl bits have been removed.
  EXPECT_EQ(bit_subset.MakeSubset(*gfp), keep_fixed.size() + sparse.size());
  IWDYFP& fp = (*gfp)[0];
  EXPECT_EQ(fp.nset(), 0);
  const Sparse_Fingerprint& sfp = gfp->sparse_fingerprint(0);
  EXPECT_EQ(sfp.nset(), 0);
}

TEST_F(TestBitSubset, SomeBitsRemoved) {
  std::vector<uint32_t> fixed(64);
  const std::vector<uint32_t> keep_fixed {0, 10, 20, 30, 40, 50, 60};
  SetBits(keep_fixed, fixed);
  std::unordered_map<uint32_t, uint32_t> sparse { {1, 1}, {2, 22}, {30, 1}, {55, 55}, {1024, 1}};

  std::optional<IW_General_Fingerprint> gfp = BuildGfp(fixed, sparse);
  ASSERT_TRUE(gfp);

  // Retain some of the bits from both the sparse and dense fingerprint.
  const std::string proto_text = R"pb(
params {
}
xref {
  key: "FPD"
  value {
    bits: 10
    bits: 20
    bits: 60
  }
}
xref {
  key: "NCS"
  value {
    bits: 2
    bits: 30
    bits: 1024
  }
}
)pb";

  GfpBitSubset::GfpBitSubset proto;

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(proto_text, &proto));

  bit_subset::BitSubset bit_subset;
  ASSERT_TRUE(bit_subset.Build(proto));
  ASSERT_TRUE(bit_subset.InitialiseGfpKnown(*gfp));

  EXPECT_EQ(bit_subset.MakeSubset(*gfp), 6);
  IWDYFP& fp = (*gfp)[0];
  EXPECT_EQ(fp.nset(), 3);
  auto iter = proto.xref().find("FPD");
  ASSERT_NE(iter, proto.xref().end());
  for (const auto b : iter->second.bits()) {
    EXPECT_TRUE(fp.is_set(b));
  }

  const Sparse_Fingerprint& sfp = gfp->sparse_fingerprint(0);
  EXPECT_EQ(sfp.nbits(), 3);

  iter = proto.xref().find("NCS");
  ASSERT_NE(iter, proto.xref().end());
  for (const auto b : iter->second.bits()) {
    EXPECT_TRUE(sfp.is_set(b));
    EXPECT_EQ(sfp.count_for_bit(b), sparse[b]);
  }
}
}  // namespace