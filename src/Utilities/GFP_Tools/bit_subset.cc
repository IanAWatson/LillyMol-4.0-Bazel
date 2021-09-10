#include "Foundational/iwmisc/proto_support.h"

#include "Utilities/GFP_Tools/bit_subset.h"

namespace bit_subset {

using std::cerr;

BitSubset::BitSubset() {
  _flatten_sparse_counted = 0;

  _nfixed = 0;
  _fixed_mask = nullptr;
  _nsparse = 0;
  _sparse_subset = nullptr;

}

BitSubset::~BitSubset() {
  if (_fixed_mask != nullptr) {
    delete [] _fixed_mask;
    _nfixed = 0;
  }

  if (_sparse_subset != nullptr) {
    delete [] _sparse_subset;
    _nsparse = 0;
  }
}

void
FormMask(std::unordered_set<uint32_t>& subset,
         const int nbits,
         IWDYFP& mask) {
  mask.allocate_space_for_bits(nbits);
  mask.clear();
  for (const auto iter : subset) {
    mask.set(iter);
  }
}

int
BitSubset::InitialiseGfpKnown(const IW_General_Fingerprint& gfp) {
  _nfixed = number_fingerprints();
  _nsparse = number_sparse_fingerprints();
  if (_nfixed == 0 && _nsparse == 0) {
    cerr << "BitSubset::InitialiseGfpKnown:no fingerprints\n";
    return 0;
  }

  _fixed_mask = new IWDYFP[_nfixed];
  _sparse_subset = new Subset*[_nsparse];

  for (int i = 0; i < _nfixed; ++i) {
    const IWString& tag = fixed_fingerprint_tag(i);
    const auto iter =  _tag_to_bit_subset.find(tag);
    if (iter == _tag_to_bit_subset.end()) {
      cerr << "BitSubset::InitialiseGfpKnown:no fixed '" << tag << "' tag, nfixed " << _nfixed << '\n';
      return 0;
    }
    const IWDYFP& fp = gfp[i];
    FormMask(iter->second, fp.nbits(), _fixed_mask[i]);
  }

  for (int i = 0; i < _nsparse; ++i) {
    const IWString& tag = sparse_fingerprint_tag(i);
    auto iter = _tag_to_bit_subset.find(tag);
    if (iter == _tag_to_bit_subset.end()) {
      cerr << "BitSubset::InitialiseGfpKnown:no sparse '" << tag << "' tag, nsparse " << _nsparse << '\n';
      return 0;
    }
    _sparse_subset[i] = &(iter->second);
  }
  
  return 1;
}

int
BitSubset::Build(IWString& fname) {
  std::optional<GfpBitSubset::GfpBitSubset> proto =
                   iwmisc::ReadBinaryProto<GfpBitSubset::GfpBitSubset>(fname);
  if (! proto) {
    cerr << "BitSubset::Build:cannot read '" << fname << "'\n";
    return 0;
  }

  return Build(*proto);
}

int
BitSubset::Build(const GfpBitSubset::GfpBitSubset& proto) {
  for (const auto& [key, value] : proto.xref()) {
    IWString tag(key);
    Subset subset;

    for (const auto iter : value.bits()) {
      subset.insert(iter);
    }
    _tag_to_bit_subset.emplace(tag, std::move(subset));
  }

  _flatten_sparse_counted = proto.params().flatten_sparse_fingerprints();

  return _tag_to_bit_subset.size();
}

BitXref::BitXref() {
  _flatten_sparse_counted = 0;

  _nfixed = 0;
  _fixed = nullptr;
  _nsparse = 0;
  _sparse = nullptr;
}

BitXref::~BitXref() {
  if (_fixed != nullptr) {
    delete [] _fixed;
    _nfixed = 0;
  }

  if (_sparse != nullptr) {
    delete [] _sparse;
    _nsparse = 0;
  }
}

int
BitXref::Build(IWString& fname) {
  std::optional<GfpBitSubset::GfpBitToFeature> proto = iwmisc::ReadBinaryProto<GfpBitSubset::GfpBitToFeature>(fname);
  if (! proto) {
    cerr << "BitXref::Build:cannot read '" << fname << "'\n";
    return 0;
  }

  return Build(*proto);
}

int
BitXref::Build(const GfpBitSubset::GfpBitToFeature& proto) {
  for (const auto& [key, value] : proto.xref()) {
    IWString tag(key);
    BitToFeature bit_to_feature;

    for (auto [key, value] : value.bit_to_feature()) {
      bit_to_feature.emplace(static_cast<gfp_bit_type>(key), static_cast<uint32_t>(value));
    } 
    _tag_to_bit_to_feature.emplace(tag, std::move(bit_to_feature));
  }

  _flatten_sparse_counted = proto.params().flatten_sparse_fingerprints();
  return 1;
}

struct FeatureCount {
  uint32_t feature;
  uint32_t count;
};

// Return the number of features set in `gfp`.
int
NumberFeatures(const IW_General_Fingerprint& gfp) {
  int nfixed = number_fingerprints();
  int nsparse = number_sparse_fingerprints();
  int nfeatures = 0;
  for (int i = 0; i < nfixed; ++i) {
    nfeatures += gfp[i].nset();
  }
  for (int i = 0; i < nsparse; ++i) {
    nfeatures += gfp.sparse_fingerprint(i).nset();
  }

  return nfeatures;
}

int
BitSubset::MakeSubset(IW_General_Fingerprint& gfp) {
  if (_nfixed == 0 && _nsparse == 0 && ! InitialiseGfpKnown(gfp)) {
    cerr << "BitSubset::MakeSubset:cannot initialise\n";
    return -1;  // Special return value indicating fatal error.
  }

  int rc = 0;
  for (int i = 0; i < _nfixed; ++i) {
    IWDYFP& fp = gfp[i];
    const int initial_nset = fp.nset();
    int notused = 0;
    fp.iwand(_fixed_mask[i], notused);
    rc += initial_nset - fp.nset();
  }

  for (int i = 0; i < _nsparse; ++i) {
    Sparse_Fingerprint & sparse_fingerprint = gfp.sparse_fingerprint(i);

    rc += sparse_fingerprint.ReduceToSubset(*_sparse_subset[i]);
  }

  return rc;
}


int
BitXref::WriteSvmlFeatures(const IW_General_Fingerprint& gfp,
                        IWString_and_File_Descriptor& output) {
  if (_nfixed == 0 && _nsparse == 0 && ! InitialiseGfpKnown(gfp)) {
    cerr << "BitXref::WriteSvmlFeatures:cannot initialise\n";
    return -1;  // Special return value indicating fatal error.
  }

  // The maximum number of features.
  const int nfeatures = NumberFeatures(gfp);

  int ndx = 0;
  std::unique_ptr<FeatureCount[]> feature_counts(new FeatureCount[nfeatures]);

  for (int i = 0; i < _nfixed; ++i) {
    const BitToFeature* b2f = _fixed[i];
    const IWDYFP & fp = gfp[i];
    int j = 0;
    while (int b = fp.next_on_bit(j)) {
      const auto iter = b2f->find(b);
      if (iter == b2f->cend()) {
        continue;
      }
      feature_counts[ndx].feature = iter->second;
      feature_counts[ndx].count = 1;
      ndx++;
    }
  }

  for (int i = 0; i < _nsparse; ++i) {
    const BitToFeature* b2f = _sparse[i];
    const Sparse_Fingerprint& fp = gfp.sparse_fingerprint(i);
    int j = 0;
    uint32_t b = 0;
    int c = 0;
    while (fp.next_bit_set(j, b, c)) {
      const auto iter = b2f->find(b);
      if (iter == b2f->end()) {
        continue;
      }
      feature_counts[ndx].feature = iter->second;
      feature_counts[ndx].count = 1;
      ndx++;
    }
  }

  std::sort(feature_counts.get(), feature_counts.get() + ndx,
            [](const FeatureCount& fc1, const FeatureCount& fc2) {
    return fc1.feature < fc2.feature;
  });

  for (int i = 0; i < ndx; ++i) {
    output << ' ' << feature_counts[i].feature << ':' << feature_counts[i].count;
  }

  return ndx;  // The number of features written.
}

int
BitXref::InitialiseGfpKnown(const IW_General_Fingerprint& gfp) {
  _nfixed = number_fingerprints();
  _nsparse = number_sparse_fingerprints();
  _fixed = new BitToFeature*[_nfixed];
  _sparse = new BitToFeature*[_nsparse];
  for (int i = 0; i < _nfixed; ++i) {
    const IWString tag = fixed_fingerprint_tag(i);
    const auto iter = _tag_to_bit_to_feature.find(tag);
    if (iter == _tag_to_bit_to_feature.cend()) {
      cerr << "BitXref::InitialiseGfpKnown:no fixed '" << tag << "' fingerprint in input\n";
      return 0;
    }
    _fixed[i] = &(iter->second);
  }

  for (int i = 0; i < _nsparse; ++i) {
    const IWString tag = sparse_fingerprint_tag(i);
    const auto iter = _tag_to_bit_to_feature.find(tag);
    if (iter == _tag_to_bit_to_feature.cend()) {
      cerr << "BitXref::InitialiseGfpKnown:no sparse '" << tag << "' fingerprint in input\n";
      for (auto& [k, v] : _tag_to_bit_to_feature) {
        cerr << "key '" << k << "'\n";
      }
      return 0;
    }
    _sparse[i] = &(iter->second);
  }

  return 1;
}

}  // namespace bit_subset
