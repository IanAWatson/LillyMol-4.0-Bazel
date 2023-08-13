// Compute modal fingerprint from a collection

#include <algorithm>
#include <iostream>
#include <memory>

#include "absl/container/flat_hash_map.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "gfp.h"

namespace gfp_model_fingerprint {

using std::cerr;

void
Usage(int rc) {
  cerr << "Compute modal fingerprint from a collection\n";
  cerr << " -p fixed    compute bit pairs as well as individual bits in fixed width fingerprints\n";
  cerr << " -p sparse   compute bit pairs as well as individual bits in sparse fingerprints\n";
  cerr << " -v          verbose output\n";

  ::exit(rc);
}

class Options {
  private:
    int _verbose;

    resizable_array_p<IW_General_Fingerprint> _pool;

    int _do_fixed_fp_pairwise;
    int _do_sparse_fp_pairwise;

    int _log10_counts;

  // private functions
    int ReadFingerprints(const char* fname);
    int ReadFingerprints(iwstring_data_source& input);

    int FixedFingerprint(int fingerprint, IWString_and_File_Descriptor& output);
    void FixedFingerprint(const IWDYFP& fp, uint32_t* count) const;
    void Normalise(uint32_t* count, uint32_t nbits) const;
    int PairwiseFixedFingerprint(int fingerprint,
        uint32_t* count,
        uint32_t nbits,
        IWString_and_File_Descriptor& output) const;

    void PairwiseFixedFingerprint(
        const IWDYFP& fp,
        uint32_t nbits,
        uint32_t* count,
        absl::flat_hash_map<uint32_t, uint32_t>& pairs) const;
    void PairwiseFixedFingerprint(
        const IWDYFP& fp,
        uint32_t nbits,
        uint32_t* count,
        uint32_t* pairs) const;

    int SparseFingerprint(int fingerprint, IWString_and_File_Descriptor& output);
    void SparseFingerprint(const Sparse_Fingerprint& fp,
                absl::flat_hash_map<uint32_t, uint32_t> & count) const;
    int PairwiseSparseFingerprint(int fingerprint, IWString_and_File_Descriptor& output);
    void PairwiseSparseFingerprint(const Sparse_Fingerprint& fp, absl::flat_hash_map<uint32_t, uint32_t>& count);
    int PairwiseSparseFingerprint(int fingerprint,
                const absl::flat_hash_map<uint32_t, uint32_t>& count,
                IWString_and_File_Descriptor& output);

    void HashToFingerprint(const absl::flat_hash_map<uint32_t, uint32_t>& count,
                        Sparse_Fingerprint_Creator& sfc);
    void HashToFingerprintLog10(const absl::flat_hash_map<uint32_t, uint32_t>& count,
                        Sparse_Fingerprint_Creator& sfc) const;
    void HashToFingerprintLinear(const absl::flat_hash_map<uint32_t, uint32_t>& count,
                        Sparse_Fingerprint_Creator& sfc) const;
  public:
    Options();

    int Initialise(Command_Line& cl);

    int Process(IWString_and_File_Descriptor& output);
};

Options::Options() {
  _do_fixed_fp_pairwise = 0;
  _do_sparse_fp_pairwise = 0;

  _log10_counts = 0;
}

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_present('v');

  if (cl.option_present('p')) {
    const_IWSubstring p;
    for (int i = 0; cl.value('p', p, i); ++i) {
      if (p == "fixed") {
        _do_fixed_fp_pairwise = 1;
      } else if (p == "sparse") {
        _do_sparse_fp_pairwise = 1;
      } else {
        cerr << "Unrecognised -p qualifier '" << p << "'\n";
        return 0;
      }
    }
  }

  if (cl.option_present('l')) {
    _log10_counts = 1;
    if (_verbose) {
      cerr << "Will use log10 of count data\n";
    }
  }

  for (const char* fname: cl) {
    if (! ReadFingerprints(fname)) {
      cerr << "Options::Initialise:cannot read '" << fname << "'\n";
      return 0;
    }
  }

  if (_verbose) {
    cerr << "Read " << _pool.size() << " fingerprints\n";
  }

  return 1;
}

int
Options::ReadFingerprints(const char* fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::ReadFingerprints:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadFingerprints(input);
}

int
Options::ReadFingerprints(iwstring_data_source& input) {
  IW_TDT tdt;
  while (tdt.next(input)) {
    std::unique_ptr<IW_General_Fingerprint> fp = std::make_unique<IW_General_Fingerprint>();
    int fatal;
    if (fp->construct_from_tdt(tdt, fatal)) {
    } else if (! fatal) {
      break;
    } else {
      cerr << "Options::ReadFingerprints:cannot parse tdt " << tdt << '\n';
      return 0;
    }

    _pool << fp.release();
  }

  return _pool.size();
}

int
Options::Process(IWString_and_File_Descriptor& output) {
  for (int i = 0; i < number_fingerprints(); ++i) {
    FixedFingerprint(i, output);
  }

  for (int i = 0; i < number_sparse_fingerprints(); ++i) {
    SparseFingerprint(i, output);
  }

  return 1;
}

int
Options::FixedFingerprint(int fingerprint, IWString_and_File_Descriptor& output) {
  const uint32_t nbits = (*_pool[0])[fingerprint].nbits();
  std::unique_ptr<uint32_t[]> count = std::make_unique<uint32_t[]>(nbits);
  std::fill_n(count.get(), nbits, 0);

  if (_verbose) {
    cerr << "Options::FixedFingerprint:fingerprint contains " << nbits << " bits\n";
  }

  for (const IW_General_Fingerprint* fp : _pool) {
    FixedFingerprint(fp->item(fingerprint), count.get());
  }

  Normalise(count.get(), nbits);
  uint64_t tot = 0;
  for (uint32_t i = 0; i < nbits; ++i) {
    tot += count[i];
  }

  double float_tot = static_cast<float>(tot);

  Sparse_Fingerprint_Creator sfc;
  for (uint32_t i = 0; i < nbits; ++i) {
    if (count[i] == 0) {
      continue;
    }

    int c = static_cast<int>(static_cast<double>(count[i]) / float_tot + 0.49999);
    if (c == 0) {
      c = 1;
    }
    sfc.hit_bit(i, c);
  }

  IWString tag;
  tag << "NCMOD" << fixed_fingerprint_tag(fingerprint) << fingerprint << '<';

  sfc.write_fingerprint(tag, output);

  if (_do_fixed_fp_pairwise) {
    PairwiseFixedFingerprint(fingerprint, count.get(), nbits, output);
  }

  return 1;
}

// Return the maximum value from a k/v map.
uint32_t
MaxValue(const absl::flat_hash_map<uint32_t, uint32_t>& count) {
  uint32_t result = 0;
  for (const auto [_, c] : count) {
    if (c > result) {
      result = c;
    }
  }

  return result;
}


void
Options::HashToFingerprint(const absl::flat_hash_map<uint32_t, uint32_t>& count,
                Sparse_Fingerprint_Creator& sfc) {
  if (_log10_counts) {
    HashToFingerprintLog10(count, sfc);
  } else {
    HashToFingerprintLinear(count, sfc);
  }
}

void
Options::HashToFingerprintLinear(const absl::flat_hash_map<uint32_t, uint32_t>& count,
                        Sparse_Fingerprint_Creator& sfc) const {
  const uint32_t maxcount = MaxValue(count);

  for (const auto [bit, c] : count) {
    const int s = static_cast<int>(static_cast<float>(c) * 255.0 / static_cast<float>(maxcount) + 0.4999);
    // Make this an optional thing.
    if (s == 0) {
      continue;
    }
    sfc.hit_bit(bit, s);
  }

  if (_verbose) {
    cerr << "Linear scaling starts with " << count.size() << " bits, result has " << sfc.nbits() << '\n';
  }
}

void
Options::HashToFingerprintLog10(const absl::flat_hash_map<uint32_t, uint32_t>& count,
                        Sparse_Fingerprint_Creator& sfc) const {
  const uint32_t maxcount = MaxValue(count);
  double logmax = log10(static_cast<double>(maxcount));

  for (const auto [bit, c] : count) {
    if (c == 1) {
      continue;
    }

    double logc = log10(c);

    const int s = static_cast<int>(logc * 255.0 / logmax + 0.4999);
    // Make this an optional thing.
    if (s == 0) {
      continue;
    }
    sfc.hit_bit(bit, s);
  }
  if (_verbose) {
    cerr << "Log10 scaling starts with " << count.size() << " bits, result has " << sfc.nbits() << '\n';
  }
}

void
Options::FixedFingerprint(const IWDYFP& fp, uint32_t* count) const {
  int ndx = 0;
  int bit;
  while ((bit = fp.NextOnBit(ndx)) >= 0) {
    ++count[bit];
  }
}

void
Options::PairwiseFixedFingerprint(
        const IWDYFP& fp,
        uint32_t nbits,
        uint32_t* count,
        absl::flat_hash_map<uint32_t, uint32_t>& pairs) const {
  // first identify the bits that are set.
  std::fill_n(count, nbits, 0);
  int ndx = 0;
  int bit;

  while ((bit = fp.NextOnBit(ndx)) >= 0) {
    count[bit] = 1;
  }

  for (uint32_t i = 0; i < nbits; ++i) {
    if (count[i] == 0) {
      continue;
    }
    for (uint32_t j = i + 1; j < nbits; ++j) {
      if (count[j] == 0) {
        continue;
      }

      uint32_t b = i * nbits + j;
      auto iter = pairs.find(b);
      if (iter == pairs.end()) {
        pairs[b] = 1;
      } else {
        ++iter->second;
      }
    }
  }
}

void
Options::PairwiseFixedFingerprint(
        const IWDYFP& fp,
        uint32_t nbits,
        uint32_t* count,
        uint32_t* pairs) const {
  // first identify the bits that are set.
  std::fill_n(count, nbits, 0);
  int ndx = 0;
  int bit;

  while ((bit = fp.NextOnBit(ndx)) >= 0) {
    count[bit] = 1;
  }

  for (uint32_t i = 0; i < nbits; ++i) {
    if (count[i] == 0) {
      continue;
    }
    for (uint32_t j = i + 1; j < nbits; ++j) {
      if (count[j] == 0) {
        continue;
      }

      uint32_t ndx = i * nbits + j;
      ++pairs[ndx];
    }
  }
}
int
Options::PairwiseFixedFingerprint(int fingerprint,
        uint32_t* count,
        uint32_t nbits,
        IWString_and_File_Descriptor& output) const {
  uint32_t nb2 = nbits * nbits;
  std::unique_ptr<uint32_t[]> pairs = std::make_unique<uint32_t[]>(nb2);

  for (const IW_General_Fingerprint* fp : _pool) {
    PairwiseFixedFingerprint(fp->item(fingerprint), nbits, count, pairs.get());
  }

  uint32_t maxval = *std::max_element(pairs.get(), pairs.get() + nb2);

  double log10_maxval = std::log10(maxval);

  Sparse_Fingerprint_Creator sfc;
  uint32_t nonzero = 0;
  for (uint32_t b = 0; b < nb2; ++b) {
    if (pairs[b] == 0) {
      continue;
    }

    ++nonzero;

    int c;
    if (_log10_counts) {
      c = static_cast<int>(log10(pairs[b]) * 255.0 / log10_maxval + 0.4999);
    } else {
      c = static_cast<int>(static_cast<float>(pairs[b]) * 255.0 / static_cast<float>(maxval) + 0.4999);
    }
    // cerr << " bit " << b << " count " << pairs[b] << " count in fp " << c << '\n';
    // Make this an optional thing.
    if (c == 0) {
      continue;
    }
    sfc.hit_bit(b, c);
  }

  if (_verbose) {
    cerr << "Fixed fingerprint " << fingerprint << " contains " << nonzero <<
        " pairs, max count " << maxval << ", fingerprint contains " << sfc.nbits() << " bits\n";
  }

  IWString tag;
  tag << "NCMODPAIR" << fixed_fingerprint_tag(fingerprint) << fingerprint << "<";

  return sfc.write_fingerprint(tag, output);
}

int
Options::SparseFingerprint(int fingerprint, IWString_and_File_Descriptor& output) {
  absl::flat_hash_map<uint32_t, uint32_t> count;
  for (const IW_General_Fingerprint* fp : _pool) {
    SparseFingerprint(fp->sparse_fingerprint(fingerprint), count);
  }

  if (_verbose) {
    uint32_t maxval = MaxValue(count);
    cerr << "Sparse fingerprint " << fingerprint << " contains " << count.size() << " bits, max count " << maxval << '\n';
  }

  Sparse_Fingerprint_Creator sfc;
  HashToFingerprint(count, sfc);

  IWString tag;
  tag << "NCMODSPMOD" << fingerprint << "<";
  sfc.write_fingerprint(tag, output);

  if (_do_sparse_fp_pairwise) {
    PairwiseSparseFingerprint(fingerprint, output);
  }

  return 1;
}

int
Options::PairwiseSparseFingerprint(int fingerprint,
                IWString_and_File_Descriptor& output) {
  absl::flat_hash_map<uint32_t, uint32_t> count;
  for (const IW_General_Fingerprint* fp : _pool) {
    PairwiseSparseFingerprint(fp->sparse_fingerprint(fingerprint), count);
  }

  if (_verbose) {
    const uint32_t maxval = MaxValue(count);
    cerr << "Sparse fingerprint " << fingerprint << " pairs have " << count.size() << " bits, max count " << maxval << '\n';
  }

  Sparse_Fingerprint_Creator sfc;
  HashToFingerprint(count, sfc);

  IWString tag;
  tag << "NCMODSPR" << fingerprint << '<';
  return sfc.write_fingerprint(tag, output);
}

void
Options::PairwiseSparseFingerprint(const Sparse_Fingerprint& fp,
                absl::flat_hash_map<uint32_t, uint32_t>& count) {
  // We don't know the 'width' of a sparse fingerprint, so we use an arbitrary
  // number to form the pairwise bits. The calculations involved almost
  // certainly overflow, which is fine.
  static constexpr uint32_t kMagic = 123456;

  int ndx1 = 0;
  uint32_t b1;
  int count1;
  while (fp.next_bit_set(ndx1, b1, count1)) {
    int ndx2 = ndx1 + 1;
    uint32_t b2;
    int count2;
    while (fp.next_bit_set(ndx2, b2, count2)) {
      if (b1 == b2) {
        continue;
      }
      uint32_t key;
      if (b1 < b2) {
        key = b1 * kMagic + b2;
      } else {
        key = b2 * kMagic + b1;
      }
      auto iter = count.find(key);
      if (iter == count.end()) {
        count[key] = 1;
      } else {
        ++iter->second;
      }
    }
  }
}

void
Options::SparseFingerprint(const Sparse_Fingerprint& fp,
                absl::flat_hash_map<uint32_t, uint32_t> & count) const {
  
  int ndx = 0;
  uint32_t bit;
  int notused;
  // cerr << "Sparse fp has " << fp.nset() << " bits\n";
  while (fp.next_bit_set(ndx, bit, notused)) {
    auto iter = count.find(bit);
    if (iter == count.end()) {
      count[bit] = 1;
    } else {
      ++iter->second;
    }
  }
}

void
Options::Normalise(uint32_t* count, uint32_t nbits) const {
  uint32_t maxval = *std::max_element(count, count + nbits);

  for (uint32_t i = 0; i < nbits; ++i) {
    if (count[i] == 0) {
      continue;
    }
    const int c = static_cast<int>(static_cast<float>(count[i]) / static_cast<float>(maxval) + 0.4999);
    if (c == 0) {
      continue;
    }
    count[i] = c;
  }
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vF:P:p:l");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (need_to_call_initialise_fingerprints(cl)) {
    if (! initialise_fingerprints(cl, verbose)) {
      cerr << "Cannot initialise GFP options\n";
      Usage(1);
    }
  }

  if (cl.empty()) {
    cerr << "Insuffient arguments\n";
    Usage(1);
  }

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise\n";
    return 1;
  }

  IWString_and_File_Descriptor output(1);

  if (! options.Process(output)) {
    cerr << "Modal fingerprint failed\n";
    return 1;
  }

  return 0;
}

}  // namespace gfp_model_fingerprint

int
main(int argc, char ** argv)
{
  return gfp_model_fingerprint::Main(argc, argv);
}
