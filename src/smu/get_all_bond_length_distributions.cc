// Generate bond length distributions for all bonded types.

#include <stdlib.h>

#include <algorithm>
#include <iostream>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"

namespace get_bond_length_distribution {

using std::cerr;

void
Usage(int rc) {
  ::exit(rc);
}

// Copied from get_bond_length_distribution.cc
struct BucketisedDistances
{
  // From bucketised distance to count.
  std::map<int, int> bond_length_data;

  void Extra(float distance);

  int Write(const char sep, IWString_and_File_Descriptor& output) const;
  // int WriteProto(Smu::BondLengthDistribution& proto, IWString& fname) const;
};

void
BucketisedDistances::Extra(float distance)
{
  const int b = static_cast<int>(distance * 1000.0f);
  bond_length_data[b]++;
}

int
BucketisedDistances::Write(const char sep,
                           IWString_and_File_Descriptor& output) const
{
  Accumulator_Int<uint64_t> acc;
  for (auto [bucket, count] : bond_length_data) {
    output << bucket / 1000.0f << sep << count << '\n';
    acc.extra(bucket, count);
  }

  cerr << acc.n() << " values btw " << acc.minval() << " and " << acc.maxval() << " mean " << acc.average() << '\n';

  return 1;
}

#ifdef BUCKETS_HAVE_PROTO
int
BucketisedDistances::WriteProto(Smu::BondLengthDistribution& proto,
                                IWString& fname) const
{
  Accumulator_Int<int> acc;
  int max_count = 0;
  int ndx_max_count = -1;
  for (auto [bucket, count] : bond_length_data) {
    acc.extra(bucket, count);
    if (count > max_count) {
      max_count = count;
      ndx_max_count = bucket;
    }
  }

  proto.set_min(acc.minval() / 1000.0f);
  proto.set_max(acc.maxval() / 1000.0f);
  proto.set_mean(ndx_max_count / 1000.0f);

  return iwmisc::WriteTextProto(proto, fname);
}
#endif


class BondLengthDistribution {
  private:
    // Map from a hash of two atomic numbers and a bond type to a BucketisedDistances.
    std::unordered_map<int, BucketisedDistances> _bld;

  // private functions
    // A hash from two atomic numbers and a bond type to a unique integer.
    int Hash(int atnum1, int btype, int atnum2) const;

  public:
    int AccumulateDistances(Molecule& m);

    int WriteDistances(const IWString& stem);
};

// The hash function works because all atomic numbers are < 10. That way
// we can generate numbers that can be converted back and forward.
int
BondLengthDistribution::Hash(int atnum1, int btype, int atnum2) const {
  if (atnum1 > atnum2) {
    std::swap(atnum1, atnum2);
  }

  return atnum1 * 100 + btype * 10 + atnum2;
}

std::tuple<int, int, int>
HashToComponents(int h) {
  int atnum1 = h / 100;
  h %= 100;
  int btype = h / 10;
  h %= 10;
  return std::tuple<int, int, int>(atnum1, btype, h);
}

int
BtypeToNumber(const Bond* b) {
  if (b->is_single_bond()) {
    return 1;
  }
  if (b->is_double_bond()) {
    return 2;
  }
  if (b->is_triple_bond()) {
    return 3;
  }

  cerr << "Unrecognised bond type\n";  // cannot happen
  return 0;
}

int
BondLengthDistribution::AccumulateDistances(Molecule& m) {
  const int nedges = m.nedges();
  for (int i = 0; i < nedges; ++i) {
    const Bond * b = m.bondi(i);
    const Atom * a1 = m.atomi(b->a1());
    const Atom * a2 = m.atomi(b->a2());
    const atomic_number_t atnum1 = a1->atomic_number();
    const atomic_number_t atnum2 = a2->atomic_number();
    const float d = m.distance_between_atoms(b->a1(), b->a2());

    int h = Hash(atnum1, BtypeToNumber(b), atnum2);

    auto iter = _bld.find(h);
    if (iter != _bld.end()) {
      iter->second.Extra(d);
    } else {
      auto [iter, _] = _bld.emplace(h, BucketisedDistances());
      iter->second.Extra(d);
    }
  }

  return 1;
}

int
BondLengthDistribution::WriteDistances(const IWString& stem) {
  for (const auto& [h, bld] : _bld) {
    auto [atnum1, btype, atnum2] = HashToComponents(h);
    IWString fname;
    fname << stem << '.' << atnum1 << '.' << btype << '.' << atnum2;

    IWString_and_File_Descriptor output;
    if (! output.open(fname.null_terminated_chars())) {
      cerr << "BondLengthDistribution::WriteDistances:cannot open '" << fname << "'\n";
      return 0;
    }
    bld.Write(' ', output);
  }

  return 1;
}

int
GetBondLengthDistribution(BondLengthDistribution& bond_length_distribution,
            data_source_and_type<Molecule>& input) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! bond_length_distribution.AccumulateDistances(*m)) {
      return 0;
    }
  }

  return 1;
}

int
GetBondLengthDistribution(const char * fname, 
                          FileType input_type,
                          BondLengthDistribution& bond_length_distribution) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return GetBondLengthDistribution(bond_length_distribution, input);
}

int
GetBondLengthDistribution(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:E:S:i:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (! process_standard_aromaticity_options(cl, verbose, 'A')) {
    cerr << "Cannot process aromaticity options\n";
    return 1;
  }

  if (! process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process standard elements options (-E)\n";
    return 1;
  }

  if (! cl.option_present('S')) {
    cerr << "Must specify output file name stem via the -S option\n";
    Usage(1);
  }

  const IWString output_stem = cl.string_value('S');

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (! cl.option_present('i')) {
    if (1 == cl.number_elements() && 0 == strcmp("-", cl[0])) {
      input_type = FILE_TYPE_SMI;
    } else if (! all_files_recognised_by_suffix(cl)) {
      cerr << "Cannot auto-detect file types\n";
      return 1;
    }
  } else if (! process_input_type(cl, input_type)) {
    cerr << "cannot discern input type\n";
    return 1;
  }

  BondLengthDistribution bld;

  for (const char * fname : cl) {
    if (! GetBondLengthDistribution(fname, input_type, bld)) {
      cerr << "Fatal error processing '" << fname << "\n";
      return 1;
    }
  }

  bld.WriteDistances(output_stem);

  return 0;
}

}  // namespace get_bond_length_distribution

int
main(int argc, char** argv) {
  const int rc = get_bond_length_distribution::GetBondLengthDistribution(argc, argv);

  return rc;
}
