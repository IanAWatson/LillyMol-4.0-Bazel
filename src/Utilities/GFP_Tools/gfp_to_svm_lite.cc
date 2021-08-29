// Converts GFP files to svm lite.
// Primary function is to assign a unique feature number to each
// bit in the gfp file(s).

#include <algorithm>
#include <memory>
#include <unordered_map>

#include "google/protobuf/text_format.h"
#include "google/protobuf/io/zero_copy_stream.h"
#include "google/protobuf/io/zero_copy_stream_impl.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#define ACTIVITY_DATA_IMPLEMENATION_H
#include "Foundational/iwmisc/activity_data_from_file.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "Utilities/GFP_Tools/gfp.h"
#include "Utilities/GFP_Tools/gfp_to_svm_lite.pb.h"

namespace gfp_to_svm_lite {

using std::cerr;

// The type of features that are used as feature numbers in the svm_lite input.
using feature_type_t = uint32_t;
// The type of bit numbers that arise from gfp files.
using gfp_bit_type_t = uint32_t;

int verbose = 0;

IWString smiles_tag = "$SMI<";
IWString identifier_tag = "PCN<";

// The activity can be part of the name.
int activity_column = -1;

IWDigits bit_number;
IWDigits count;

char sep = ' ';

int support = 0;

int flatten_sparse_fingerprints = 0;

// We operate in one of two modes. Creating or using
// a proto file.
bool create_bit_map = false;
// If we are using an existing proto map file, it can only
// be processed after we have read the first gfp item.
bool use_bit_map = false;
IWString input_proto_fname;

void
Usage(int rc) {
  exit(rc);
}

// Class to enable conversion from GFP bit numbers to feature
// numbers used by svm-lite.
// It operates in two phases.
// 1. Discern the bits in the input and a conversion to feature number.
// 2. Use that conversion table to do the conversion.
// Step 1 can be done either by looking at the input, the ProfileBits functions
// below (create_bit_map), or by reading a precomputed conversion table (use_bit_map).
class GfpToFeature {
  private:
    // As new GFP bits are encountered, we assign a globally unique feature
    // number;
    feature_type_t _next_feature;
    // 
    resizable_array<std::unordered_map<gfp_bit_type_t, feature_type_t>> _dense_gfp_to_bit;
    resizable_array<std::unordered_map<gfp_bit_type_t, feature_type_t>> _sparse_gfp_to_bit;

    // The names of the fingerprints - used when writing the proto.
    resizable_array<IWString> _dense_gfp_name;
    resizable_array<IWString> _sparse_gfp_name;

    extending_resizable_array<int> _molecules_containing_feature;

    // If a support requirement is imposed, we store it and write with
    // the proto.
    int _support;

    // private functions

    // Generic function that tries to lookup `gfp_bit` in `xref`.
    // If successful, it returns xref[gfp_bit]. If not, it increments
    // _next_feature, sets xref[gfp_bit] = _next_feature.
    feature_type_t _FeatureForBit(std::unordered_map<gfp_bit_type_t, feature_type_t>& xref, gfp_bit_type_t gfp_bit);
    std::optional<feature_type_t> _FeatureForBit(const std::unordered_map<gfp_bit_type_t, feature_type_t>& xref, gfp_bit_type_t gfp_bit) const;

    GfpBitToFeatureMap::BitXref MakeProto(const std::unordered_map<gfp_bit_type_t, feature_type_t>& xref) const;

    void MakeXref(const GfpBitToFeatureMap::BitXref& proto, std::unordered_map<gfp_bit_type_t, feature_type_t>& xref);

  public:
    GfpToFeature();

    // Uses the globally available counts of dense and sparse fingerprints
    // to initialize.
    int Initialize();

    // Building from protos.
    int FromProto(const char * fname);
    int FromProto(const GfpBitToFeatureMap::GfpBitXref& proto);

    int DebugPrint(std::ostream& output) const;

    // Given a feature coming from either dense fingerprint `ndx` or
    // sparse fingerprint `ndx`, return the corresponding bit number.
    // Updates internal hashes
    feature_type_t FetchOrAssignFeatureForDenseBit(int ndx, gfp_bit_type_t gfp_bit);
    feature_type_t FetchOrAssignFeatureForSparseBit(int ndx, gfp_bit_type_t gfp_bit);

    // Const version that only returns a valid feature number if
    // `gfp_bit` is in the hashes.
    std::optional<feature_type_t> FeatureForDenseBit(int ndx, gfp_bit_type_t gfp_bit) const;
    std::optional<feature_type_t> FeatureForSparseBit(int ndx, gfp_bit_type_t gfp_bit) const;

    // Remove all features that occur in fewer than `support` molecules.
    // Returns the number of features discarded.
    int ImposeSupport(int support);

    // Convert to proto form.
    GfpBitToFeatureMap::GfpBitXref ToProto() const;

    // Write the bit cross reference data in proto form.
    int WriteProto(const char * fname) const;
};

GfpToFeature::GfpToFeature() {
  _next_feature = 0;
  _support = 0;
}

int
GfpToFeature::Initialize() {
  const int nfixed = number_fingerprints();
  const int nsparse = number_sparse_fingerprints();

  _dense_gfp_to_bit.extend(nfixed, std::unordered_map<gfp_bit_type_t, feature_type_t>());
  _sparse_gfp_to_bit.extend(nsparse, std::unordered_map<gfp_bit_type_t, feature_type_t>());
  _dense_gfp_name.extend(nfixed, IWString());
  _sparse_gfp_name.extend(nsparse, IWString());
  for (int i = 0; i < nfixed; ++i) {
    _dense_gfp_name[i] = fixed_fingerprint_tag(i);
  }
  for (int i = 0; i < nsparse; ++i) {
    _sparse_gfp_name[i] = sparse_fingerprint_tag(i);
  }
  if (verbose) {
    cerr << "Initialised for " << nfixed << " fixed and " << nsparse << " sparse fingerprints\n";
  }
  return 1;
}

int
GfpToFeature::DebugPrint(std::ostream& output) const {
  output << "GfpToFeature: with " << _dense_gfp_to_bit.number_elements() << " fixed and " << _sparse_gfp_to_bit.number_elements() << " fingerprints\n";
  output << "_molecules_containing_feature contains " << _molecules_containing_feature.number_elements() << " items\n";
  for (int i = 0; i < _molecules_containing_feature.number_elements(); ++i) {
    if (_molecules_containing_feature[i] > 0) {
      cerr << "  count " << i << '\n';
    }
  }
  return output.good();
}

int
GfpToFeature::FromProto(const char * fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "GfpToFeature::FromProto:cannot open '" << fname << "'\n";
    return 0;
  }
 
  google::protobuf::io::FileInputStream file_input_stream(input.fd());

  GfpBitToFeatureMap::GfpBitXref proto;

  if (! google::protobuf::TextFormat::Parse(&file_input_stream, &proto)) {
    cerr << "GfpToFeature::FromProto:cannot parse '" << fname << "'\n";
    return 0;
  }

  return FromProto(proto);
}

// Given a mapping from gfp bit number in `proto`, fill in the
// corresponding unordered_map `xref`. It is really just a copy.
void
GfpToFeature::MakeXref(const GfpBitToFeatureMap::BitXref& proto,
                               std::unordered_map<gfp_bit_type_t, feature_type_t>& xref) {
  for (const auto& [key, value] : proto.gfp_to_svm()) {
    xref.emplace(key, value);
    _molecules_containing_feature[value] = 1;
  }
}

int
GfpToFeature::FromProto(const GfpBitToFeatureMap::GfpBitXref& proto) {
  Initialize();
  const int nfixed = number_fingerprints();
  const int nsparse = number_sparse_fingerprints();
  for (int i = 0; i < nfixed; ++i) {
    const std::string tag(_dense_gfp_name[i].data(), _dense_gfp_name[i].length());
    auto iter = proto.xref().find(tag);
    if (iter == proto.xref().end()) {
      cerr << "GfpToFeature::FromProto:no data for " << _dense_gfp_name[i] << "'\n";
      return 0;
    }
    MakeXref(iter->second, _dense_gfp_to_bit[i]);
  }

  for (int i = 0; i < nsparse; ++i) {
    const std::string tag(_sparse_gfp_name[i].data(), _sparse_gfp_name[i].length());
    auto iter = proto.xref().find(tag);
    if (iter == proto.xref().end()) {
      cerr << "GfpToFeature::FromProto:no data for " << _sparse_gfp_name[i] << "'\n";
      return 0;
    }
    MakeXref(iter->second, _sparse_gfp_to_bit[i]);
  }

  return 1;
}

feature_type_t
GfpToFeature::_FeatureForBit(std::unordered_map<gfp_bit_type_t, feature_type_t>& xref,
                               gfp_bit_type_t gfp_bit) {
  const auto iter = xref.find(gfp_bit);
  if (iter != xref.end()) {
    _molecules_containing_feature[iter->second]++;
    return iter->second;
  }

  _next_feature++;
  xref.emplace(gfp_bit, _next_feature);
  _molecules_containing_feature[_next_feature]++;
  return _next_feature;
}

std::optional<feature_type_t>
GfpToFeature::_FeatureForBit(const std::unordered_map<gfp_bit_type_t, feature_type_t>& xref,
                               gfp_bit_type_t gfp_bit) const {
  const auto iter = xref.find(gfp_bit);
  if (iter == xref.end()) {
    return std::nullopt;
  }

  if (_molecules_containing_feature[iter->second] == 0) {
    return std::nullopt;
  }

  return iter->second;
}

feature_type_t
GfpToFeature::FetchOrAssignFeatureForDenseBit(int ndx, gfp_bit_type_t gfp_bit) {
  return _FeatureForBit(_dense_gfp_to_bit[ndx], gfp_bit);
}

feature_type_t
GfpToFeature::FetchOrAssignFeatureForSparseBit(int ndx, gfp_bit_type_t feature) {
  return _FeatureForBit(_sparse_gfp_to_bit[ndx], feature);
}

std::optional<feature_type_t>
GfpToFeature::FeatureForDenseBit(int ndx, gfp_bit_type_t gfp_bit) const {
  return _FeatureForBit(_dense_gfp_to_bit[ndx], gfp_bit);
}

std::optional<feature_type_t>
GfpToFeature::FeatureForSparseBit(int ndx, gfp_bit_type_t feature) const {
  return _FeatureForBit(_sparse_gfp_to_bit[ndx], feature);
}

int
GfpToFeature::ImposeSupport(int support) {
  _support = support;
  int rc = 0;
  for (int f = 0; f < _molecules_containing_feature.number_elements(); ++f) {
    const int c = _molecules_containing_feature[f];
    if (c >= support) {
      continue;
    }
    _molecules_containing_feature[f] = 0;
    rc++;
  }
  return rc;
}

// Convert an unordered_map to a BitXref proto.
GfpBitToFeatureMap::BitXref
GfpToFeature::MakeProto(const std::unordered_map<gfp_bit_type_t, feature_type_t>& xref) const {
  GfpBitToFeatureMap::BitXref result;
  for (const auto& [key, value] : xref) {
    if (_molecules_containing_feature[value] == 0) {
      continue;
    }
    google::protobuf::MapPair<uint64_t, uint64_t> to_insert(key, value);
    result.mutable_gfp_to_svm()->insert(to_insert);
  }
  return result;
}

GfpBitToFeatureMap::GfpBitXref
GfpToFeature::ToProto() const {
  GfpBitToFeatureMap::GfpBitXref result;
  const int nfixed = number_fingerprints();
  const int nsparse = number_sparse_fingerprints();
  for (int i = 0; i < nfixed; ++i) {
    GfpBitToFeatureMap::BitXref xref = MakeProto(_dense_gfp_to_bit[i]);

    const std::string tag_name(_dense_gfp_name[i].data(), _dense_gfp_name[i].length());
    google::protobuf::MapPair<std::string, GfpBitToFeatureMap::BitXref> to_insert(tag_name, xref);

    result.mutable_xref()->insert(to_insert);
  }

  for (int i = 0; i < nsparse; ++i) {
    GfpBitToFeatureMap::BitXref xref = MakeProto(_sparse_gfp_to_bit[i]);

    const std::string tag_name(_sparse_gfp_name[i].data(), _sparse_gfp_name[i].length());
    google::protobuf::MapPair<std::string, GfpBitToFeatureMap::BitXref> to_insert(tag_name, xref);

    result.mutable_xref()->insert(to_insert);
  }

  if (_support > 0) {
    result.set_support(_support);
  }

  return result;
}

int
GfpToFeature::WriteProto(const char * fname) const {
  const GfpBitToFeatureMap::GfpBitXref proto = ToProto();
  IWString_and_File_Descriptor output;
  if (! output.open(fname)) {
    cerr << "GfpToFeature::WriteProto:cannot open '" << fname << "'\n";
    return 0;
  }

  using google::protobuf::io::ZeroCopyOutputStream;
  using google::protobuf::io::FileOutputStream;
  std::unique_ptr<ZeroCopyOutputStream> zero_copy_output(new FileOutputStream(output.fd()));
  if (! google::protobuf::TextFormat::Print(proto, zero_copy_output.get())) {
    cerr << "GfpToFeature::WriteProto:cannot write " << fname << "\n";
    return 0;
  }

  return 1;
}

struct FeatureCount {
  feature_type_t feature;
  int count;
};

void
CollectFeatures(const IWDYFP& fp,
                int fpnum,
                const GfpToFeature& bit_xref,
                FeatureCount* feature_count,
                int & ndx) {
  int b = 0;
  while (fp.next_on_bit(b)) {
    std::optional<feature_type_t> feature = bit_xref.FeatureForDenseBit(fpnum, b);
    if (! feature) {
      continue;
    }
    feature_count[ndx].feature = *feature;
    feature_count[ndx].count = 1;
    ndx++;
  }
}

void
CollectFeatures(const IWDYFP& fp,
                int ndx,
                const GfpToFeature& bit_xref,
                extending_resizable_array<int>& feature_count) {
  int b = 0;
  while (fp.next_on_bit(b)) {
    std::optional<feature_type_t> feature = bit_xref.FeatureForDenseBit(ndx, b);
    if (! feature) {
      continue;
    }
    feature_count[*feature] = 1;
  }
}

void
CollectFeatures(const Sparse_Fingerprint& fp,
                int fpnum,
                const GfpToFeature& bit_xref,
                FeatureCount* feature_count,
                int& ndx) {
  int i = 0;
  gfp_bit_type_t b = 0;
  int c = 0;
  while (fp.next_bit_set(i, b, c)) {
    std::optional<feature_type_t> feature = bit_xref.FeatureForSparseBit(fpnum, b);
    if (! feature) {
      continue;
    }
    feature_count[ndx].feature = *feature;
    if (flatten_sparse_fingerprints) {
      feature_count[ndx].count = 1;
    } else {
      feature_count[ndx].count = c;
    }
    ndx++;
  }
}

void
CollectFeatures(const Sparse_Fingerprint& fp,
                int ndx,
                const GfpToFeature& bit_xref,
                extending_resizable_array<int>& feature_count) {
  int i = 0;
  gfp_bit_type_t b = 0;
  int c = 0;
  while (fp.next_bit_set(i, b, c)) {
    std::optional<feature_type_t> feature = bit_xref.FeatureForSparseBit(ndx, b);
    if (! feature) {
      continue;
    }
    if (flatten_sparse_fingerprints) {
      feature_count[*feature] = 1;
    } else {
      feature_count[*feature] += c;
    }
  }
}

template <typename Activity>
int
GfpToSvmLite(const IW_General_Fingerprint& gfp,
             const IWString& id,
             Activity activity,
             GfpToFeature& bit_xref,
             IWString_and_File_Descriptor& output) {
  static bool first_call = true;
  if (first_call) {
    first_call = false;
    if (use_bit_map) {
      if (! bit_xref.FromProto(input_proto_fname)) {
        cerr << "GfpToSvmLite:cannot read proto " << input_proto_fname << "\n";
        return 0;
      }
    }
  }

  const int nfixed = number_fingerprints();
  const int nsparse = number_sparse_fingerprints();
  int nfeatures = 0;

  for (int i = 0; i < nfixed; ++i) {
    nfeatures += gfp[i].nset();
  }
  for (int i = 0; i < nsparse; ++i) {
    nfeatures += gfp.sparse_fingerprint(i).nset();
  }

  std::unique_ptr<FeatureCount[]> features(new FeatureCount[nfeatures]);
  int ndx = 0;
  for (int i = 0; i < nfixed; ++i) {
    CollectFeatures(gfp[i], i, bit_xref, features.get(), ndx);
  }
  for (int i = 0; i < nsparse; ++i) {
    CollectFeatures(gfp.sparse_fingerprint(i), i, bit_xref, features.get(), ndx);
  }
  std::sort(features.get(), features.get() + ndx, [](const FeatureCount& fc1, const FeatureCount& fc2) {
    return fc1.feature < fc2.feature;
  });

  output << activity;
  for (int i = 0; i < ndx; ++i) {
    bit_number.append_number(output, features[i].feature);
    count.append_number(output, features[i].count);
  }

  output << " # " << id << '\n';
  output.write_if_buffer_holds_more_than(8192);

  return 1;
}

void
ProfileBits(const IWDYFP& fp, int ndx, GfpToFeature& bit_xref) {
  int bit = 0;
  while (fp.next_on_bit(bit)) {
    bit_xref.FeatureForDenseBit(ndx, bit);
  }
}

void
ProfileBits(const Sparse_Fingerprint& fp, int ndx, GfpToFeature& bit_xref) {
  int i = 0;
  gfp_bit_type_t b = 0;
  int c = 0;
  while (fp.next_bit_set(i, b, c)) {
    bit_xref.FetchOrAssignFeatureForSparseBit(ndx, b);
  }
}

int
ProfileBits(const IW_General_Fingerprint& gfp,
             GfpToFeature& bit_xref) {
  const int nfixed = number_fingerprints();
  for (int i = 0; i < nfixed; ++i) {
    ProfileBits(gfp[i], i, bit_xref);
  }

  const int nsparse = number_sparse_fingerprints();
  for (int i = 0; i < nsparse; ++i) {
    ProfileBits(gfp.sparse_fingerprint(i), i, bit_xref);
  }

  return 1;
}

int
EchoItem(const IW_TDT& tdt,
         const IWString& tag,
         IWString_and_File_Descriptor& output) {
  IWString s;
  if (!tdt.dataitem_value(tag, s, 0)) {
    cerr << "Cannot extract " << tag << " from " << tdt << '\n';
    return 0;
  }
  output << tag << s << ">\n";
  return 1;
}

int
EchoSmilesID(const IW_TDT& tdt,
             IWString_and_File_Descriptor& output)  {
  if (! EchoItem(tdt, smiles_tag, output) ||
      ! EchoItem(tdt, identifier_tag, output)) {
    cerr << "Cannot echo smiles/id\n";
    return 0;
  }

  return 1;
}

template <typename Activity>
std::optional<Activity>
ActivityFromColumn(const IWString& id,
                   int activity_column) {
  int i = 0;
  const_IWSubstring token;
  for (int col = 0; id.nextword(token, i); ++col) {
    if (col != activity_column) {
      continue;
    }
    Activity result;
    if (! token.numeric_value(result)) {
      cerr << "ActivityFromColumn:invalid numeric " << id << "'\n";
      return std::nullopt;
    }
    return result;
  }
  cerr << "ActivityFromColumn:did not encounter column " << activity_column << '\n';
  return std::nullopt;
}

// Get the activity for `id` from either `activity` or from
// column `activity_column` of `id`;
template <typename Activity>
std::optional<Activity>
GetActivity(const Activity_Data_From_File<Activity>& activity,
            const IWString& id,
            int activity_column) {
  if (activity_column > 0)
    return ActivityFromColumn<Activity>(id, activity_column);

  const auto iter = activity.find(id);
  if (iter != activity.end()) {
    return iter->second;
  }

  if (id.nwords() == 1) {
    cerr << "No activity for '" << id << "'\n";
    return std::nullopt;
  }

  IWString tmp(id);
  tmp.truncate_at_first(' ');
  return GetActivity(activity, tmp, activity_column);
}

template <typename Activity>
int
GfpToSvmLite(IW_TDT& tdt,
             const Activity_Data_From_File<Activity>& activity,
             GfpToFeature& bit_xref,
             IWString_and_File_Descriptor& output) {
  IW_General_Fingerprint gfp;
  int fatal = 0;
  if (! gfp.construct_from_tdt(tdt, fatal)) {
    if (fatal) {
      return 0;
    }
    return 1;
  }

  IWString id;
  if (! tdt.dataitem_value(identifier_tag, id, 0) || id.empty()) {
    cerr << "Cannot extract" << identifier_tag << " from " << tdt << '\n';
    return 0;
  }

  std::optional<Activity> act = GetActivity(activity, id, activity_column);
  if (! act) {
    return 0;
  }

  return GfpToSvmLite(gfp, id, *act, bit_xref, output);
}

int
ProfileBits(IW_TDT& tdt,
            GfpToFeature& bit_xref) {
  IW_General_Fingerprint gfp;
  int fatal = 0;
  if (! gfp.construct_from_tdt(tdt, fatal)) {
    if (fatal) {
      return 0;
    }
    return 1;
  }

  static bool first_call = true;
  if (first_call) {
    first_call = false;
    bit_xref.Initialize();
  }

  return ProfileBits(gfp, bit_xref);
}

template <typename Activity>
int
GfpToSvmLite(iwstring_data_source& input,
             const Activity_Data_From_File<Activity>& activity,
             GfpToFeature& bit_xref,
             IWString_and_File_Descriptor& output) {
  IW_TDT tdt;
  while (tdt.next(input)) {
    if (! GfpToSvmLite(tdt, activity, bit_xref, output)) {
      cerr << "Cannot process " << tdt << '\n';
      return 0;
    }
  }

  return 1;
}

int
ProfileBits(iwstring_data_source& input,
            GfpToFeature& bit_xref) {
  IW_TDT tdt;
  while (tdt.next(input)) {
    if (! ProfileBits(tdt, bit_xref)) {
      cerr << "Cannot process " << tdt << '\n';
      return 0;
    }
  }

  return 1;
}

template <typename Activity>
int
GfpToSvmLite(const char * fname,
             const Activity_Data_From_File<Activity>& activity,
             GfpToFeature& bit_xref,
             IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Cannot open " << fname << '\n';
    return 0;
  }

  return GfpToSvmLite(input, activity, bit_xref, output);
}

int
ProfileBits(const char * fname,
             GfpToFeature& bit_xref) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Cannot open " << fname << '\n';
    return 0;
  }

  return ProfileBits(input, bit_xref);
}

int
GfpToSvmLite(int argc, char** argv) {
  Command_Line cl(argc, argv, "vF:P:C:U:A:p:fc:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  verbose = cl.option_count('v');

  if (need_to_call_initialise_fingerprints(cl))
  {
    if (! initialise_fingerprints(cl, verbose))
    {
      cerr << "Cannot initialise fingerprints\n";
      return 1;
    }
  }

  Activity_Data_From_File<float> activity;
  if (cl.option_present('c') && cl.option_present('A')) {
    cerr << "Can use just one of the -c or -A options\n";
    Usage(1);
  }

  if (cl.option_present('c')) {
    if (! cl.value('c', activity_column) || activity_column < 1) {
      cerr << "Invalid activity column value (-c)\n";
      return 1;
    }
    if (verbose)
      cerr << "Activity values extracted from column " << activity_column << '\n';
    activity_column--;
  } else if (cl.option_present('A')) {
    if (! activity.construct_from_command_line(cl, 'A', verbose)) {
      cerr << "Cannot read activity data (-A)\n";
      return 1;
    }
  } else {
    cerr << "MUst specify a source of activity via either -A or -c options\n";
    Usage(1);
  }

  if (cl.option_present('p')) {
    if (! cl.value('p', support) || support < 1) {
      cerr << "Invalid support value (-p)\n";
      Usage(1);
    }

    if (verbose)
      cerr << "WIll only output features found in " << support << " molecules\n";
  }

  if (cl.option_present('f')) {
    flatten_sparse_fingerprints = 1;
    if (verbose)
      cerr << "Sparse fingerprints flattened to binary\n";
  }

  if (cl.empty()) {
    cerr << "Must specify input file(s)\n";
    Usage(1);
  }

  GfpToFeature bit_xref;

  // Branching depending on whether we are creating or using a bit
  // to feature cross reference.
  IWString output_proto_fname;
  if (cl.option_present('C')) {
    cl.value('C', output_proto_fname);
    create_bit_map = true;
    for (const char * fname : cl) {
      cerr << "Profiling " << fname << '\n';
      if (! ProfileBits(fname, bit_xref)) {
        cerr << "Cannot profile bits in " << fname << '\n';
        return 1;
      }
    }
  } else if (cl.option_present('U')) {
    cl.value('U', input_proto_fname);
    use_bit_map = true;
  } else {
    cerr << "Must specify whether creating (-C) or using (-U) a cross reference\n";
    Usage(1);
  }

  if (use_bit_map && support > 0) {
    cerr << "Specifying support and using an existing profile not supported\n";
    Usage(1);
  }
  if (create_bit_map && support > 0) {
    int features_removed = bit_xref.ImposeSupport(support);
    if (verbose)
      cerr << "Support requirement " << support << " suppressed " << features_removed << " features\n";
  }

  IWString_and_File_Descriptor output(1);

  bit_number.set_include_leading_space(1);
  bit_number.initialise(10000);
  count.set_leading_string(':');
  count.initialise(256);

  for (const char * fname : cl) {
    cerr << "Begin processing " << fname << '\n';
    if (! GfpToSvmLite(fname, activity, bit_xref, output)) {
      cerr << "Fatal error processing " << fname << '\n';
      return 1;
    }
  }

  if (create_bit_map)
    bit_xref.WriteProto(output_proto_fname);

  return 0;
}

}  // namespace gfp_to_svm_lite

int
main(int argc, char ** argv) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  return gfp_to_svm_lite::GfpToSvmLite(argc, argv);
}
