// filter for gfp fingerprints

#include <stdlib.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>

#include "re2/re2.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "Utilities/GFP_Tools/gfp.h"
#include "Utilities/GFP_Tools/gfp_standard.h"

#include "smiles_id_dist.h"

namespace gfp_distance_filter_standard {

using std::cerr;

constexpr const char* kSmilesTag = "$SMI<";
constexpr const char * kIdentifierTag = "PCN<";
constexpr const char * kDistanceTag = "DIST<";

void Usage(int rc) {
  cerr << "Neighbour finding and distance filtering\n";
  cerr << "Uses 'needles' and 'haystack' concept. Needles specified via the -p options, haystack as input\n";
  cerr << "Calculation type governed by the -C option\n";
  cerr << " -C 1             identify needles that are within -t distance of any haystack member\n";
  cerr << " -t <dist>        distance threshold - needles must be within <dist> of a haystack member\n";
  cerr << " -C 2             for each haystack member, identify the nearest needle\n";
  cerr << " -p <fname>       file containing fingerprints of the needles\n";
  cerr << " -a               find the closest match - not just the first within threshold\n";
  cerr << " -G <tag>         specify guard fingerprint tag\n";
  cerr << " -G thr=<dist>    the distance associated with the guard fingerprint\n";
  cerr << " -w <diff>        only compare fp's where the atom count differs by <diff> or less\n";
  cerr << " -n               create output for nnprocess rather than tabular\n";
  cerr << " -r <n>           report progress every <n> haystack members processed\n";
  cerr << " -X <fname>       write needles not matched to <fname>\n";
  cerr << " -v               verbose output\n";
  ::exit(rc);
}

class Options {
  private:
    int _verbose;

    // We can do several different kinds of calculations.
    enum class CalculationType {
      // Identify needles that are within a given distance of any member of the haystack.
      // Uses _threshold as the distance.
      kIdentifyNeedlesNearHaystack,

      // For each member of the haystack, identify the nearest neighbour among _needles
      // Things beyond _threshold are ignored.
      kFindNearestNeedleNeighbour,
    };

    CalculationType _calculation_type;

    // The -t option.
    float _threshold;

    GFP_Standard* _needles;
    int _number_needles;

    int _atom_count_window;
    int _calculations_bypassed_by_atom_count_window;

    // If we have a guard fingerprint, the tag.
    IWString _guard_fp_tag;
    // And a FixedBitVector for each fingerprint in _needles.
    // I thought of extending the GFP_Standard class, but this is likely
    // more efficient.
    IWDYFP* _guard_fp;

    // And the distance associated with the guard fp
    float _guard_fp_distance;

    // We can assess the impact of the guard fingerprint.
    int _calculations_attempted;
    int _calculations_bypassed_by_guard;

    // GFP_Standard do not have a name field. Create one.
    Smiles_ID_Dist* _smiles_id_dist;

    // Once a needle has been identified as within range, we can optionally
    // stop searching it.
    int * _already_written;
    // If we are stopping searching, we keep track of the number of active needles
    // and can stop when that reaches zero.
    int _active_needles;

    // this does not work. When looking for needles that are within range of a haystack
    // item, it is written as soon as found with that initial distance. Not sure this
    // is needed. It will result in no output till all the haystack has been searched.
    int _find_all_matches;

    Report_Progress _report_progress;

    uint32_t _tdts_read;

    // Of the calculations that return a value in range
    Accumulator<double> _distance;
    int _exact_matches;

    // For items that are not matched by the haystack, we can write them.
    IWString_and_File_Descriptor _stream_for_not_matched;

    int _non_matches;

    int _tabular_output;

    uint32_t _for_testing_stop;

    char _output_separator;

  // private functions
    int ReadNeedles(IWString& fname);
    int ReadNeedles(iwstring_data_source& input);
    int HandleWriteNonMatch(const IWString& id, const IW_TDT& tdt);
    int CanBeCompared(const GFP_Standard& fp1, const GFP_Standard& fp2);
    int DoOutput(const IWString& smiles1, const IWString& id1,
         const IWString& smiles2, const IWString& id2,
         float distance,
         IWString_and_File_Descriptor& output);

    int WithinRange(IWDYFP& g1, IWDYFP& g2);

    void UpdateDistanceStats(float distance);
    int IdentifyNeedlesNearHaystack(GFP_Standard& gfp,
                 IWDYFP& guard_fp,
                 const IWString& id,
                 const IW_TDT& tdt,
                 IWString_and_File_Descriptor& output);
    int FindNearestNeedleNeighbour(GFP_Standard& gfp,
                 IWDYFP& guard_fp,
                 const IWString& id,
                 const IW_TDT& tdt,
                 IWString_and_File_Descriptor& output);
    int ReportProgress(std::ostream& output) const;

  public:
    Options();
    ~Options();

    int Initialise(Command_Line& cl);

    int Process(const char * fname, IWString_and_File_Descriptor& output);
    int Process(iwstring_data_source& input, IWString_and_File_Descriptor& output);
    int Process(const IW_General_Fingerprint& gfp,
                const IW_TDT& tdt,
                IWString_and_File_Descriptor& output);
    int Process(GFP_Standard& gfp,
                IWDYFP& guard_fp,
                const IWString& id,
                const IW_TDT& tdt,
                IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;

    int WriteNonMatches();
};

Options::Options() {
  _threshold = 0.0f;

  _verbose = 0;
  _needles = nullptr;
  _guard_fp = nullptr;
  _already_written = nullptr;
  _active_needles = 0;
  _guard_fp_distance = 0.0f;
  _calculations_attempted = 0;
  _calculations_bypassed_by_guard = 0;
  _atom_count_window = -1;
  _calculations_bypassed_by_atom_count_window = 0;
  _smiles_id_dist = nullptr;
  _number_needles = 0;
  _tdts_read = 0;
  _find_all_matches = 0;
  _tabular_output = 1;
  _non_matches = 0;
  _exact_matches = 0;
  _output_separator = ' ';
  _calculation_type = CalculationType::kFindNearestNeedleNeighbour;
  _for_testing_stop = std::numeric_limits<uint32_t>::max();
}

Options::~Options() {
  if (_needles != nullptr) {
    delete [] _needles;
  }
  if (_guard_fp != nullptr) {
    delete [] _guard_fp;
  }
  if (_already_written != nullptr) {
    delete [] _already_written;
  }
  if (_smiles_id_dist != nullptr) {
    delete [] _smiles_id_dist;
  }
}

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (! cl.option_present('p')) {
    cerr << "Must specify needles via the -p option\n";
    Usage(1);
  }

  if (cl.option_present('t')) {
    if (! cl.value('t', _threshold) || _threshold <= 0.0 || _threshold >= 1.0) {
      cerr << "Options::Initialise:invalid threshold (-t)\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Threshold set to " << _threshold << '\n';
    }
  }

  if (cl.option_present('C')) {
    const_IWSubstring c = cl.option_value('C');
    if (c == '1') {
        _calculation_type = CalculationType::kIdentifyNeedlesNearHaystack;
    } else if (c == '2') {
        _calculation_type = CalculationType::kFindNearestNeedleNeighbour;
    } else {
      cerr << "Unrecognised calculation type '" << c << "'\n";
      return 0;
    }
  }

  if (cl.option_present('n')) {
    _tabular_output = 0;
    if (_verbose) {
      cerr << "Will generate output for nnprocess\n";
    }
  }

  // If we are looking for nearest neighbours of each haystack fingerprint,
  // and there is no threshold, set it to 1.0.
  if (_calculation_type == CalculationType::kFindNearestNeedleNeighbour &&
      _threshold == 0.0f) {
    _threshold = 1.0f;
  }

  // Make sure to parse this before reading the needles.
  if (cl.option_present('G')) {
    const_IWSubstring g;
    for (int i = 0; cl.value('G', g, i); ++i) {
      if (g.starts_with("thr=")) {
        g.remove_leading_chars(4);
        if (! g.numeric_value(_guard_fp_distance) || _guard_fp_distance <= 0.0 ||
            _guard_fp_distance > 1.0f) {
          cerr << "Options::Initialise:invalid guard fingerprint distance '" << g << "'\n";
          return 0;
        }
      } else {
        _guard_fp_tag = g;
      }
    }

    if (_guard_fp_tag.empty()) {
      cerr << "Guard fingerprint tag not specified\n";
      return 0;
    }

    if (_guard_fp_distance == 0.0f) {
      cerr << "Guard fingerprint distance threshold not specified\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Guard fingerprint in " << _guard_fp_tag << " distance " << _guard_fp_distance << '\n';
    }
  }

  if (cl.option_present('p')) {
    IWString fname = cl.string_value('p');
    if (! ReadNeedles(fname)) {
      cerr << "Cannot read needles from '" << fname << "'\n";
      return 0;
    }
  }

  if (cl.option_present('w')) {
    if (! cl.value('w', _atom_count_window) || _atom_count_window < 0) {
      cerr << "The atom count window (-w) must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will only compare fp's where the atom count differs by " << _atom_count_window << " or less\n";
    }
  }

  if (cl.option_present('r')) {
    if (! _report_progress.initialise(cl, 'r', _verbose)) {
      cerr << "Options::Initialise:cannot initialise progress reporter (-r)\n";
      return 0;
    }
  }

  if (cl.option_present('a')) {
    _find_all_matches = 1;
    if (_verbose) {
      cerr << "Will find all matches, not just closest\n";
    }
  }

  if (cl.option_present('o')) {
    IWString tmp = cl.string_value('o');
     
    if (! char_name_to_char(tmp)) {
      cerr << "Options::Initialise:unrecobnised char '" << tmp << "'\n";
      return 0;
    }
    _output_separator = tmp[0];
  }

  if (cl.option_present('b')) {
    if (! cl.value('b', _for_testing_stop) || _for_testing_stop == 0) {
      cerr << "The for testing stop value (-b) must be a non zero whole number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will stop processing after " << _for_testing_stop << " items read\n";
    }
  }

  if (cl.option_present('X')) {
    const char * fname = cl.option_value('X');
    if (! _stream_for_not_matched.open(fname)) {
      cerr << "Options::Initialise:cannot open -X file '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will write ids not matched to '" << fname << "'\n";
    }
  }

  return 1;
}

int
Options::ReadNeedles(IWString& fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::ReadNeedles:cannot open '" << fname << "'\n";
    return 0;
  }

  std::unique_ptr<RE2> rx = std::make_unique<RE2>("^PCN");
  _number_needles = input.grep(*rx);
  if (_number_needles <= 0) {
    cerr << "Options::ReadNeedles:no data in '" << fname << "'\n";
    return 0;
  }

  cerr << "Grep determined " << _number_needles << " needles in '" << fname << "'\n";
  _needles = new GFP_Standard[_number_needles];
  _already_written = new_int(_number_needles);
  _smiles_id_dist = new Smiles_ID_Dist[_number_needles];
  _active_needles = _number_needles;

  cerr << "Guard tag " << _guard_fp_tag << '\n';
  if (!_guard_fp_tag.empty()) {
    _guard_fp = new IWDYFP[_number_needles];
    cerr << "Allocated guard " << _number_needles << '\n';
  }

  return ReadNeedles(input);
}

void
SetSmiles(const IW_TDT& tdt,
          Smiles_ID_Dist& sid) {
  IWString smiles;
  tdt.dataitem_value(kSmilesTag, smiles);
  sid.set_smiles(smiles);
}

int
Build(const IW_General_Fingerprint& gfp,
      GFP_Standard& destination) {
  destination.build_molecular_properties(gfp.molecular_properties_integer());

  static bool first_call = true;
  static std::array<int, 3> stdfp_index;

  using gfp::StdFpIndex;

  if (first_call) {
    if (! gfp::GetStandardFingerprintIndices(stdfp_index)) {
      return 0;
    }
    first_call = false;
  }

#ifdef DEBUG_BUILD
  cerr << "FP 0 " << gfp[0].nbits() << '\n';
  cerr << "FP 1 " << gfp[1].nbits() << '\n';
  cerr << "FP 2 " << gfp[2].nbits() << '\n';
  cerr << "FP 3 " << gfp[3].nbits() << '\n';
#endif
  destination.build_molecular_properties(gfp.molecular_properties_integer());
  destination.build_mk(gfp[stdfp_index[StdFpIndex::kMK]]);
  destination.build_mk2(gfp[stdfp_index[StdFpIndex::kMK2]]);
  destination.build_iw(gfp[stdfp_index[StdFpIndex::kIWfp]]);

  return 1;
}

int
BuildGuardFp(const IW_General_Fingerprint& gfp,
             IWDYFP& destination) {
  destination = gfp[3];  // Hard coded for now...
  return 1;
}

int
Options::ReadNeedles(iwstring_data_source& input) {
  IW_TDT tdt;
  int ndx = 0;
  for ( ; tdt.next(input); ++ndx) {
    IW_General_Fingerprint gfp;

    int fatal;
    if (!gfp.construct_from_tdt(tdt, fatal)) {
      cerr << "Cannot read previously selected fingerprint\n";
      return 0;
    }

    if (! Build(gfp, _needles[ndx])) {
      cerr << "Options::ReadNeedles:cannot parse tdt\n";
      cerr << tdt;
      return 0;
    }
    if (!_guard_fp_tag.empty() && ! BuildGuardFp(gfp, _guard_fp[ndx])) {
      cerr << "Cannot construct guard fingerprint\n";
      cerr << tdt;
      return 0;
    }
    _smiles_id_dist[ndx].set_id(gfp.id());
    SetSmiles(tdt, _smiles_id_dist[ndx]);
  }

  if (ndx != _number_needles) {
    cerr << "Options::ReadNeedles:mismatch btween grep " << _number_needles << " and what is read " << ndx << '\n';
    return 0;
  }

  if (_verbose) {
    cerr << "Read " << _number_needles << " needles\n";
  }

  return 1;
}

int
Options::Process(const char * fname,
                 IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::Process:cannot open '" << fname << "'\n";
    return 0;
  }

  return Process(input, output);
}

// Read haystack members from `input` and compare with _needles.
int
Options::Process(iwstring_data_source& input,
                 IWString_and_File_Descriptor& output) {
  IW_TDT tdt;
  while (tdt.next(input)) {
    ++_tdts_read;

    IW_General_Fingerprint gfp;

    int fatal;
    if (!gfp.construct_from_tdt(tdt, fatal)) {
      cerr << "Cannot read fingerprint\n";
      cerr << tdt;
      return 0;
    }

    if (! Process(gfp, tdt, output)) {
      return 0;
    }

    if (_tdts_read >= _for_testing_stop) {
      cerr << "Processing terminated via request at " << _tdts_read << '\n';
      return 1;
    }
  }
  return 1;
}

int
Options::Process(const IW_General_Fingerprint& gfp,
                 const IW_TDT& tdt,
                 IWString_and_File_Descriptor& output) {
  GFP_Standard std_fp;
  if (! Build(gfp, std_fp)) {
    cerr << "Options::Process:invalid fingerprint\n";
    return 0;
  }
  IWDYFP guard_fp;
  if (! _guard_fp_tag.empty() && ! BuildGuardFp(gfp, guard_fp)) {
    cerr << "Cannot build guard fp\n";
    return 0;
  }

  if (! Process(std_fp, guard_fp, gfp.id(), tdt, output)) {
    return 0;
  }

  if (_report_progress()) {
     ReportProgress(cerr);
  }

  return 1;
}

int
Options::ReportProgress(std::ostream& output) const {
  output << "Read " << _tdts_read << " tdts ";
  switch(_calculation_type) {
    case CalculationType::kIdentifyNeedlesNearHaystack:
      output << _active_needles << " of " << _number_needles << " active\n";
      break;
    case CalculationType::kFindNearestNeedleNeighbour:
      output << _non_matches << " non matches " << iwmisc::Fraction<float>(_non_matches, static_cast<int>((_tdts_read - _non_matches))) << '\n';
      break;
    default:
      return 0;
  }

  return 1;
}

const_IWSubstring
GetSmiles(const IW_TDT& tdt) {
  const_IWSubstring result;
  tdt.dataitem_value(kSmilesTag, result);
  return result;
}

void
Options::UpdateDistanceStats(float distance) {
  _distance.extra(distance);
  if (distance == 0.0f) {
    ++_exact_matches;
  }
}

int
Options::CanBeCompared(const GFP_Standard& fp1,
                       const GFP_Standard& fp2) {
  if (_atom_count_window < 0) {
    return 1;
  }

  int delta = std::abs(fp1.natoms() - fp2.natoms());
  if (delta <= _atom_count_window) {
    return 1;
  }

  ++_calculations_bypassed_by_atom_count_window;

  return 0;
}

int
Options::Process(GFP_Standard& gfp,
                 IWDYFP& guard_fp,
                 const IWString& id,
                 const IW_TDT& tdt,
                 IWString_and_File_Descriptor& output) {
  switch (_calculation_type) {
    case CalculationType::kIdentifyNeedlesNearHaystack:
      return IdentifyNeedlesNearHaystack(gfp, guard_fp, id, tdt, output);

    case CalculationType::kFindNearestNeedleNeighbour:
      return FindNearestNeedleNeighbour(gfp, guard_fp, id, tdt, output);

    default:
      cerr << "What kind of calculation is being done\n";
      return 0;
  }
}

int
Options::DoOutput(const IWString& smiles1, const IWString& id1,
         const IWString& smiles2, const IWString& id2,
         float distance,
         IWString_and_File_Descriptor& output) {
  UpdateDistanceStats(distance);

  if (_tabular_output) {
      output << smiles1 << _output_separator << id1 << _output_separator <<
                smiles2 << _output_separator << id2 << _output_separator
                << distance << '\n';
  } else {
    output << kSmilesTag << smiles1 << ">\n";
    output << kIdentifierTag << id1 << ">\n";
    output << kSmilesTag << smiles2 << ">\n";
    output << kIdentifierTag << id2 << ">\n";
    output << kDistanceTag << distance << ">\n";
    output << "|\n";
  }

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
Options::WithinRange(IWDYFP& g1,
                     IWDYFP& g2) {
  float distance = 1.0f - g1.tanimoto(g2);
  ++_calculations_attempted;
  if (distance <= _guard_fp_distance) {
    return 1;
  }

  ++_calculations_bypassed_by_guard;

  return 0;
}

// Here we identify items in _needles that are close to a member of
// the haystack, `gfp`.
int
Options::IdentifyNeedlesNearHaystack(GFP_Standard& gfp,
                 IWDYFP& guard_fp,
                 const IWString& id,
                 const IW_TDT& tdt,
                 IWString_and_File_Descriptor& output) {
  for (int i = 0; i < _number_needles; ++i) {
    if (_already_written[i]) {
      continue;
    }

#ifdef OPTIMIZATIONS_DO_NOT_WORK
    if (! CanBeCompared(_needles[i], gfp)) {
      continue;
    }
    if (_guard_fp_distance > 0.0f && ! WithinRange(_guard_fp[i], guard_fp)) {
      continue;
    }
#endif

    std::optional<float> dist = gfp.tanimoto_distance_if_less(_needles[i], _threshold);
    if (! dist) {
      continue;
    }

    _already_written[i] = 1;
    --_active_needles;
    DoOutput(_smiles_id_dist[i].smiles(), _smiles_id_dist[i].id(),
             GetSmiles(tdt), id, *dist, output);
    if (_find_all_matches) {
      continue;
    }

    if (_active_needles <= 0) {
      break;
    }
  }

  return 1;
}

// Identify the needle that is closest to `gfp`.
int
Options::FindNearestNeedleNeighbour(GFP_Standard& gfp,
                 IWDYFP& guard_fp,
                 const IWString& id,
                 const IW_TDT& tdt,
                 IWString_and_File_Descriptor& output) {
  float min_distance = 1.0f;
  int index_of_closest = -1;

  for (int i = 0; i < _number_needles; ++i) {
#ifdef OPTIMIZATIONS_DO_NOT_WORK
    if (! CanBeCompared(_needles[i], gfp)) {
      continue;
    }
    if (_guard_fp_distance > 0.0f && ! WithinRange(_guard_fp[i], guard_fp)) {
      continue;
    }
#endif

    float threshold = std::min(min_distance, _threshold);
    std::optional<float> dist = gfp.tanimoto_distance_if_less(_needles[i], threshold);
    if (! dist) {
      continue;
    }
    min_distance = *dist;
    index_of_closest = i;
  }

  if (index_of_closest < 0) {
    return HandleWriteNonMatch(id, tdt);
  }

  DoOutput(GetSmiles(tdt), id, _smiles_id_dist[index_of_closest].smiles(),
           _smiles_id_dist[index_of_closest].id(), min_distance, output);

  return 1;
}

// If we are doing a calculation on needles, write any that are unselected.
int
Options::WriteNonMatches() {
  if (_calculation_type != CalculationType::kIdentifyNeedlesNearHaystack) {
    return 1;
  }

  for (int i = 0; i < _number_needles; ++i) {
    if (_already_written[i]) {
      continue;
    }

    _stream_for_not_matched << _smiles_id_dist[i].smiles() << _output_separator <<
                               _smiles_id_dist[i].id() << '\n';
    _stream_for_not_matched.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

int
Options::HandleWriteNonMatch(const IWString& id, const IW_TDT& tdt) {
  ++_non_matches;

  if (_stream_for_not_matched.active()) {
    _stream_for_not_matched << id << '\n';
    _stream_for_not_matched.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _tdts_read << " tdts, number not matched " << _active_needles << '\n';
  output << _active_needles << " needles still active\n";
  output << _non_matches << " haystack items not matched\n";
  output << _exact_matches << " exact matches\n";
  if (_atom_count_window >= 0) {
    output << _calculations_bypassed_by_atom_count_window << " calculations bypassed by atom count window " << _atom_count_window << '\n';
  }
  if (! _guard_fp_tag.empty()) {
    output << _calculations_bypassed_by_guard << " calculations bypassed by guard " << _guard_fp_tag << " thr " << _guard_fp_distance << '\n';
    output << "Fraction bypassed " << iwmisc::Fraction<float>(_calculations_bypassed_by_guard, _calculations_attempted) << '\n';
  }
  if (_distance.n() > 0) {
    output << "For items selected (threshold " << _threshold << ")\n";
    output << "Distances btw " << _distance.minval() << " and " << _distance.maxval() << " mean " << _distance.average() << '\n';
  }
  return 1;
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vC:p:r:t:X:ao:w:G:b:n");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  if (cl.size() != 1) {
    cerr << "Sorry, haystack file must be the only argument\n";
    return 1;
  }

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  IWString_and_File_Descriptor output(1);
  if (! options.Process(cl[0], output)) {
    cerr << "Fatal error\n";
    return 1;
  }

  output.flush();

  options.WriteNonMatches();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace gfp_distance_filter_standard

int
main(int argc, char** argv)
{
  int rc = gfp_distance_filter_standard::Main(argc, argv);

  return rc;
}
