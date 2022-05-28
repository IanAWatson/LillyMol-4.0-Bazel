// Distance filter for gfp fingerprints

#include <stdlib.h>

#include <algorithm>
#include <iostream>
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
  cerr << " -r <n>           report progress every <n> haystack members processed\n";
  cerr << " -X <fname>       write needles not matched to <fname>\n";
  cerr << " -v               verbose output\n";
  ::exit(rc);
}

class Options {
  private:
    // We can do several different kinds of calculations.
    enum class CalculationType {
      // Identify needles that are within a given distance of any member of the haystack.
      // Uses _threshold as the distance.
      kIdentifyNeedlesNearHaystack,

      // For each member of the haystack, identify the nearest neighbour among _needles
      kFindNearestNeedleNeighbour,
    };

    CalculationType _calculation_type;

    float _threshold;

    int _verbose;

    GFP_Standard* _needles;
    int _number_needles;

    // GFP_Standard do not have a name field. Create one.
    IWString* _name;
    Smiles_ID_Dist* _smiles_id_dist;

    // Once a needle has been identified as within range, we can optionally
    // stop searching it.
    int * _already_written;
    // If we are stopping searching, we keep track of the number of active needles
    // and can stop when that reaches zero.
    int _active_needles;

    int _find_all_matches;

    Report_Progress _report_progress;

    int _tdts_read;

    // Of the calculations that return a value in range
    Accumulator<double> _distance;
    int _exact_matches;

    // For items that are not matched by the haystack, we can write them.
    IWString_and_File_Descriptor _stream_for_not_matched;

    int _non_matches;

    int _tabular_output;

    char _output_separator;

  // private functions
    int ReadNeedles(IWString& fname);
    int ReadNeedles(iwstring_data_source& input);
    int HandleWriteNonMatch(const IWString& id, const IW_TDT& tdt);
    int DoOutput(const IWString& smiles1, const IWString& id1,
         const IWString& smiles2, const IWString& id2,
         float distance,
         IWString_and_File_Descriptor& output);

    void UpdateDistanceStats(float distance);
    int IdentifyNeedlesNearHaystack(GFP_Standard& gfp,
                 const IWString& id,
                 const IW_TDT& tdt,
                 IWString_and_File_Descriptor& output);
    int FindNearestNeedleNeighbour(GFP_Standard& gfp,
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
    int Process(GFP_Standard& gfp,
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
  _already_written = nullptr;
  _active_needles = 0;
  _name = nullptr;
  _smiles_id_dist = nullptr;
  _number_needles = 0;
  _tdts_read = 0;
  _find_all_matches = 0;
  _tabular_output = 1;
  _non_matches = 0;
  _exact_matches = 0;
  _output_separator = ' ';
  _calculation_type = CalculationType::kFindNearestNeedleNeighbour;
}

Options::~Options() {
  if (_needles != nullptr) {
    delete [] _needles;
  }
  if (_already_written != nullptr) {
    delete [] _already_written;
  }
  if (_name != nullptr) {
    delete [] _name;
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

  // If we are looking for nearest neighbours of each haystack fingerprint,
  // and there is no threshold, set it to 1.0.
  if (_calculation_type == CalculationType::kFindNearestNeedleNeighbour &&
      _threshold == 0.0f) {
    _threshold = 1.0f;
  }

  if (cl.option_present('p')) {
    IWString fname = cl.string_value('p');
    if (! ReadNeedles(fname)) {
      cerr << "Cannot read needles from '" << fname << "'\n";
      return 0;
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
  destination.build_iw(gfp[0]);
  destination.build_mk(gfp[1]);
  destination.build_mk2(gfp[2]);

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

    GFP_Standard std_fp;
    if (! Build(gfp, std_fp)) {
      cerr << "Options::Process:invalid fingerprint\n";
      return 0;
    }

    if (! Process(std_fp, gfp.id(), tdt, output)) {
      return 0;
    }

    if (_report_progress()) {
      ReportProgress(cerr);
    }
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
      output << _non_matches << " non matches\n";
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
Options::Process(GFP_Standard& gfp,
                 const IWString& id,
                 const IW_TDT& tdt,
                 IWString_and_File_Descriptor& output) {
  switch (_calculation_type) {
    case CalculationType::kIdentifyNeedlesNearHaystack:
      return IdentifyNeedlesNearHaystack(gfp, id, tdt, output);

    case CalculationType::kFindNearestNeedleNeighbour:
      return FindNearestNeedleNeighbour(gfp, id, tdt, output);

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

// Here we identify items in _needles that are close to a member of
// the haystack, `gfp`.
int
Options::IdentifyNeedlesNearHaystack(GFP_Standard& gfp,
                 const IWString& id,
                 const IW_TDT& tdt,
                 IWString_and_File_Descriptor& output) {
  for (int i = 0; i < _number_needles; ++i) {
    if (_already_written[i]) {
      continue;
    }

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
                 const IWString& id,
                 const IW_TDT& tdt,
                 IWString_and_File_Descriptor& output) {
  float min_distance = 1.0f;
  int index_of_closest = -1;

  for (int i = 0; i < _number_needles; ++i) {
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
  if (_distance.n() > 0) {
    output << "For items selected (threshold " << _threshold << ")\n";
    output << "Distances btw " << _distance.minval() << " and " << _distance.maxval() << " mean " << _distance.average() << '\n';
  }
  return 1;
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vC:p:r:t:X:ao:");
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
