/*
  Spread implementation with fixed fingerprints MPI
*/
#include <stdlib.h>

#include <algorithm>
#include <iostream>
#include <random>

#include <mpich/mpi.h>

#define REPORT_PROGRESS_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/numeric_data_from_file.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "gfp_standard.h"
#include "smiles_id_dist.h"
#include "sparse_collection.h"

namespace gfp_spread_standard_mpi {

using std::cerr;

void
Usage(int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
  cerr << "Spread implementation requiring MPR IW MK MK2\n";
  cerr << " -n <nsel>       number of items to select\n";
  cerr << " -S rand         choose first item randomly\n";
  cerr << " -s <size>       size of fingerprint file\n";
  cerr << " -t <dist>       stop selection once distance drops below <dist>\n";
  cerr << " -r <n>          report progress every <n> items selected\n";
  cerr << " -q <n>          squeeze out already selected items evern <n> selections (def 1000)\n";
  cerr << " -b              brief output (smiles id dist)\n";
  cerr << " -v              verbose output\n";

  ::exit(rc);
}

IWString smiles_tag("$SMI<");
IWString identifier_tag("PCN<");
IWString distance_tag("DIST<");
IWString mk_tag("FPMK<");
IWString mk2_tag("FPMK2<");
IWString iw_tag("FPIW<");

Fraction_as_String fraction_as_string;

/*
  Each item consists of a fingerprint and information about its state.
*/

class SSpread_Item : public GFP_Standard {
 private:
  // Whether or not this has been selected.
  int _selected;

  // Closest distance to something previously selected.
  float _distance;

  // The _initial_ndx of the closest previously selected item.
  int _nearest_previously_selected;

  // Within the input file, the index. This is used for
  // indexing into the _smiles and _pcn arrays.
  int _initial_ndx;

 public:
  SSpread_Item();

  int selected() const {
    return _selected;
  }
  void set_selected(int s) {
    _selected = s;
  }

  float distance() const {
    return _distance;
  }

  int nearest_previously_selected() const {
    return _nearest_previously_selected;
  }

  void IsFirstSelected();

  void SetDistance(SSpread_Item& sel) {
    _distance = GFP_Standard::tanimoto_distance(sel);
    _nearest_previously_selected = sel.initial_ndx();
  }

  SSpread_Item& operator=(const SSpread_Item& rhs);
  // In retrospect, there seems to be no advantage to having
  // a move assignment operator.
  SSpread_Item& operator=(SSpread_Item&& rhs);

  void
  set_initial_ndx(int s) {
    _initial_ndx = s;
  }
  int
  initial_ndx() const {
    return _initial_ndx;
  }

  int SelectedOrZeroDistance() const;

  int NofitySelected(const SSpread_Item& rhs) {
    float d = GFP_Standard::tanimoto_distance(rhs);
    if (d >= _distance) {
      return 0;
    }
    _distance = d;
    _nearest_previously_selected = rhs.initial_ndx();
    return 1;
  }

  int NotifySelectedUseThreshold(const SSpread_Item& rhs) {
    std::optional<float> d = GFP_Standard::tanimoto_distance(rhs, _distance);
    if (! d) {
      return 0;
    }
    _distance = *d;
    _nearest_previously_selected = rhs.initial_ndx();
    return 1;
  }


  const IWString&
  smiles(const IWString* smiles) const {
    return smiles[_initial_ndx];
  }
  const IWString&
  pcn(const IWString* pcn) const {
    return pcn[_initial_ndx];
  }

  int DebugPrint(std::ostream& output) const;
};

SSpread_Item::SSpread_Item() {
  _selected = 0;
  _distance = 1.0;
  _initial_ndx = -1;
  _nearest_previously_selected = -1;
}

SSpread_Item&
SSpread_Item::operator=(const SSpread_Item& rhs) {
  _selected = rhs._selected;
  _distance = rhs._distance;
  _nearest_previously_selected = rhs._nearest_previously_selected;
  _initial_ndx = rhs._initial_ndx;
  GFP_Standard* me = this;
  const GFP_Standard* other = &rhs;
  *me = *other;
  return *this;
}

SSpread_Item&
SSpread_Item::operator=(SSpread_Item&& rhs) {
  _initial_ndx = rhs._initial_ndx;
  _selected = rhs._selected;
  _nearest_previously_selected = rhs._nearest_previously_selected;
  _distance = rhs._distance;
  GFP_Standard* me = this;
  const GFP_Standard* other = &rhs;
  *me = std::move(*other);
  return *this;
}

void
SSpread_Item::IsFirstSelected() {
  _selected = 1;
  _distance = 1.0f;
  _nearest_previously_selected = -1;
}

int
SSpread_Item::SelectedOrZeroDistance() const {
  if (_selected || 0.0f == _distance) {
    return 1;
  }

  return 0;
}

int
SSpread_Item::DebugPrint(std::ostream& output) const {
  output << "SSpread_Item:item " << _initial_ndx << " sel? " << _selected << 
            " dist " << _distance << " nsn " << _nearest_previously_selected << '\n';
  return 1;
}

// A structure to send info from leader to workers
// For each chunk, the ndx of the farthest item.
// Note that we make assumptions below that this is size = 3 * sizeof(int).

struct NdxDist {
  int ndx = -1;
  float distance = -1.0f;
  int nearest_previously_selected = -1;
};

class Options {
  private:
    int _verbose;
    int _choose_first_point_randomly;
    int _nsel;
    float _stop_once_distance_drops_below;

    int _squeeze = 0;
    int _next_squeeze = 0;

    int _brief_output = 0;

    Report_Progress _report_progress;

    int _pool_size;
    SSpread_Item* _pool;

    IWString* _smiles;
    IWString* _pcn;

    // Private functions;

    int AllocatePool();
    int ReadPool(iwstring_data_source& input);
    int ReadFingerprints(iwstring_data_source& input);
    int WriteSelectedBrief( int selected, IWString_and_File_Descriptor& output);
    int WriteSelected(int selected,
                      const struct NdxDist& data,
                      IWString_and_File_Descriptor& output);
    int DoSqueeze(int next_selected);

  public:
    Options();
    ~Options();

    int Initialise(Command_Line& cl, int world_rank);

    int ReadFingerprints(const char * fname);

    int Spread(IWString_and_File_Descriptor& output);

    int DoWorker();
};

Options::Options() {
  _verbose = 0;
  _choose_first_point_randomly = 0;
  _nsel = std::numeric_limits<int>::max();
  _stop_once_distance_drops_below = 0.0;

  _brief_output = 0;

  _pool_size = 0;
  _pool = nullptr;

  _smiles = nullptr;
  _pcn = nullptr;
}

Options::~Options() {
  if (_smiles != nullptr) {
    delete [] _smiles;
  }
  if (_pcn != nullptr) {
    delete [] _pcn;
  }
  if (_pool != nullptr) {
    delete [] _pool;
  }
}

int
Options::Initialise(Command_Line& cl, int world_rank) {
  _verbose = cl.option_count('v');

  if (cl.option_present('S')) {
    const_IWSubstring s = cl.string_value('S');

    if ("rand" == s) {
      _choose_first_point_randomly = 1;

      if (_verbose)
        cerr << "Will choose first item randomly\n";
    } else {
      cerr << "Unrecognised -S qualifier '" << s << "'\n";
      return 0;
    }
  }

  if (cl.option_present('n')) {
    if (!cl.value('n', _nsel) || _nsel < 2) {
      cerr << "The number if fingerprints to select must be a whole +ve number\n";
      Usage(2);
    }
    if (_verbose) {
      cerr << "Will select " << _nsel << " items\n";
    }
  }

  if (cl.option_present('t')) {
    if (!cl.value('t', _stop_once_distance_drops_below) || 
                  _stop_once_distance_drops_below < 0.0f ||
                  _stop_once_distance_drops_below >= 1.0f) {
      cerr << "The stop selection distance option (-t) must be a valid distance\n";
      return 0;
    }

    if (_verbose)
      cerr << "Will stop selection once distance drops below " << _stop_once_distance_drops_below
           << '\n';
  }

  if (cl.option_present('q')) {
    if (!cl.value('q', _squeeze) || _squeeze < 0) {
      cerr << "The squeeze every option (-q) must be a whole non negative number\n";
      Usage(1);
    }

    if (_verbose)
      cerr << "Will squeeze selected items every " << _squeeze << " items selected\n";

    if (_squeeze == 0) {
      _next_squeeze = -1;
    } else {
      _next_squeeze = _squeeze;
    }
  } else {
    _squeeze = 1000;
    _next_squeeze = 1000;
  }

  if (cl.option_present('s')) {
    if (!cl.value('s', _pool_size) || _pool_size < 2) {
      cerr << "The pool size specification (-s) must be a whole +ve number\n";
      Usage(2);
    }

    if (_verbose)
      cerr << "Problem sized for " << _pool_size << " fingerprints\n";
  }

  if (! AllocatePool()) {
    return 0;
  }

  if (cl.option_present('b')) {
    _brief_output = 1;

    if (_verbose)
      cerr << "Brief output\n";
  }

  if (world_rank == 0 && cl.option_present('r')) {
    if (!_report_progress.initialise(cl, 'r', _verbose)) {
      cerr << "Cannot initialise report progress option (-r)\n";
      Usage(2);
    }
  }

  return 1;
}


int
RandomPoolMember(int pool_size)
{
  std::random_device rd;
  std::default_random_engine generator(rd());
  std::uniform_int_distribution<int> u(0, pool_size - 1);
  return u(generator);
}

int
Options::DoSqueeze(int next_selected) {
  if (_pool_size < 1000) {
    return next_selected;
  }

  int ndx = 0;
  int rc = -1;
  for (int i = 0; i < _pool_size; ++i) {
    if (_pool[i].selected()) {
      continue;
    }

    if (i == next_selected) {
      rc = ndx;
    }

    if (i == ndx) {
      ++ndx;
      continue;
    }

    _pool[ndx] = _pool[i];
    ++ndx;
  }
  
  _pool_size = ndx;
  _next_squeeze += _squeeze;

  return rc;
}

int
Options::WriteSelectedBrief( int selected,
                   IWString_and_File_Descriptor& output) {
  const SSpread_Item& sel = _pool[selected];
  output << sel.smiles(_smiles) << ' ' << sel.pcn(_pcn) << ' ' << sel.distance() << '\n';
  output.write_if_buffer_holds_more_than(4096);
  return 1;
}

int
Options::WriteSelected(int selected,
              const struct NdxDist& data,
              IWString_and_File_Descriptor& output) {
  if (selected < 0) {
    return 1;
  }
  if (_brief_output)
    return WriteSelectedBrief(selected, output);

  const SSpread_Item& sel = _pool[selected];

  output << smiles_tag << sel.smiles(_smiles) << ">\n";
  output << identifier_tag << sel.pcn(_pcn) << ">\n";

  if (data.nearest_previously_selected >= 0) {
    int ndx = data.nearest_previously_selected;
    output << smiles_tag << _smiles[ndx] << ">\n";
    output << identifier_tag << _pcn[ndx] << ">\n";
    output << distance_tag << data.distance << ">\n";
  } else {
    output << smiles_tag << "*>\n";
    output << identifier_tag << "*>\n";
    output << distance_tag << "1>\n";
  }

  output << "|\n";

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
Options::Spread(IWString_and_File_Descriptor& output) {
  int sel;

  if (_choose_first_point_randomly) {
    sel = RandomPoolMember(_pool_size);
  } else {
    sel = 0;
  }

  constexpr int kManager = 0;

  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Scope here for efficiency.
  const int tag = 0;
  NdxDist data;

  for (int items_selected = 0; items_selected < _nsel; items_selected++) {
    WriteSelected(sel, data, output);
#ifdef DEBUG_MANAGER
    cerr << "Manager broadcasgint selected item " << sel << '\n';
#endif
    if (auto rc = MPI_Bcast(&sel, 1, MPI_INT, kManager, MPI_COMM_WORLD);
        rc != MPI_SUCCESS) {
      cerr << "Options::Spread:bcast failed " << rc << '\n';
      return 0;
    }
#ifdef DEBUG_MANAGER
    cerr << "Broadcast sent\n";
#endif
    if (sel < 0) {
      break;
    }

    float max_distance = -1.0f;
    int nexts = -1;
    for (int i = 1; i < world_size; ++i) {
      if (auto rc = MPI_Recv(&data, 3, MPI_INT, i, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          rc != MPI_SUCCESS) {
        cerr << "Options::Spread:Recv from " << i << " failed " << rc << '\n';
        return 0;
      }
      if (data.distance > max_distance) {
        nexts = data.ndx;
        max_distance = data.distance;
      }
    }
    sel = nexts;
    if (items_selected == _next_squeeze) {
      sel = DoSqueeze(sel);
    }
  }

  return 1;
}


int
Options::ReadFingerprints(const char * fname) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Options::ReadFingerprints:cannot open fingerprint file '" << fname << "'\n";
    return 0;
  }

  return ReadFingerprints(input);
}

int
Options::ReadFingerprints(iwstring_data_source& input) {
  if (0 == _pool_size) {
    _pool_size = input.count_records_starting_with(identifier_tag);

    if (0 == _pool_size) {
      cerr << "No occurrences of " << identifier_tag << "' in input\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Job automatically sized to " << _pool_size << " fingerprints\n";
    }
  }

  if (!AllocatePool()) {
    return 0;
  }

  return ReadPool(input);
}

int
Options::ReadPool(iwstring_data_source& input) {
  IW_TDT tdt;

  int ndx = 0;

  for (; tdt.next(input) && ndx < _pool_size; ndx++) {
    _pool[ndx].set_initial_ndx(ndx);

    tdt.dataitem_value(smiles_tag, _smiles[ndx]);

    tdt.dataitem_value(identifier_tag, _pcn[ndx]);

    IW_General_Fingerprint gfp;

    int fatal;
    if (!gfp.construct_from_tdt(tdt, fatal)) {
      cerr << "Cannot read fingerprint\n";
      return 0;
    }

    if (0 == ndx) {
      if (!standard_fingerprints_present())
        return 0;
    }

    _pool[ndx].build_molecular_properties(gfp.molecular_properties_integer());
    _pool[ndx].build_iw(gfp[0]);
    _pool[ndx].build_mk(gfp[1]);
    _pool[ndx].build_mk2(gfp[2]);
  }


  for (int i = 0; i < ndx; ++i) {
    if (_pool[i].initial_ndx() != i) {
      cerr << "Initial index wrong\n";
    }
  }

  if (ndx < 2) {
    cerr << "Yipes, did not read enough fingerprints\n";
    return 0;
  }

  _pool_size = ndx;

  return _pool_size;
}

int
Options::AllocatePool() {
  _pool = new SSpread_Item[_pool_size];
  _smiles = new IWString[_pool_size];
  _pcn = new IWString[_pool_size];

  if (nullptr == _pool || nullptr == _smiles || nullptr == _pcn) {
    cerr << "Yipes, could not allocate " << _pool_size << " fingerprints\n";
    return 0;
  }

  return 1;
}

int
Options::DoWorker() {

  NdxDist data;  // scope here for efficiency.
  constexpr int kManager = 0;

  for (int number_selected = 0; ;++number_selected) {
    int sel;
#ifdef DEBUG_WORKER
    cerr << "Worker waiting for identify of next sel\n";
#endif
    if (auto rc = MPI_Bcast(&sel, 1, MPI_INT, kManager, MPI_COMM_WORLD);
        rc != MPI_SUCCESS) {
      cerr << "Options::DoWorker:cannot receive next selected " << rc << '\n';
      return 0;
    }
#ifdef DEBUG_WORKER
    cerr << "Woroker got notice thta " << sel << " is selected\n";
#endif

    if (sel < 0) {  // Signal from manager that we are done.
      break;
    }

    SSpread_Item& selected = _pool[sel];
    selected.set_selected(1);

    int next_sel = -1;
    float max_distance = -1.0f;
    for (int i = 0; i < _pool_size; ++i) {
      if (_pool[i].SelectedOrZeroDistance()) {
        continue;
      }

      _pool[i].NotifySelectedUseThreshold(selected);
      float d = _pool[i].distance();
      if (d > max_distance) {
        max_distance = d;
        next_sel = i;
      }
    }

    // Sending negative values back is OK, signals we are done.
    data.ndx = next_sel;
    data.distance = max_distance;
    data.nearest_previously_selected = _pool[next_sel].nearest_previously_selected();

#ifdef DEBUG_WORKER
    cerr << "Worker sending " << data.ndx << " as next sel, dist " << max_distance << '\n';
#endif
    if (auto rc = MPI_Send(&data, 3, MPI_INT, kManager, 0, MPI_COMM_WORLD);
        rc != MPI_SUCCESS) {
      cerr << "Options::DoWorker:cannot send data to manager " << rc << '\n';
      return 0;
    }

    if (number_selected >= _next_squeeze) {
      sel = DoSqueeze(sel);  // update to `s` not needed, just to preserve api.
    }
  }

  return 1;
}

void
DoWorker(Options& options) {
  options.DoWorker();
  MPI_Finalize();
}

int
GFPSpreadStandard(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vA:s:r:n:bh:p:S:N:t:q:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  MPI_Init(NULL, NULL);

  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  if (world_size < 2) {
    cerr << "World size must be greater than 1, got " << world_size << '\n';
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  Options options;
  if (! options.Initialise(cl, world_rank)) {
    cerr << "Cannot initialise\n";
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(2);
  }

  if (cl.number_elements() > 1) {
    cerr << "Sorry, cannot handle multiple input files, concatenate them and try again\n";
    Usage(2);
  }

  set_report_fingerprint_status(0);

  if (! options.ReadFingerprints(cl[0])) {
    cerr << "Cannot read fingerprints '" << cl[0] << "'\n";
    return 1;
  }

  if (world_rank != 0) {
    DoWorker(options);
    return 0;
  }

  fraction_as_string.initialise(0.0, 1.0, 4);
  fraction_as_string.set_leading_string("DIST<");
  fraction_as_string.append_to_each_stored_string(">\n");

  IWString_and_File_Descriptor output(1);

  if (! options.Spread(output)) {
    cerr << "Error during selection\n";
    return 1;
  }

  delete_gfp_file_scope_static_objects();

  MPI_Finalize();

  return 0;
}

}  // namespace gfp_spread_standard_mpi

int
main(int argc, char** argv)
{
  int rc = gfp_spread_standard_mpi::GFPSpreadStandard(argc, argv);

  return rc;
}
