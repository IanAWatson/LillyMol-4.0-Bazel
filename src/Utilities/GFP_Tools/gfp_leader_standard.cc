/*
  Leader implementation with fixed fingerprints
*/

#include <stdlib.h>
#include <iostream>
#include <limits>
#include <memory>

#include <omp.h>

#include "google/protobuf/text_format.h"
#include "google/protobuf/io/zero_copy_stream.h"
#include "google/protobuf/io/zero_copy_stream_impl.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/report_progress.h"

#include "Utilities/GFP_Tools/nearneighbors.pb.h"

#include "gfp_standard.h"
#include "smiles_id_dist.h"
#include "sparse_collection.h"

using std::cerr;
using std::endl;
using std::cout;

using iw_tf_data_record::TFDataWriter;

const char * prog_name = nullptr;

static int verbose = 0;

static int clusters_to_form = std::numeric_limits<int>::max();

static Report_Progress report_progress;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString distance_tag("DIST<");
static IWString mk_tag("FPMK<");
static IWString mk2_tag("FPMK2<");
static IWString iw_tag("FPIW<");

static float threshold = 0.0;

class Leader_Item : public Smiles_ID_Dist
{
  private:
    const GFP_Standard * _gfp;

  public:
    Leader_Item();

    void set_gfp(const GFP_Standard * s) { _gfp = s;}
};

Leader_Item::Leader_Item()
{
  _gfp = nullptr;

  return;
}

static int pool_size = 0;

static Leader_Item * leader_item = nullptr;
static GFP_Standard * fingerprints = nullptr;
static int * selected = nullptr;
static float * distances = nullptr;

static void
usage(int rc)
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Leader implementation requiring MPR IW MK MK2\n";
  cerr << " -t <dist>      distance threshold\n";
  cerr << " -n <number>    number of clusters to form\n";
  cerr << " -A <fname>     file of previously selected fingerprints\n";
  cerr << " -s <size>      size of fingerprint file\n";
  cerr << " -p <n>         process <n> leaders at once (max 2 for now)\n";
  cerr << " -h <n>         maximum number of OMP threads to use\n";
  cerr << " -r <n>         report progress every <n> clusters formed\n";
  cerr << " -Y ...         output as proto <ascii,bin>\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

static int
choose_next_leader(const int * selected,
                   int n)
{
  for (int i = 0; i < n; i++)
  {
//      cerr << " CNL i = " << i << " sel " << selected[i] << endl;
    if (! selected[i])
      return i;
  }

  return -1;
}

static int
choose_next_leaders(const int * selected,
                    int n,
                    int * leaders,
                    int leaders_requested)
{
  int leader_count = 0;

  for (int i = 0; i < n; i++)
  {
    if (selected[i])
      continue;

    leaders[leader_count] = i;
    leader_count++;

    if (leader_count == leaders_requested)
      return leaders_requested;
  }

  return leader_count;
}

static int
leader()
{
  int cluster_id = 0;

  int icentre;

  int start = 0;

#ifdef DEBUG_LEADER1
  for (int i = 0; i < pool_size; i++)
  {
     cerr << " i = " << i << " selected " << selected[i] << endl;
  }
#endif

  while ((icentre = choose_next_leader(selected + start, pool_size - start)) >= 0)
  {
    cluster_id++;

    if (report_progress())
      cerr << "Formed " << cluster_id << " clusters, " << count_non_zero_occurrences_in_array(selected, pool_size) << " items in clusters\n";

#ifdef DEBUG_LEADER1
    cerr << "Selected " << icentre << " index " << (icentre+start) << endl;
#endif

    icentre += start;

    selected[icentre] = cluster_id;
    distances[icentre] = 0.0f;

    start = icentre + 1;

#pragma omp parallel for schedule(dynamic,256)
#ifdef DEBUG_LEADER1
    cerr << "leader was " << icentre << " scanning from " << start << " to " << pool_size << endl;
#endif
    for(int i= start;i<pool_size;++i)
    {
      if (!selected[i])
      {
        distances[i] = fingerprints[icentre].tanimoto_distance(fingerprints[i]);
        if (distances[i] <= threshold) selected[i] = cluster_id;
      }
    }

    if (cluster_id >= clusters_to_form)
      break;
  }

#ifdef DEBUG_LEADER1
  for (int i = 0; i < pool_size; i++)
  {
    cerr << " i = " << i << " selected " << selected[i] << endl;
  }
#endif

  return cluster_id;
}

static int
leader(int * leaders,
       int leaders_requested,
       float * cand_distances)
{
  int cluster_id = 0;

  int leaders_found;
  int start = 0;

  while ((leaders_found = choose_next_leaders(selected + start, pool_size - start, leaders, leaders_requested)))
  {
    for (int i = 0; i < leaders_found; i++)
    {
      leaders[i] += start;
    }
    if (1 == leaders_found || fingerprints[leaders[0]].tanimoto_distance(fingerprints[leaders[1]]) <= threshold)
    {
      cluster_id++;

      int icentre = leaders[0];
      start = icentre + 1;
      selected[icentre] = cluster_id;
      distances[icentre] = 0.0f;
#pragma omp parallel for schedule(dynamic,256)
      for(int i= start;i<pool_size;++i)
      {
        if (!selected[i])
        {
          distances[i] = fingerprints[icentre].tanimoto_distance(fingerprints[i]);
          if (distances[i] <= threshold) selected[i] = cluster_id;
        }
      }
    }
    else
    {
      // Both of these will end up being actual centers so amortize FPDB load costs by comparing
      // everything else to both
      const int c0id = cluster_id+1;
      const int c1id = cluster_id+2;
      cluster_id += 2;
      start = leaders[1] + 1;
      //cout << "S"<<c0id << " " << pool.id(potl_centers[0]) << " " << 0.0f<<endl;
      //cout << "S"<<c1id << " " << pool.id(potl_centers[1]) << " " << 0.0f<<endl;
      selected[leaders[0]] = c0id;
      selected[leaders[1]] = c1id;
      distances[leaders[0]] = 0.0f;
      distances[leaders[1]] = 0.0f;
      #pragma omp parallel for schedule(dynamic,256) private(cand_distances)
      for (int i = start; i < pool_size; i++)
      {
        if (!selected[i])
        {
          fingerprints[i].tanimoto_distance_2(fingerprints[leaders[0]], fingerprints[leaders[1]], cand_distances);
          if (cand_distances[0] <= threshold)
          {
            selected[i]  = c0id;
            distances[i] = cand_distances[0];
          }
          else if (cand_distances[1] <= threshold)
          {
            selected[i]  = c1id;
            distances[i] = cand_distances[1];
          }
        }
      }
    }
  }

  return cluster_id;
}

static int
do_previously_selected(const IW_General_Fingerprint & gfp,
                       float threshold)
{
  GFP_Standard sgfp;

  sgfp.build_molecular_properties(gfp.molecular_properties_integer());
  sgfp.build_mk(gfp[kMkFp]);
  sgfp.build_mk2(gfp[kMk2Fp]);
  sgfp.build_iw(gfp[kIwFp]);

  for (int i = 0; i < pool_size; i++)
  {
    if (selected[i])
      continue;

    float d = 1.0f - sgfp.tanimoto(fingerprints[i]);

    if (d > threshold)
      continue;

    selected[i] = -1;
  }

  return 1;
}

static int
do_previously_selected(iwstring_data_source & input,
                      float threshold)
{
  IW_TDT tdt;

  while (tdt.next(input))
  {
    IW_General_Fingerprint gfp;

    int fatal;
    if (! gfp.construct_from_tdt(tdt, fatal))
    {
      cerr << "Cannot read previously selected fingerprint\n";
      return 0;
    }

    do_previously_selected(gfp, threshold);
  }

  return 1;
}

static int
do_previously_selected(const char * fname,
                       float threshold)
{
  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open previously selected file '" << fname << "'\n";
    return 0;
  }

  return do_previously_selected(input, threshold);
}

void
TransferSmilesAndId(const IW_TDT& tdt,
                    Leader_Item& destination) {
  IWString tmp;

  tdt.dataitem_value(smiles_tag, tmp);
  destination.set_smiles(tmp);

  tdt.dataitem_value(identifier_tag, tmp);
  destination.set_id(tmp);
}

static int
read_pool(iwstring_data_source & input)
{
  IW_TDT tdt;

  IWString tmp;

  int ndx = 0;

  std::array<int, 3> stdfp_index;
  int fatal;

  using gfp::StdFpIndex;

  for (;tdt.next(input) && ndx < pool_size; ndx++) {
    leader_item[ndx].set_gfp(fingerprints+ndx);

    TransferSmilesAndId(tdt, leader_item[ndx]);

    IW_General_Fingerprint gfp;
    if (! gfp.construct_from_tdt(tdt, fatal)) {
      cerr << "Cannot read fingerprint\n";
      return 0;
    }

    if (0 == ndx && ! gfp::GetStandardFingerprintIndices(stdfp_index)) {
        return 0;
    }

    fingerprints[ndx].build_molecular_properties(gfp.molecular_properties_integer());
    fingerprints[ndx].build_mk(gfp[stdfp_index[StdFpIndex::kMK]]);
    fingerprints[ndx].build_mk2(gfp[stdfp_index[StdFpIndex::kMK2]]);
    fingerprints[ndx].build_iw(gfp[stdfp_index[StdFpIndex::kIWfp]]);
  }

  if (ndx < 2) {
    cerr << "Yipes, did not read enough fingerprints\n";
    return 0;
  }

  pool_size = ndx;

  return ndx;
}

static int
allocate_pool(int s)
{
  leader_item = new Leader_Item[s];
  selected = new_int(s);
  fingerprints = new GFP_Standard[s];
  distances = new float[s];

  if (nullptr == leader_item || nullptr == selected || nullptr == fingerprints || nullptr == distances)
  {
    cerr << "Yipes, could not allocate " << pool_size << " fingerprints\n";
    return 0;
  }

  pool_size = s;

  return 1;
}

int
WriteTextProto(NearNeighbors::Neighbors& proto,
               IWString_and_File_Descriptor& output) {
  google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);
  std::string destination;
  printer.PrintToString(proto, &destination);
  output << destination << '\n';
#ifdef OLD
  using google::protobuf::io::ZeroCopyOutputStream;
  using google::protobuf::io::FileOutputStream;
  std::unique_ptr<ZeroCopyOutputStream> zero_copy_output(new FileOutputStream(output.fd()));
  if (! google::protobuf::TextFormat::Print(proto, zero_copy_output.get(), true)) {
    std::cerr << "WriteProtoAsText:cannot write\n";
    return 0;
  }
#endif

  return 1;
}

int
WriteBinaryProto(NearNeighbors::Neighbors& proto,
                 TFDataWriter& writer) {
  return writer.WriteSerializedProto(proto);
}

static int
gfp_leader_standard(int argc, char ** argv)
{
  Command_Line cl(argc, argv, "vt:A:s:p:r:h:n:Y:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  set_report_fingerprint_status(0);

  if (! cl.option_present('t'))
  {
    cerr << "Must specify cluster threshold via the -t option\n";
    usage(2);
  }

  if (cl.option_present('t'))
  {
    if (! cl.value('t', threshold) || threshold <= 0.0 || threshold >= 1.0)
    {
      cerr << "The threshold value must be a valid distance\n";
      usage(2);
    }

    if (verbose)
      cerr << "Threshold set to " << threshold << endl;
  }

  if (cl.option_present('s'))
  {
    if (! cl.value('s', pool_size) || pool_size < 2)
    {
      cerr << "The pool size specification (-s) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Problem sized for " << pool_size << " fingerprints\n";
  }

  if (cl.option_present('n'))
  {
    if (! cl.value('n', clusters_to_form) || clusters_to_form < 1)
    {
      cerr << "The number of clusters to form (-n) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose)
      cerr << "Will stop processing after forming " << clusters_to_form << " clusters\n";
  }

  if (cl.option_present('r'))
  {
    if (! report_progress.initialise(cl, 'r', verbose))
    {
      cerr << "Cannot initialise report progress option (-r)\n";
      usage(2);
    }
  }

  int write_ascii_protos = 0;
  int write_binary_protos = 0;
  IWString proto_fname;

  if (cl.option_present('Y')) {
    IWString y;
    for (int i = 0; cl.value('Y', y, i); ++i) {
      if (y == "ascii") {
        write_ascii_protos = 1;
      } else if (y == "bin") {
        write_binary_protos = 1;
      } else {
        proto_fname = y;
      }
    }
  }

  if (0 == cl.number_elements())
  {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.number_elements() > 1)
  {
    cerr << "Sorry, cannot handle multiple input files, concatenate them and try again\n";
    usage(2);
  }

  const char * fname = cl[0];

  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << "Cannot open fingerprint file '" << fname << "'\n";
    return 2;
  }

  if (0 == pool_size)
  {
    pool_size = input.count_records_starting_with(identifier_tag);

    if (0 == pool_size)
    {
      cerr << "No occurrences of " << identifier_tag << "' in input\n";
      return 0;
    }

    if (! allocate_pool(pool_size))
      return 3;

    cerr << "Job automatically sized to " << pool_size << " fingerprints\n";
  }
  else
  {
    if (! allocate_pool(pool_size))
      return 3;
  }

  if (! read_pool(input))
  {
    cerr << "Cannot read fingerprints\n";
    return 0;
  }

  if (cl.option_present('h'))
  {
    int h;
    if (! cl.value('h', h) || h < 0)
    {
      cerr << "The maximum number of threads to use (-h) must be a valid whole +ve number\n";
      usage(2);
    }

    omp_set_num_threads(h);
  }

  if (cl.option_present('A'))
  {
    const char * fname = cl.option_value('A');

    if (! do_previously_selected(fname, threshold))
    {
      cerr << "Cannot process previously selected file '" << fname << "'\n";
      return 0;
    }
  }

  int clusters_formed;

  if (cl.option_present('p')) {
    int p;
    if (! cl.value('p', p) || p <= 0 || p > 2) {
      cerr << "The only valid values for the -p option are 1 and 2\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Will process " << p << " leaders at each step\n";
    }

    int * leaders = new int[p]; std::unique_ptr<int[]> free_leaders(leaders);
    float * cand_distances = new float[p]; std::unique_ptr<float[]> free_cand_distances(cand_distances);

    clusters_formed = leader(leaders, p, cand_distances);
  } else {
    clusters_formed = leader();
  }

  if (verbose) {
    cerr << "Clustered " << pool_size << " fingerprints into " << clusters_formed << " clusters\n";
  }

  std::unique_ptr<IWString_and_File_Descriptor> ascii_proto_output;
  std::unique_ptr<TFDataWriter> binary_proto_output;
  if (write_ascii_protos) {
    ascii_proto_output = std::make_unique<IWString_and_File_Descriptor>();
    if (! ascii_proto_output->open(proto_fname)) {
      cerr << "Cannot open file for ascii proto output " << proto_fname << "\n";
      return 1;
    }
  } else if (write_binary_protos) {
    binary_proto_output = std::make_unique<TFDataWriter>();
    if (! binary_proto_output->Open(proto_fname)) {
      cerr << "Cannot open file for binary proto output " << proto_fname << "\n";
      return 1;
    }
  }

  extending_resizable_array<int> cluster_sizes;

  int start = 0;
  for (int c = 1; c <= clusters_formed; c++) {
    int items_in_cluster = 0;

    std::unique_ptr<NearNeighbors::Neighbors> proto;
    if (write_ascii_protos || write_binary_protos) {
      proto.reset(new NearNeighbors::Neighbors());
    }

    for (int j = start; j < pool_size; j++) {
      if (selected[j] != c) {
        continue;
      }

      if (0  == items_in_cluster) {
        start = j;
      }

      items_in_cluster++;
    }

    ++cluster_sizes[items_in_cluster];

    if (proto) {
      proto->set_id(leader_item[start].id().AsString());
      proto->set_smi(leader_item[start].smiles().AsString());
      proto->set_cnum(c - 1);
    } else {
      cout << smiles_tag << leader_item[start].smiles() << ">\n";
      cout << identifier_tag << leader_item[start].id() << ">\n";
      cout << "CLUSTER<" << (c-1) << ">\n";
      cout << "CSIZE<" << items_in_cluster << ">\n";
    }

    start++;

    if (start == pool_size) {
      if (! proto) {
        cout << "|\n";
      }
      break;
    }

    for (int j = start; j < pool_size; j++) {
      if (c != selected[j]) {
        continue;
      }

      if (proto) {
        NearNeighbors::Neighbor * n = proto->add_neighbor();
        n->set_id(leader_item[j].id().AsString());
        n->set_smi(leader_item[j].smiles().AsString());
        n->set_d(distances[j]);
      } else {
        cout << smiles_tag << leader_item[j].smiles() << ">\n";
        cout << identifier_tag << leader_item[j].id() << ">\n";
        cout << distance_tag << distances[j] << ">\n";
      }
    }

    if (ascii_proto_output) {
      WriteTextProto(*proto, *ascii_proto_output);
    } else if (binary_proto_output) {
        WriteBinaryProto(*proto, *binary_proto_output);
    } else {
      cout << "|\n";
    }
  }

  if (verbose) {
    for (uint32_t i = 0; i < cluster_sizes.size(); ++i) {
      if (cluster_sizes[i]) {
        cerr << cluster_sizes[i] << " clusters of size " << i << '\n';
      }
    }
  }

  return 0;
}

int
main(int argc, char ** argv)
{
  prog_name = argv[0];

  int rc = gfp_leader_standard(argc, argv);

  return rc;
}
