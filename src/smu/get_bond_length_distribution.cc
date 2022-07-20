// Extracts bond length distributions from 3D files

#include <stdlib.h>

#include <iostream>
#include <map>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/proto_support.h"

#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

#include "smu/bond_length_distribution.pb.h"

namespace get_bond_length_distribution {

using std::cerr;

// Atoms cannot be this close.
constexpr float kVeryShort = 0.8;

void
Usage(int rc) {
  cerr << "Extracts bond length distributions from molecules\n";
  cerr << " -S <stem>      name stem for bond length distribution data\n";
  cerr << " -q <query>     query(s) identifying bonds to be measured\n";
  cerr << " -b             matched atoms in a query match must be bonded\n";
  cerr << " -v             verbose output\n";

  ::exit(rc);
}

struct BucketisedDistances
{
  // From bucketised distance to count.
  std::map<int, int> bond_length_data;

  void Extra(float distance);

  int Write(const char sep, IWString_and_File_Descriptor& output) const;
  int WriteProto(Smu::BondLengthDistribution& proto, IWString& fname) const;
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
  Accumulator_Int<int> acc;
  for (auto [bucket, count] : bond_length_data) {
    output << bucket / 1000.0f << sep << count << '\n';
    acc.extra(bucket, count);
  }

  cerr << acc.n() << " values btw " << acc.minval() << " and " << acc.maxval() << " mean " << acc.average() << '\n';

  return 1;
}

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

struct Options
{
  int verbose = 0;

  int molecules_read = 0;

  int no_3d_data = 0;

  // The queries that define the bonds to be measured.
  // Note that we do not check whether or not the atoms are
  // bonded, so that increases the generalizability of this tool.
  resizable_array_p<Substructure_Query> queries;

  // For each query, bond length data.

  std::unique_ptr<BucketisedDistances[]> distances;

  // How many queries hit each input molecule.
  extending_resizable_array<int> queries_hit;

  bool matched_atoms_must_be_bonded = false;

  int non_bonded_pairs_skipped = 0;

  // Output separator.
  char sep = ' ';

  private:
    int AccumulateDistances(const Molecule& m, const Substructure_Results& sresults, BucketisedDistances& destination);

    int Write(int ndx, IWString_and_File_Descriptor& output) const;

  public:

    int Build(Command_Line& cl);

    int Report(std::ostream& output) const;

    int Process(Molecule& m);

    int Write(const IWString& name_stem);
    int WriteProto(const IWString& name_stem) const;
};

int
Options::Build(Command_Line& cl)
{
  verbose = cl.option_count('v');

  if (! cl.option_present('q')) {
    cerr << "Must specify one or more queries via the -q option\n";
    Usage(1);
  }

  if (! process_queries(cl, queries, verbose, 'q')) {
    cerr << "Cannot read queries\n";
    return 0;
  }

  matched_atoms_must_be_bonded = cl.option_present('b');

  distances.reset(new BucketisedDistances[queries.number_elements()]);

  return 1;
}

int
Options::Report(std::ostream& output) const
{
  output << "Read " << molecules_read << " molecules\n";
  output << no_3d_data << " had no 3d data\n";
  if (matched_atoms_must_be_bonded) {
    output << non_bonded_pairs_skipped << " non bonded matched atoms skipped\n";
  }
  for (int i = 0; i < queries_hit.number_elements(); ++i) {
    if (queries_hit[i] > 0) {
      cerr << queries_hit[i] << " molecules hit " << i << " queries\n";
    }
  }

  return output.good();
}

int
Options::Process(Molecule& m)
{
  if (m.highest_coordinate_dimensionality() < 3) {
    ++no_3d_data;
    cerr << "Options::Process:no 3D info in " << m.name() << '\n';
    return 1;
  }

  Molecule_to_Match target(&m);

  int current_molecule_queries_hit = 0;
  for (int i = 0; i < queries.number_elements(); ++i) {
    Substructure_Results sresults;
    const int nhits = queries[i]->substructure_search(target, sresults);
    if (nhits == 0) {
      continue;
    }
    ++current_molecule_queries_hit;
    AccumulateDistances(m, sresults, distances[i]);
  }

  queries_hit[current_molecule_queries_hit]++;
  if (current_molecule_queries_hit == 0) {
    cerr << "No queries hit " << m.smiles() << ' ' << m.name() << '\n';
  }

  return 1;
}

int
Options::AccumulateDistances(const Molecule& m,
                             const Substructure_Results& sresults,
                             BucketisedDistances& destination)
{
  for (const Set_of_Atoms * e : sresults.embeddings()) {
    if (e->number_elements() < 2) {
      cerr << "Options::AccumulateDistances:only " << e->number_elements() << " matched atoms\n";
      return 0;
    }
    const atom_number_t a1 = e->item(0);
    const atom_number_t a2 = e->item(1);

    if (matched_atoms_must_be_bonded && ! m.are_bonded(a1, a2)) {
      ++non_bonded_pairs_skipped;
      continue;
    }

    const float d = m.distance_between_atoms(a1, a2);
    if (d < kVeryShort) {
      cerr << "Options::AccumulateDistances:distance too short " << d << '\n';
      return 0;
    }
    destination.Extra(d);
  }

  return 1;
}

int
Options::Write(const IWString& name_stem) {
  for (int i = 0; i < queries.number_elements(); ++i) {
    IWString fname;
    fname << name_stem << queries[i]->comment() << ".txt";
    IWString_and_File_Descriptor output;
    if (! output.open(fname.null_terminated_chars())) {
      cerr << "Options::Write:cannot open '" << fname << "'\n";
      return 0;
    }
    if (! Write(i, output)) {
      cerr << "Options::Write:error processing '" << fname << "'\n";
      return 0;
    }
  }

  return 1;
}

// Write the distances associated with query `ndx` to `output`.
int
Options::Write(int ndx, IWString_and_File_Descriptor& output) const
{
  output << queries[ndx]->comment() << sep << "count\n";
  return distances[ndx].Write(sep, output);
}

int
CharToAtomicNumber(char c,
                   int& z)
{
  switch (c) {
    case 'H':
      z = 1;
      return 1;
    case 'C':
      z = 6;
      return 1;
    case 'N':
      z = 7;
      return 1;
    case 'O':
      z = 8;
      return 1;
    case 'F':
      z = 9;
      return 1;
    default:
      cerr << "Unrecognised atomic symbol '" << c << "'\n";
      z = 0;
      return 0;
  }
}

int
CharToBondType(char b,
               SubstructureSearch::BondType& bond)
{
  switch (b) {
    case '1':
      bond = SubstructureSearch::SS_SINGLE_BOND;
      return 1;
    case '2':
      bond = SubstructureSearch::SS_DOUBLE_BOND;
      return 1;
    case '3':
      bond = SubstructureSearch::SS_TRIPLE_BOND;
      return 1;
    case '4':
      bond = SubstructureSearch::SS_AROMATIC_BOND;
      return 1;
    default:
      cerr << "CharToBondType:unrecognised bond type '" << b << "'\n";
      return 0;
  }
}

int
Options::WriteProto(const IWString& name_stem) const {
  for (int i = 0; i < queries.number_elements(); ++i) {
    const IWString& name = queries[i]->comment();
    if (name.length() != 3) {
      cerr << "Options::WriteProto:invalid query name '" << name << "'\n";
      return 0;
    }
    int z1, z2;
    SubstructureSearch::BondType btype;
    if (! CharToAtomicNumber(name[0], z1) ||
        ! CharToBondType(name[1], btype) ||
        ! CharToAtomicNumber(name[2], z2)) {
      cerr << "Cannot parse query name '" << name << "'\n";
      return 0;
    }

    IWString fname;
    fname << name_stem << queries[i]->comment() << ".textproto";

    Smu::BondLengthDistribution proto;
    proto.set_at1(z1);
    proto.set_at2(z2);
    proto.set_btype(btype);
    if (! distances[i].WriteProto(proto, fname)) {
      cerr << "Options::WriteProto:cannot write '" << fname << "'\n";
      return 0;
    }
  }

  return 1;
}

void
Preprocess(Molecule& m) {
  return;
}

int
GetBondLengthDistribution(Options& options,
                          data_source_and_type<Molecule>& input,
                          IWString_and_File_Descriptor& output)
{
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    options.molecules_read++;

    Preprocess(*m);

    if (! options.Process(*m)) {
      return 0;
    }
  }

  return 1;
}

int
GetBondLengthDistribution(Options& options,
                          FileType input_type,
                          const char * fname,
                          IWString_and_File_Descriptor& output)
{
  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "Cannot open " << fname << "'\n";
    return 0;
  }

  return GetBondLengthDistribution(options, input, output);
}

int
GetBondLengthDistribution(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vE:A:i:bq:S:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  if (cl.empty()) {
    Usage(1);
  }

  Options options;

  if (! options.Build(cl)) {
    cerr << "Cannot initialise\n";
    return 1;
  }

  FileType input_type = FILE_TYPE_SMI;
  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot establish input type (-i)\n";
      return 1;
    }
  }

  IWString_and_File_Descriptor output(1);

  for (const char * fname: cl) {
    if (! GetBondLengthDistribution(options, input_type, fname, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  options.Report(cerr);

  IWString output_name_stem;
  if (cl.option_present('S')) {
    cl.value('S', output_name_stem);
  } else {
    output_name_stem = "bond_length_distribution_";
  }

  options.Write(output_name_stem);
  options.WriteProto(output_name_stem);

  return 0;
}

}

int
main(int argc, char** argv) {
  const int rc = get_bond_length_distribution::GetBondLengthDistribution(argc, argv);

  return rc;
}
