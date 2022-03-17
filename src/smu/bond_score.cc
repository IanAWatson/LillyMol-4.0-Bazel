// Given a mean bond length for each bond type, compute
// a score of how different the bonds in a molecule are.

#include <stdlib.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <unordered_map>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwmisc/proto_support.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

#include "smu/bond_length_distribution.pb.h"
#include "smu/dataset.pb.h"

#include "support.h"

namespace bond_score {

using std::cerr;

void
Usage(int rc)
{
  ::exit(rc);
}

struct AtomicNumberBtypeAtomicNumber {
  atomic_number_t z1 = 0;
  int btype = 0;
  atomic_number_t z2 = 0;
};

int
BondToNumber(const Bond& b)
{
  if (b.is_aromatic()) {
    return 4;
  }

  if (b.is_single_bond()) {
    return 1;
  }

  if (b.is_double_bond()) {
    return 2;
  }

  if (b.is_triple_bond()) {
    return 3;
  }

  return 0;
}

int
BondToNumber(SubstructureSearch::BondType btype)
{
  switch (btype) {
    case SubstructureSearch::SS_AROMATIC_BOND:
      return 4;
    case SubstructureSearch::SS_SINGLE_BOND:
      return 1;
    case SubstructureSearch::SS_DOUBLE_BOND:
      return 2;
    case SubstructureSearch::SS_TRIPLE_BOND:
      return 3;
    case SubstructureSearch::SS_BOND:  // should never happen.
      cerr << "BondToNumber:invalid bond type\n";
      return 0;
  }

  // Should never come here.
  return 0;
}

uint32_t
HashValue(atomic_number_t z1,
          int btype,
          atomic_number_t z2) {
  if (z1 > z2) {
    std::swap(z1, z2);
  }

  return 1000000 * z1 + 100 * btype + z2;
}

// Hash function for hashing two atomic numbers (assumed < 10) and
// an integer bond type - in the range [1-4].
uint32_t
HashValue(atomic_number_t z1,
          const Bond& bond,
          atomic_number_t z2)
{
  return HashValue(z1, BondToNumber(bond), z2);
}

uint32_t
HashValue(const Smu::BondLengthDistribution& proto)
{
  atomic_number_t z1 = proto.at1();
  atomic_number_t z2 = proto.at2();

  return HashValue(z1, BondToNumber(proto.btype()), z2);
}

struct Options
{
  int verbose = 0;

  int molecules_read = 0;

  // Molecules with no 3d data.
  int no_3d_data = 0;

  int molecules_with_unclassified_bonds = 0;

  // Accumulate stats on the scores assigned.
  Accumulator<double> acc_score;

  // Output separator.
  char sep = ' ';

  // A mapping from a hashed pair of atoms and a distribution.
  std::unordered_map<uint32_t, Smu::BondLengthDistribution> dist;

  // If reading TFDataRecords

  int no_optimised_geometry = 0;

  private:
    int BondScore(Conformer& conformer, iw_tf_data_record::TFDataWriter& output);
    void ComputeDeviationFromReference(const Geometry& geom, BondTopology& bt);

  public:
    Options() {
    }

    int Build(Command_Line& cl);

    int Process(Molecule& m, IWString_and_File_Descriptor& output);

  // Reading dataset protos.
    int Process(iw_tf_data_record::TFDataReader& input, iw_tf_data_record::TFDataWriter& output);

    int Report(std::ostream& output) const;
};

int
Options::Build(Command_Line& cl) {
  verbose = cl.option_count('v');

  int ndist = cl.option_count('D');
  if (ndist == 0) {
    cerr << "Must specify BondLengthDistribution proto via the -D option\n";
    Usage(1);
  }

  IWString d;
  for (int i = 0; cl.value('D', d, i); ++i) {
    std::optional<Smu::BondLengthDistribution> proto = iwmisc::ReadTextProto<Smu::BondLengthDistribution>(d);
    if (! proto) {
      cerr << "Options::Build:cannot process '" << d << "'\n";
      return 0;
    }
    const uint32_t hash = HashValue(*proto);
    dist.emplace(hash, *proto);
  }

  return 1;
}

int
Options::Report(std::ostream& output) const
{
  output << "Read " << molecules_read << " molecules\n";
  output << molecules_with_unclassified_bonds << " molecules with unclassified bonds\n";
  output << acc_score.n() << " scores\n";
  if (acc_score.n() == 0) {
    return 1;
  }
  output << "btw " << acc_score.minval() << " and " << acc_score.maxval() << " mean " << acc_score.average() << '\n';

  return 1;
}

float
ComputeScore(float distance,
             const Smu::BondLengthDistribution& dist)
{
  float score = abs(distance - dist.mean());
  return score;
}

int
Options::Process(Molecule& m,
                 IWString_and_File_Descriptor& output)
{
  if (m.highest_coordinate_dimensionality() < 3) {
    ++no_3d_data;
    cerr << "Options::Process:no 3D info in " << m.name() << '\n';
    return 1;
  }

  float score = 0.0f;
  int unclassified_bonds = 0;
  for (const Bond * b : m.bond_list()) {
    const atom_number_t a1 = b->a1();
    const atom_number_t a2 = b->a2();
    const atomic_number_t z1 = m.atomic_number(a1);
    const atomic_number_t z2 = m.atomic_number(a2);
    const uint32_t h = HashValue(z1, *b, z2);
    const auto iter = dist.find(h);
    if (iter == dist.end()) {
      ++unclassified_bonds;
      continue;
    }
    score += ComputeScore(m.distance_between_atoms(a1, a2), iter->second);
  }

  output << m.smiles() <<
            sep << m.name() <<
            sep << unclassified_bonds <<
            sep << score <<
            sep << (score / (m.nedges() - unclassified_bonds)) <<
            '\n';

  if (unclassified_bonds) {
    ++molecules_with_unclassified_bonds;
  }
  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
Preprocess(Molecule& m)
{
  return 1;
}

int
BondScore(Options& options,
          data_source_and_type<Molecule>& input,
          IWString_and_File_Descriptor& output)
{
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    options.molecules_read++;

    Preprocess(*m);

    if (! options.Process(*m, output)) {
      return 0;
    }
  }

  return 1;
}

int
SmuBtypeToNumber(BondTopology::BondType bt)
{
  switch (bt) {
    case BondTopology::BOND_SINGLE:
      return 1;
    case BondTopology::BOND_DOUBLE:
      return 2;
    case BondTopology::BOND_TRIPLE:
      return 3;
    case BondTopology::BOND_UNDEFINED:
      cerr << "BOND_UNDEFINED encountered\n";
      return 0;
    default:
      cerr << "Unrecognised bond " << bt << '\n';
      return 0;
  }
}

void
Options::ComputeDeviationFromReference(const Geometry& geom,
                                       BondTopology& bt)
{
  float score = 0.0;
  int n = 0;
  int unclassified_bonds = 0;
  for (const auto& bond : bt.bonds()) {
    int a1 = bond.atom_a();
    int a2 = bond.atom_b();
    std::optional<int> at1 = smu::AtomTypeToAtomicNumber(bt.atoms(a1));
    std::optional<int> at2 = smu::AtomTypeToAtomicNumber(bt.atoms(a2));
    if (! at1 || ! at2) {
      cerr << "Cannot discern atomic numbers\n";
      return;
    }
    int btype = SmuBtypeToNumber(bond.bond_type());
    const uint32_t hash = HashValue(*at1, btype, *at2);
    const auto iter = dist.find(hash);
    if (iter == dist.end()) {
      ++unclassified_bonds;
      continue;
    }
    const double d = smu::DistanceBetweenAtoms(geom, a1, a2);
    score += ComputeScore(d, iter->second);
    ++n;
  }
  bt.set_deviation_from_reference(score / n);

  acc_score.extra(score);

  if (unclassified_bonds > 0) {
    ++molecules_with_unclassified_bonds;
  }
}

int
Options::BondScore(Conformer& conformer,
          iw_tf_data_record::TFDataWriter& output)
{
  if (conformer.optimized_geometry().atom_positions().size() == 0) {
    ++no_optimised_geometry;
    return 1;
  }

  for (BondTopology& bt : *conformer.mutable_bond_topologies()) {
    ComputeDeviationFromReference(conformer.optimized_geometry(), bt);
  }

  return output.WriteSerializedProto(conformer);
}


int
Options::Process(iw_tf_data_record::TFDataReader& input,
                 iw_tf_data_record::TFDataWriter& output)
{
  while (true) {
    std::optional<Conformer> data = input.ReadProto<Conformer>();
    if (! data) {
      cerr << "Did not read data\n";
      return 1;
    }

    if (! BondScore(*data, output)) {
      return 0;
    }
  }
}

int
BondScore(Options& options,
          IWString& fname,
          iw_tf_data_record::TFDataWriter& output) 
{
  iw_tf_data_record::TFDataReader input(fname);
  if (! input.good()) {
    cerr << "BondScore:cannot open '" << fname << "'\n";
    return 0;
  }

  cerr << "OPened '" << fname << "' continuing\n";
  return options.Process(input, output);
}

int
BondScore(Options& options,
          FileType input_type,
          const char * fname,
          IWString_and_File_Descriptor& output)
{
  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "Cannot open " << fname << "'\n";
    return 0;
  }

  return BondScore(options, input, output);
}

int
BondScore(int argc, char** argv)
{
  Command_Line cl(argc, argv, "vA:Ei:D:I:O:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  FileType input_type = FILE_TYPE_SMI;
  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot establish input type (-i)\n";
      return 1;
    }
  }

  Options options;

  if (! options.Build(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }


  if (cl.option_present('I')) {
    if (! cl.option_present('O')) {
      cerr << "When reading TFDataRecords, must specify output file via the -O option\n";
      Usage(1);
    }
    IWString input_fname = cl.option_value('I');
    IWString output_fname = cl.option_value('O');
    iw_tf_data_record::TFDataWriter output;
    if (! output.Open(output_fname)) {
      cerr << "Cannot open TFDataRecord output " << output_fname << "'\n";
      return 1;
    }
    if (! BondScore(options, input_fname, output)) {
      cerr << "Failure processing '" << input_fname << "'\n";
      return 1;
    }
  } else {
    if (cl.empty()) {
      Usage(1);
    }

    IWString_and_File_Descriptor output(1);

    for (const char * fname: cl) {
      if (! BondScore(options, input_type, fname, output)) {
        cerr << "Fatal error processing '" << fname << "'\n";
        return 1;
      }
    }
  }

  options.Report(cerr);

  return 0;
}

}  // namespace bond_score

int
main(int argc, char** argv) {
  const int rc = bond_score::BondScore(argc, argv);

  return rc;
}
