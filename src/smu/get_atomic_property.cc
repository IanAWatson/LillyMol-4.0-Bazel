// Do a substructure search for an atomic property.
// Scan molecules, convert to Molecule, do a substructure
// search, and report the atomic property of matched atoms.
#include <iostream>
#include <memory>

#include <google/protobuf/descriptor.h>

#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

#include "smu/dataset.pb.h"
#include "smu/support.h"

namespace get_atomic_property {

using std::cerr;

using GoogleSmu::BondTopology;

struct Options {
  int verbose = 0;

  int conformers_read = 0;

  int no_molecule = 0;

  int no_bond_topologies = 0;

  int no_starting_bond_topology = 0;

  int molecules_searched = 0;

  resizable_array_p<Substructure_Query> queries;

  // For each query, accumulate the atomic properties found.
  Accumulator<double> *acc;

  // And the raw values.
  resizable_array<double> *raw_values;

  std::string property_name;

  // If writing distributions, a file name stem for those
  // output files.
  IWString raw_data_stem;

  Report_Progress report_progress;

  extending_resizable_array<int> hits_to_query;

  // We can ignore N-F and O-F bonds if we wish.
  // For two atomic numbers, set the value MAX(z1,z1) + z2
  // in this array.
  extending_resizable_array<int> ignored;
  // We need a fast means of checking whether or not an
  // atomic number is set in the `ignored` array.
  extending_resizable_array<int> need_to_check;

  // We can specify a min and max value for the atomic property
  // and if the value is outside this range, it will be written.
  double minval = - std::numeric_limits<double>::max();
  double maxval = std::numeric_limits<double>::max();

  // If set, write a molecule that has isotopic values for every
  // atom
  int isotope_on_each_atom = 0;

  // In order to convert floating point values to integer isotopes, we need a multiplier.
  float isotope_multiplier = 1.0f;

  IWString_and_File_Descriptor stream_for_extrema;

  int extreme_values_written = 0;

  char output_separator = ' ';

  private:
    int BuildIgnoreDirective(const const_IWSubstring& s);

    int WriteRawData(const resizable_array<double>& raw_values, 
                     const IWString& column_name, IWString& fname) const;
    int WriteRawData(const resizable_array<double>& raw_values,
                     IWString_and_File_Descriptor& output) const;

  public:
    Options();
    ~Options();

    int Build(Command_Line& cl);

    // Given two atomic numbers, if those atoms are bonded,
    // are these atoms ignored. O-F, N-F typically.
    int Ignore(int z1, int z2);

    // A result has been recorded for a given atom. `ndx` is the 
    // query number providing the match.
    int Extra(const Options& options, Molecule& m, atom_number_t zatom, int ndx, double value);

    // If `value` is outside the range of [minval,maxval] write it
    // to stream_for_extrema.
    int WriteExtremeValue(Molecule& m, int query_number, double value);

    int Report(std::ostream& output) const;
};

Options::Options() {
  acc = nullptr;
  raw_values = nullptr;
}

Options::~Options() {
  if (acc != nullptr) {
    delete [] acc;
  }
  if (raw_values != nullptr) {
    delete [] raw_values;
  }
}

// A simple hash function for a pair of atomic numbers.
// Only valid if both are less than 10.
int
AtomPairHash(int z1, int z2) {
  if (z1 < z2) {
    return 10 * z2 + z1;
  }
  return 10 * z1 + z2;
}

int
Options::Ignore(int z1, int z2) {
  return ignored[AtomPairHash(z1, z2)];
}

// Looks like 'at1-at2'.
int
Options::BuildIgnoreDirective(const const_IWSubstring& s) {
  const_IWSubstring s1, s2;
  if (! s.split(s1, '-', s2) || s1.empty() || s2.empty()) {
    cerr << "Options::BuildIgnoreDirective:invalid form '" << s << "'\n";
    return 0;
  }
  int z1, z2;
  if (! s1.numeric_value(z1) || z1 <= 0 ||
      ! s2.numeric_value(z2) || z2 <= 0) {
    cerr << "Options::BuildIgnoreDirective:invald atomic numbers '" << s << "'\n";
    return 0;
  }

  ignored[AtomPairHash(z1, z2)] = 1;
  need_to_check[z1] = 1;
  need_to_check[z2] = 1;

  return 1;
}

int
Options::Build(Command_Line& cl) {
  verbose = cl.option_count('v');

  if (cl.option_present('q')) {
    if (! process_queries(cl, queries, verbose, 'q')) {
      cerr << "Cannot discern queries\n";
      return 0;
    }
  }

  if (! cl.option_present('p')) {
    cerr << "Must specify the property to fetch via the -p option\n";
    return 0;
  }

  if (cl.option_present('p')) {
    IWString p = cl.string_value('p');
    property_name = p.null_terminated_chars();
    if (verbose) {
      cerr << "Loopking for property " << property_name << '\n';
    }
  }

  if (cl.option_present('m')) {
    if (! cl.value('m', isotope_multiplier) || isotope_multiplier <= 0.0f) {
      cerr << "The isotope multiplier must be a +ve value\n";
      return 0;
    }

    if (verbose) {
      cerr << "Will multiply float value by " << isotope_multiplier << " for conversion to isotopes\n";
    }
  }

  if (cl.option_present('r')) {
    if (! report_progress.initialise(cl, 'r', verbose)) {
      cerr << "Cannot initialise progress reporting\n";
      return 0;
    }
  }

  if (cl.option_present('S')) {
    cl.value('S', raw_data_stem);
    if (verbose) {
      cerr << "Distributions written to files '" << raw_data_stem << "'\n";
    }
  }

  if (cl.option_present('z')) {
    if (! cl.value('z', minval)) {
      cerr << "options::Build:invalid -z qualifier\n";
      return 0;
    }
    if (verbose) {
      cerr << "Will write atoms that are below " << minval << '\n';
    }

    if (! cl.option_present('X')) {
      cerr << "Must also specify stream for extreme values (-X)\n";
      return 0;
    }
  }

  if (cl.option_present('e')) {
    isotope_on_each_atom = 1;
    if (verbose) {
      cerr << "Will write an isotopically labelled molecule with all values\n";
    }
  }

  if (cl.option_present('Z')) {
    if (! cl.value('Z', maxval)) {
      cerr << "options::Build:invalid -Z qualifier\n";
      return 0;
    }
    if (verbose) {
      cerr << "Will write atoms that are above " << maxval << '\n';
    }
    if (! cl.option_present('X')) {
      cerr << "Must also specify stream for extreme values (-X)\n";
      return 0;
    }
  }

  if (cl.option_present('X')) {
    IWString fname = cl.string_value('X');
    if (! stream_for_extrema.open(fname.null_terminated_chars())) {
      cerr << "Options::Build::cannot open -X file '" << fname << "'\n";
      return 0;
    }
    if (verbose) {
      cerr << "Extreme values written to " << fname << "\n";
    }
  }

  if (cl.option_present('I')) {
    const_IWSubstring token;
    for (int i = 0; cl.value('I', token, i); ++i) {
      if (! BuildIgnoreDirective(token)) {
        cerr << "Invalid ignore directive '" << token << "'\n";
        return 0;
      }
    }
  }

  if (cl.option_present('s')) {
    IWString s;
    for (int i = 0; cl.value('s', s, i); ++i) {
      std::unique_ptr<Substructure_Query> q = std::make_unique<Substructure_Query>();
      if (! q->create_from_smarts(s)) {
        cerr << "Invalid smarts '" << s << "'\n";
        return 0;
      }
      queries << q.release();
    }
  }

  if (queries.empty()) {
    cerr << "Options::Build:no queries\n";
    return 0;
  }

  const int nq = queries.number_elements();

  acc = new Accumulator<double>[nq];
  raw_values = new resizable_array<double>[nq];
  for (int i = 0; i < nq; ++i) {
    raw_values[i].reserve(4000000);
  }

  return queries.number_elements();
}

// Write quantile information about `raw_values` to `output`.
// `name` is the name of the feature.
int
ReportDistribution(resizable_array<double>& raw_values,
                   const IWString& name,
                   std::ostream& output) {
  if (raw_values.empty()) {
    cerr << "No values for '" << name << "'\n";
    return 0;
  }

  raw_values.iwqsort_lambda([](const double d1, double d2) {
    if (d1 < d2) {
      return -1;
    } else if (d1 > d2) {
      return 1;
    } else {
      return 0;
    }
  });

  const unsigned int n = raw_values.size();
  int q1 = n / 100;
  int q5 = 5 * n / 100;
  int q10 = 10 * n / 100;
  int q50 = 50 * n / 100;
  int q90 = 9 * n / 10;
  int q95 = 95 * n / 100;
  int q99 = 99 * n / 100;
  constexpr char kSep = ' ';

  output << "min" << kSep << raw_values[0] << kSep <<
            "q1" << kSep << raw_values[q1] << kSep << 
            "q5" << kSep << raw_values[q5] << kSep << 
            "q10" << kSep << raw_values[q10] << kSep << 
            "q50" << kSep << raw_values[q50] << kSep << 
            "q90" << kSep << raw_values[q90] << kSep << 
            "q95" << kSep << raw_values[q95] << kSep << 
            "q99" << kSep << raw_values[q99] << kSep << 
            "max" << kSep << raw_values.back() << '\n';
  output << name << kSep <<
            raw_values[0] << kSep <<
            raw_values[q1] << kSep <<
            raw_values[q5] << kSep <<
            raw_values[q10] << kSep <<
            raw_values[q50] << kSep <<
            raw_values[q90] << kSep <<
            raw_values[q95] << kSep <<
            raw_values[q99] << kSep <<
            raw_values.back() << '\n';
  return 1;
}

int
Options::Extra(const Options& options,
               Molecule& m,
               atom_number_t zatom,
               int ndx, double value) {
  acc[ndx].extra(value);
  raw_values[ndx] << value;

  if (value > minval && value < maxval) {
    return 1;
  }

  if (value < 0.0) {
    m.set_isotope(zatom, static_cast<int>(-value * options.isotope_multiplier));
  } else {
    m.set_isotope(zatom, static_cast<int>(value * options.isotope_multiplier));
  }
  
  return 1;
}

int
Options::WriteExtremeValue(Molecule& m,
                           int query_number,
                           double value) {
  constexpr char kSep = ' ';
  if (! stream_for_extrema.is_open()) {
    return 0;
  }

  ++extreme_values_written;

  stream_for_extrema << m.smiles() << kSep << m.name()
                     << kSep << static_cast<float>(value)
                     << kSep << query_number << '\n';
  stream_for_extrema.write_if_buffer_holds_more_than(8192);
  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << conformers_read << " conformers\n";
  output << no_bond_topologies << " no bond topologies\n";
  output << no_starting_bond_topology << " no starting bond topology\n";
  output << molecules_searched << " molecules searched\n";
  output << extreme_values_written << " extreme_values_written\n";
  output << queries.number_elements() << " queries\n";
  for (int i = 0; i < queries.number_elements(); ++i) {
    output << i << ' ' << queries[i]->comment() << ' ' << acc[i].n() << " matches " <<
              '[' << acc[i].minval() << ',' << acc[i].maxval() << "] mean " << acc[i].average() << '\n';
    ReportDistribution(raw_values[i], queries[i]->comment(), output);

    if (raw_data_stem.length()) {
      IWString fname;
      fname << raw_data_stem;
      IWString column_name;
      if (queries[i]->comment().empty()) {
        fname << '_' << i;
        column_name = queries[i]->comment();
      } else {
        fname << queries[i]->comment();
        column_name = "Value";
      }
      fname << ".txt";
      WriteRawData(raw_values[i], column_name, fname);
    }
  }

  return output.good();
}

int
Options::WriteRawData(const resizable_array<double>& raw_values,
                      const IWString& column_name,
                      IWString& fname) const {
  IWString_and_File_Descriptor output;
  if (! output.open(fname.null_terminated_chars())) {
    cerr << "Options::WriteRawData:cannot open " << fname << "\n";
    return 0;
  }

  output << column_name << output_separator << "count\n";

  return WriteRawData(raw_values, output);
}

// Report the number of different values in `values`
// which is assumed to have been sorted.
template <typename T>
int
NumberDistinctValues(const resizable_array<T>& values) {
  int result = 1;
  const int n = values.number_elements();
  for (int i = 1; i < n; ++i) {
    if (values[i] != values[i - 1]) {
      ++result;
    }
  }

  return result;
}

// Write `data` as value,count pairs to `output`.
// `data` is assumed to be sorted.
template <typename T>
int
WriteValueCounts(const resizable_array<T>& data,
                 char sep,
                 IWString_and_File_Descriptor& output) {
  T current_value = data[0];
  int current_count = 1;
  const int n = data.number_elements();
  for (int i = 1; i < n; ++i) {
    if (data[i] == current_value) {
      ++current_count;
    } else {
      output << current_value << sep << current_count << '\n';
      current_value = data[i];
      current_count = 1;
      output.write_if_buffer_holds_more_than(8196);
    }
  }

  if (current_count > 0) {
    output << current_value << sep << current_count << '\n';
  }

  return output.good();
}

// Bucketise `data` and write out value,count pairs for those buckets.
template <typename T>
int
WriteDiscretizedCounts(const resizable_array<T>& data,
                       const char sep,
                       IWString_and_File_Descriptor& output) {
  constexpr int kNbuckets = 20000;  // arbitrary choice
  int * count = new_int(kNbuckets); std::unique_ptr<int[]> free_count(count);
  const int n = data.number_elements();
  T minval = data[0];
  T maxval = data.back();
  //cerr << "Discretizing " << n << " values btw " << minval << " and " << maxval << '\n';
  for (int i = 0; i < n; ++i) {
    int ndx = static_cast<int>((data[i] - minval) / (maxval - minval) * kNbuckets);
    ++count[ndx];
  }

  const T dx = (maxval - minval) / static_cast<T>(kNbuckets);

  for (int i = 0; i < kNbuckets; ++i) {
    output << (minval + i * dx) << sep << count[i] << '\n';
    output.write_if_buffer_holds_more_than(8196);
  }

  return output.good();
}

// Write value,count form of `raw_values` to `output`.
// Decide if there are a small number of discrete values or
// discretise it ourselves.
int
Options::WriteRawData(const resizable_array<double>& raw_values,
                      IWString_and_File_Descriptor& output) const {
  int distinct = NumberDistinctValues(raw_values);
  if (verbose) {
    cerr << "Data as " << distinct << " distinct values\n";
  }
  if (distinct < 1000) {
    return WriteValueCounts(raw_values, output_separator, output);
  }

  return WriteDiscretizedCounts(raw_values, output_separator, output);
}


void
Usage(int rc) {
  cerr << "Scans Conformer protos for atoms matching a substructure query\n";
  cerr << " -q <query>    query specification\n";
  cerr << " -s <smarts>   query specification as smarts\n";
  cerr << " -I <a-b>      ignore atoms involved in a-b bonds\n";
  cerr << " -S <stem>     write raw data for each query to file(s) starting with <stem>\n";
  cerr << " -x <minval>   define the min value for reporting extrema (-Z)\n";
  cerr << " -X <maxval>   define the max value for reporting extrema (-Z)\n";
  cerr << " -e            write a molecule with each atom isotopically labelled with the value to -Z file\n";
  cerr << " -Z <fname>    write molecules that have a value outside [min,max] to <fname>\n";
  cerr << " -p <property> name of a AtomicMolecularProperty field in the proto\n";
  cerr << " -r <rpt>      report progress every <rpt> items processed\n";
  cerr << " -v            verbose output\n";
  ::exit(rc);
}

// Maybe return an atomic value from `conformer` for `atom_number`.
std::optional<double>
GetValue(const GoogleSmu::Molecule& conformer, int atom_number) {
  if (! conformer.properties().has_nmr_isotropic_shielding_pbe0_aug_pcs_1()) {
    return std::nullopt;
  }
  const auto& nmr = conformer.properties().nmr_isotropic_shielding_pbe0_aug_pcs_1().values();
  return nmr[atom_number];
}

int foo(int s) {
  return s + 1;
}

// First we need to find the AtomicMolecularProperty message in proto.properties().
// Once that is found, we need to get the message, and then get the Reflection
// for that message, and the field descriptor for the 'values' field.
std::optional<double>
GetValue(const GoogleSmu::Molecule& proto,
         int atom_number,
         const std::string& property_name) {
  const /*Descriptor*/auto* descriptor = proto.properties().GetDescriptor();
  const /*Reflection*/auto* reflection = proto.properties().GetReflection();
  //const FieldDescriptor * field_descriptor = proto.FindFieldByName(property);
  const auto * field_descriptor = descriptor->FindFieldByName(property_name);
  if (field_descriptor == nullptr) {
    cerr << "No feature '" << property_name << "'\n";
    return std::nullopt;
  }

  // cerr << reflection->SpaceUsedLong(proto.properties()) << " space used\n";
  const auto& values = reflection->GetMessage(proto.properties(), field_descriptor);
  auto* r2 = values.GetReflection();
  auto* d2 = values.GetDescriptor();
  auto* f2 = d2->FindFieldByName("values");
  return r2->GetRepeatedDouble(values, f2, atom_number);
}

int
PlaceIsotopeOnEachAtom(Options& options,
                  const GoogleSmu::Molecule& conformer,
                  Molecule& m,
                  IWString_and_File_Descriptor& output) {
  constexpr char kSep = ' ';
  const int matoms = m.natoms();
  resizable_array<float> values;
  values.reserve(matoms);
  for (int i = 0; i < matoms; ++i) {
    if (m.atomic_number(i) == 1) {
      continue;
    }

    std::optional<double> maybe_value = GetValue(conformer, i, options.property_name);
    if (! maybe_value) {
      continue;
    }
    const double value = *maybe_value;
    values << static_cast<float>(value);
    if (value < 0) {
      m.set_isotope(i, static_cast<int>(-value * options.isotope_multiplier));
    } else {
      m.set_isotope(i, static_cast<int>(value * options.isotope_multiplier));
    }
  }

  output << m.smiles() <<
         kSep << m.name();
  for (float v : values) {
    output << kSep << v;
  }
  output << '\n';
  output.write_if_buffer_holds_more_than(8192);

  return 1;
}

// Do substructure searches over `m`.
int
GetAtomicProperty(Options& options,
                  const GoogleSmu::Molecule& conformer,
                  Molecule& m,
                  IWString_and_File_Descriptor& output) {
  Molecule_to_Match target(&m);

  options.molecules_searched++;

  int is_outlier = 0;

  const int nq = options.queries.number_elements();
  for (int qnum = 0; qnum < nq; ++qnum) {
    Substructure_Results sresults;
    const int nhits = options.queries[qnum]->substructure_search(target, sresults);
    if (nhits == 0) {
      continue;
    }
    options.hits_to_query[qnum]++;
    m.transform_to_non_isotopic_form();
    std::optional<double> value_to_write;

    for (int j = 0; j < nhits; ++j) {
      const Set_of_Atoms * e = sresults.embedding(j);
      std::optional<double> maybe_value = GetValue(conformer, e->item(0), options.property_name);
      if (! maybe_value) {
        continue;
      }
      const double value = *maybe_value;
      if (value_to_write == std::nullopt) {
        value_to_write = value;
      }
      options.Extra(options, m, e->item(0), qnum, value);
      if (value < options.minval || value > options.maxval) {
        ++is_outlier;
      }
    }

    if (is_outlier) {
      options.WriteExtremeValue(m, qnum, *value_to_write);
    }
  }

  return 1;
}

// Set the molecule name to `molecule_id`.
void
SetName(int molecule_id,
        Molecule& m) {
  IWString name;
  name << molecule_id;
  m.set_name(name);
}

int
GetAtomicProperty(Options& options,
                  const GoogleSmu::Molecule& conformer,
                  IWString_and_File_Descriptor& output) {
  if (conformer.bond_topologies().size() == 0) {
    options.no_bond_topologies++;
    return 1;
  }

  if (options.report_progress()) {
    cerr << "Processed " << options.conformers_read << " molecules\n";
  }

  if (conformer.properties().errors().fate() != GoogleSmu::Properties::FATE_SUCCESS) {
    return 1;
  }

  int startingbt = smu::IndexOfStartingBTid(conformer.bond_topologies());
  if (startingbt < 0) {
    options.no_starting_bond_topology++;
    startingbt = 0;
  }

  const BondTopology& bt_start = conformer.bond_topologies(startingbt);

  std::optional<Molecule> maybe_mol = smu::MoleculeFromBondTopology(bt_start);
  if (! maybe_mol) {
    options.no_molecule++;
    return 1;
  }
  SetName(conformer.molecule_id(), *maybe_mol);

  if (options.isotope_on_each_atom) {
    return PlaceIsotopeOnEachAtom(options, conformer, *maybe_mol, output);
  }

  return GetAtomicProperty(options, conformer, *maybe_mol, output);
}

int
GetAtomicProperty2(Options& options,
                   const const_IWSubstring& buffer,
                   IWString_and_File_Descriptor& output) {
  const std::string as_string(buffer.data(), buffer.length());
  GoogleSmu::Molecule conformer;
  if (! conformer.ParseFromString(as_string)) {
    cerr << "Cannot decode proto\n";
    return 0;
  }
  options.conformers_read++;

  return GetAtomicProperty(options, conformer, output);
}

int
GetAtomicProperty(Options& options,
                  iw_tf_data_record::TFDataReader& reader,
                  IWString_and_File_Descriptor& output) {
  while (reader.good() && ! reader.eof()) {
    std::optional<const_IWSubstring> data = reader.Next();
    if (! data) {
      return 1;
    }
    if (! GetAtomicProperty2(options, *data, output)) {
      return 0;
    }
  }

  return 1;
}

int
GetAtomicProperty(Options& options,
                  const char * fname,
                  IWString_and_File_Descriptor& output) {
  iw_tf_data_record::TFDataReader reader(fname);
  if (! reader.good()) {
    cerr << "GetAtomicProperty:cannto open '" << fname << "'\n";
    return 0;
  }

  return GetAtomicProperty(options, reader, output);
}

int
GetAtomicProperty(int argc, char** argv) {
  Command_Line cl(argc, argv, "vq:s:A:r:S:I:z:Z:X:ep:m:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  Options options;
  if (! options.Build(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  for (const char * fname: cl) {
    if (! GetAtomicProperty(options, fname, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  options.Report(cerr);

  return 0;
}

}  // namespace get_atomic_property

int
main(int argc, char ** argv)
{
  int rc = get_atomic_property::GetAtomicProperty(argc, argv);

  return rc;
}
