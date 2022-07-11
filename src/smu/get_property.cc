// Extract any float value from SMU data.
// A two pass process that profiles the raw values and then identifies extreme values.

#include <cmath>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"

#include "smu/dataset.pb.h"
#include "smu/support.h"

namespace get_nmr {

using std::cerr;

using GoogleSmu::BondTopology;

int
ToAtomicNumber(const BondTopology::AtomType atype) {
  switch (atype) {
    case BondTopology::ATOM_H:
      return 1;
    case BondTopology::ATOM_C:
      return 6;
    case BondTopology::ATOM_N:
      return 7;
    case BondTopology::ATOM_NPOS:
      return 7;
    case BondTopology::ATOM_O:
      return 8;
    case BondTopology::ATOM_ONEG:
      return 8;
    case BondTopology::ATOM_F:
      return 9;
    case BondTopology::ATOM_UNDEFINED:
      cerr << "ToAtomicNumber: undefined atom encountered\n";
      return -1;
    default:
      cerr << "ToAtomicNumber: undefined atom encountered\n";
      return -1;
  }

  cerr << "Unrecognised atom type " << atype << '\n';
  return -1;
}

struct JobOptions {
  int verbose = 0;

  int items_read = 0;

  // The number of Conformers that pass the inclusion tests.
  int results_obtained = 0;

  int nprocess = std::numeric_limits<int>::max();

  // We can ignore N-F and O-F bonds if we wish.
  // For two atomic numbers, set the value 10 * MAX(z1,z1) + z2
  // in this array.
  extending_resizable_array<int> ignored;

  // We need a fast means of checking whether or not an
  // atomic number is set in the `ignored` array.
  extending_resizable_array<int> need_to_check;

  Report_Progress report_progress;
  
  std::unordered_set<int> btid_to_fetch;

  // We convert all values to int, and this value is used for that conversion.
  // For something like nmr shielding, where values are large numbers, a value
  // of 1 is fine. For something like partial charges, a value of 1000 might be
  // needed;
  double scaling_factor = 1.0;

  public:

    int BuildIgnoreDirective(const const_IWSubstring& s);

    int Ignore(int z1, int z2);

    int Ignore(const BondTopology& bt, int a1, int z1);

    void SetScalingFactor(double s) {
      scaling_factor = s;
    }
    double ScalingFactor() {
      return scaling_factor;
    }
};

int
AtomPairHash(int z1, int z2) {
  constexpr int kBase = 10;
  if (z1 < z2) {
    return kBase * z2 + z1;
  }
  return kBase * z1 + z2;
}

int
JobOptions::BuildIgnoreDirective(const const_IWSubstring& s) {
  const_IWSubstring s1, s2;
  if (! s.split(s1, '-', s2) || s1.empty() || s2.empty()) {
    cerr << "JobOptions::BuildIgnoreDirective:invalid form '" << s << "'\n";
    return 0;
  }
  int z1, z2;
  if (! s1.numeric_value(z1) || z1 <= 0 ||
      ! s2.numeric_value(z2) || z2 <= 0 || z2 == z1) {
    cerr << "JobOptions::BuildIgnoreDirective:invald atomic numbers '" << s << "'\n";
    return 0;
  }

  ignored[AtomPairHash(z1, z2)] = 1;
  need_to_check[z1] = 1;
  need_to_check[z2] = 1;

  return 1;
}

int
JobOptions::Ignore(int z1, int z2) {
  return ignored[AtomPairHash(z1, z2)];
}

int
JobOptions::Ignore(const BondTopology& bt, int a1, int z1) {
  if (! need_to_check[z1]) {
    return 0;
  }

  for (const BondTopology::Bond& bond : bt.bonds()) {
    int a = bond.atom_a();
    int b = bond.atom_b();

    if (a == a1) {
      int z2 = ToAtomicNumber(bt.atoms(b));
      if (Ignore(z1, z2)) {
        return 1;
      }
    } else if (b == a1) {
      int z2 = ToAtomicNumber(bt.atoms(a));
      if (Ignore(z1, z2)) {
        return 1;
      }
    }
  }

  return 0;
}

void
Usage(int rc) {
  cerr << " -G <fname>      file containing BondTopology ids to process\n";
  cerr << " -S <fname>      file name stem for distributions\n";
  cerr << " -O <fname>      file name for outlier values\n";
  cerr << " -I <a-b>        ignore atoms involved in a-b bonds\n";
  cerr << " -r <n>          report progress every <n> Conformers\n";
  cerr << " -s <scale>      multiply raw values by <scale> to convert to int\n";
  cerr << " -v              verbose output\n";
  exit(rc);
}

int
WriteIfExtremeValue(JobOptions& options,
                    const GoogleSmu::Molecule& conformer,
                    const std::vector<std::pair<int, int>> & quantiles,
                    IWString_and_File_Descriptor& output) {
  const int startingbt = smu::IndexOfStartingBTid(conformer.bond_topologies());
  const BondTopology& bt0 = conformer.bond_topologies(startingbt);
  const auto& nmr = conformer.properties().nmr_isotropic_shielding_pbe0_aug_pcs_1().values();
  // const auto& nmr = conformer.properties().partial_charges_paboon_pbe0_aug_pc_1().values();  no values
 // const auto& nmr = conformer.properties().partial_charges_esp_fit_pbe0_aug_pc_1().values();
  const int matoms = bt0.atoms().size();
  resizable_array<int> outlier_atoms;
  resizable_array<int> outlier_values;
  for (int i = 0; i < matoms; ++i) {
    const int atomic_number = ToAtomicNumber(bt0.atoms(i));
    if (options.Ignore(bt0, i, atomic_number)) {
      continue;
    }
    const int shielding = static_cast<int>(options.ScalingFactor() * nmr[i]);
    if (shielding < quantiles[atomic_number].first) {
    } else if (shielding > quantiles[atomic_number].second) {
    } else {
      continue;
    }
    outlier_atoms << i;
    outlier_values << shielding;
  }

  // If none of the atoms are outliers, we are not interested.
  if (outlier_atoms.empty()) {
    return 1;
  }
  std::optional<Molecule> maybe_mol = smu::MoleculeFromBondTopology(bt0);
  if (! maybe_mol) {
    return 1;
  }

  // We want to identify the most extreme value in the output.
  // The value and the atomic symbol of the atom leading to that.
  int max_outlier_value = 0;
  IWString outlier_atomic_symbol;
  for (int i = 0; i < outlier_atoms.number_elements(); ++i) {
    int v = outlier_values[i];
    maybe_mol->set_isotope(outlier_atoms[i], std::abs(v));
    if (std::abs(v) > std::abs(max_outlier_value)) {
      max_outlier_value = v;
      outlier_atomic_symbol = maybe_mol->atomic_symbol(outlier_atoms[i]);
    }
  }

  constexpr char sep = ' ';

  output << maybe_mol->smiles() << sep
         << conformer.molecule_id() << sep
         << (conformer.molecule_id() / 1000) << sep
         << max_outlier_value << sep
         << outlier_atomic_symbol << sep
         << conformer.properties().single_point_energy_atomic_b5().value()
         << '\n';

  output.write_if_buffer_holds_more_than(8192);

  return 1;
}

void
CopyData(const GoogleSmu::Molecule& source,
         GoogleSmu::Molecule& destination) {
  destination.set_molecule_id(source.molecule_id());
  destination.mutable_properties()->mutable_errors()->set_fate(source.properties().errors().fate());
  for (const auto & existing_bt : source.bond_topologies()) {
    BondTopology* bt = destination.add_bond_topologies();
    *bt = existing_bt;
  }

  destination.mutable_properties()->mutable_single_point_energy_atomic_b5()->set_value(source.properties().single_point_energy_atomic_b5().value());
  //for (const double existing_nmr : source.properties().partial_charges_esp_fit_pbe0_aug_pc_1().values()) {
  //  destination.mutable_properties()->mutable_partial_charges_esp_fit_pbe0_aug_pc_1()->add_values(existing_nmr);
  //}
  *destination.mutable_properties()->mutable_nmr_isotropic_shielding_pbe0_aug_pcs_1()->mutable_values() = source.properties().nmr_isotropic_shielding_pbe0_aug_pcs_1().values();
  for (const double existing_value : source.properties().nmr_isotropic_shielding_pbe0_aug_pcs_1().values()) {
    destination.mutable_properties()->mutable_nmr_isotropic_shielding_pbe0_aug_pcs_1()->add_values(existing_value);
  }
  // single_point_energy_atomic_b5
  // nmr_isotropic_shielding_pbe0_aug_pcs_1
}

int
GetLowEnergyNmr(const GoogleSmu::Molecule& conformer,
            JobOptions& options,
            std::unordered_map<int,  GoogleSmu::Molecule>& low_energy) {
  if (conformer.properties().errors().fate() != GoogleSmu::Properties::FATE_SUCCESS) {
    return 1;
  }
  if (! conformer.properties().has_nmr_isotropic_shielding_pbe0_aug_pcs_1()) {
    return 1;
  }
  //if (! conformer.properties().has_partial_charges_esp_fit_pbe0_aug_pc_1()) {
  //  return 1;
  //}
  if (conformer.bond_topologies().size() == 0) {
    return 1;
  }

  if (options.btid_to_fetch.size() > 0) {
    int btid = conformer.molecule_id() / 1000;
    if (auto iter = options.btid_to_fetch.find(btid);
        iter == options.btid_to_fetch.end()) {
      return 1;
    }
  }
      

  options.results_obtained++;

  const int bond_topology_id = conformer.molecule_id() / 1000;
  const auto iter = low_energy.find(bond_topology_id);
  // If this is higher energy than anything previously found for this bond_topology_id,
  // not interested.
  if (iter != low_energy.end()) {
    if (conformer.properties().single_point_energy_atomic_b5().value() >
        iter->second.properties().single_point_energy_atomic_b5().value()) {
      return 1;
    }
    // Update info stored.
    CopyData(conformer, iter->second);
  } else {  // emplace new information.
    GoogleSmu::Molecule to_store;
    CopyData(conformer, to_store);
    low_energy.emplace(bond_topology_id, std::move(to_store));
  }

  return 1;
}
                
int
GetLowEnergyNmr(const const_IWSubstring& data,
            JobOptions& options,
            std::unordered_map<int,  GoogleSmu::Molecule>& low_energy) {
  const std::string as_string(data.data(), data.length());
  GoogleSmu::Molecule conformer;
  if (! conformer.ParseFromString(as_string)) {
    cerr << "Cannot decode proto\n";
    return 0;
  }
  options.items_read++;

  return GetLowEnergyNmr(conformer, options, low_energy);
}

int
GetLowEnergyNmr(iw_tf_data_record::TFDataReader& reader,
       JobOptions& options,
       std::unordered_map<int,  GoogleSmu::Molecule>& low_energy) {
  while (reader.good() && ! reader.eof()) {
    std::optional<const_IWSubstring> data = reader.Next();
    if (! data) {
      return 1;
    }
    if (! GetLowEnergyNmr(*data, options, low_energy)) {
      cerr << "Cannot process item " << reader.items_read() << '\n';
      return 0;
    }

    if (options.report_progress()) {
      cerr << "Read " << options.report_progress.times_called() << " items\n";
    }

    if (options.items_read > options.nprocess) {
      return 1;
    }
  }

  return 1;
}

int
GetLowEnergyNmr(const char * fname,
            JobOptions& options,
            std::unordered_map<int, GoogleSmu::Molecule>& low_energy) {
  iw_tf_data_record::TFDataReader reader(fname);
  if (! reader.good()) {
    cerr << "Cannot open " << fname << '\n';
    return 0;
  }

  return GetLowEnergyNmr(reader, options, low_energy);
}

int
ReadBtids(iwstring_data_source& input,
          std::unordered_set<int> & btid_to_fetch) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    int btid;
    if (! buffer.numeric_value(btid) || btid <= 0) {
      cerr << "ReadConformerIds:invalid conformer id '" << buffer << "'\n";
      return 0;
    }

    btid_to_fetch.insert(btid);
  }
  return btid_to_fetch.size();
}

int
ReadBtids(const char * fname, 
          std::unordered_set<int> & btid_to_fetch) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "ReadConformerIds:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadBtids(input, btid_to_fetch);
}

int
GetNMR(int argc, char ** argv) {
  Command_Line cl(argc, argv, "vr:O:S:G:I:s:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  JobOptions options;

  options.verbose = cl.option_count('v');

  if (cl.option_present('r')) {
    if (! options.report_progress.initialise(cl, 'r', options.verbose)) {
      cerr << "Cannot initialise progress reporting\n";
    }
  }

  if (cl.option_present('I')) {
    const_IWSubstring token;
    for (int i = 0; cl.value('I', token, i); ++i) {
      if (! options.BuildIgnoreDirective(token)) {
        cerr << "Invalid ignore directive '" << token << "'\n";
        return 0;
      }
    }
  }

  if (cl.option_present('G')) {
    const char * fname = cl.option_value('G');
    if (! ReadBtids(fname, options.btid_to_fetch)) {
      cerr << "Cannot read conformer ID's to process '" << fname << "'\n";
      return 1;
    }
  }

  if (cl.option_present('s')) {
    double s;
    if (! cl.value('s', s) || s <= 0.0) {
      cerr << "The scaling factor (-s) must be a +ve real value\n";
      return 1;
    }
    if (options.verbose) {
      cerr << "Scale to int " << s << '\n';
    }
    options.SetScalingFactor(s);
  }

  // First scan the input data and for each BondTopology get the lowest
  // energy Conformer.
  std::unordered_map<int, GoogleSmu::Molecule> low_energy;

  for (const char * fname : cl) {
    if (! GetLowEnergyNmr(fname, options, low_energy)) {
      cerr << "Fatal error processing " << fname << '\n';
      return 1;
    }
  }

  cerr << "Low energy conformers found, size " << low_energy.size() << '\n';

  constexpr char sep = ' ';
  IWString_and_File_Descriptor output(1);

  // Per atomic number accumulation of shielding values.
  std::vector<std::unordered_map<int, int>> shielding(10);

  for (const auto& [btid, conformer] : low_energy) {
    int startingbt = smu::IndexOfStartingBTid(conformer.bond_topologies());
    output << conformer.bond_topologies(startingbt).smiles() << sep 
           << conformer.molecule_id() << sep
           << conformer.properties().single_point_energy_atomic_b5().value() << sep
           << conformer.properties().errors().fate() << sep
           << conformer.bond_topologies().size() << '\n';
    const BondTopology& bt0 = conformer.bond_topologies(startingbt);
    const int natoms = bt0.atoms().size();
    //cerr << "Molecule has " << natoms << " atoms\n";
    for (int i = 0; i < natoms; ++i) {
      const int atomic_number = ToAtomicNumber(bt0.atoms(i));
      // cerr << "Checking ignore " << i << " atomic_number " << atomic_number << '\n';
      if (options.Ignore(bt0, i, atomic_number)) {
        continue;
      }
      // cerr << "Not ignored, fetching value " << i << '\n';
      double s = conformer.properties().nmr_isotropic_shielding_pbe0_aug_pcs_1().values(i);
      int j = static_cast<int>(s * options.scaling_factor);
      // cerr << " atomic_number " << atomic_number << " j " << j << '\n';
      shielding[atomic_number][j] += 1;
      // cerr << "Updated hash\n";
    }
  }
  output.flush();
  cerr << "Done processing low energy conformers\n";

  IWString stem_for_distributions = "S";
  if (cl.option_present('S')) {
    cl.value('S', stem_for_distributions);
  }

  // For each atomic number, the first and 99th quantile.
  std::vector<std::pair<int, int>> quantiles(10);

  cerr << "AtomicNumber" << sep << "values" << sep << "N" << sep << "min" << sep << "q1" << sep << "median" << sep << "q99" << sep << "max\n";
  for (unsigned int atomic_number = 0; atomic_number < shielding.size(); ++atomic_number) {
    const std::unordered_map<int, int>& shieldingi = shielding[atomic_number];
    const auto nvalues = shieldingi.size();
    if (nvalues == 0) {
      continue;
    }
    // Need to sort by shielding value for output.
    std::vector<std::pair<int, int>> values;
    values.reserve(nvalues);
    int tot = 0;
    // Iterate over shielding value and count.
    for (const auto [s, c] : shieldingi) {
      values.emplace_back(std::pair<int, int>(s, c));
      tot += c;
    }
    std::sort(values.begin(), values.end(), [](const std::pair<int, int>& p1, const std::pair<int, int>& p2) {
      return p1.first < p2.first;
    });
    int q1 = tot / 100;
    int q2 = tot - q1;
    int q1_value = -1;
    int q2_value = -1;
    int median_value = -1;
    bool q1_set = false;
    bool q2_set = false;
    bool median_set = false;
    int cumulative_count = 0;
    int median_count = tot / 2;
    IWString fname;
    fname << stem_for_distributions << atomic_number << ".txt";
    IWString_and_File_Descriptor output;
    if (! output.open(fname)) {
      cerr << "Cannot open stream for shielding distribution '" << fname << "'\n";
      return 1;
    }
    output << "Shielding" << sep << "count" << sep << "cumulative\n";
    for (unsigned int j = 0; j < nvalues; ++j) {
      cumulative_count += values[j].second;
      output << values[j].first << sep << values[j].second << sep << cumulative_count << '\n';
      if (! q1_set && cumulative_count >= q1) {
        q1_value = values[j].first;
        q1_set = true;
      }
      if (! q2_set && cumulative_count >= q2) {
        q2_value = values[j].first;
        q2_set = true;
      }
      if (! median_set && cumulative_count >= median_count) {
        median_value = values[j].first;
        median_set = true;
      }
    }
    quantiles[atomic_number].first = q1_value;
    quantiles[atomic_number].second = q2_value;
    cerr << atomic_number << sep << nvalues << sep << cumulative_count << sep << values[0].first << sep << q1_value << sep << median_value << sep << q2_value << sep << values[nvalues -1].first << '\n';
  }

  if (cl.option_present('O')) {
    IWString fname = cl.option_value('O');
    IWString_and_File_Descriptor stream_for_outliers;
    if (! stream_for_outliers.open(fname)) {
      cerr << "Cannot open stream for outliers '" << fname << "'\n";
      return 1;
    }

    for (const auto& [btid, conformer] : low_energy) {
      WriteIfExtremeValue(options, conformer, quantiles, stream_for_outliers);
    }
  }

  if (options.verbose) {
    cerr << "Read " << options.items_read << " Conformers, " << options.results_obtained << " examined\n";
  }

  return 0;
}

}  // namespace get_nmr

int
main(int argc, char ** argv)
{
  int rc = get_nmr::GetNMR(argc, argv);

  return rc;
}
