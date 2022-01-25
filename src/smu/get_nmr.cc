// Extract NMR data from SMU

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

struct JobOptions {
  int verbose = 0;

  int items_read = 0;

  // The number of Conformers that pass the inclusion tests.
  int results_obtained = 0;

  int nprocess = std::numeric_limits<int>::max();

  Report_Progress report_progress;
  
  std::unordered_set<int> btid_to_fetch;
};

void
Usage(int rc) {
  exit(rc);
}

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

int
WriteIfExtremeValue(const Conformer& conformer,
                    const std::vector<std::pair<int, int>> & quantiles,
                    IWString_and_File_Descriptor& output) {
  const BondTopology& bt0 = conformer.bond_topologies(0);
  const auto& nmr = conformer.properties().nmr_isotropic_shielding_pbe0_6_31ppgdp().values();
  const int matoms = bt0.atoms().size();
  resizable_array<int> outlier_atoms;
  resizable_array<int> outlier_values;
  for (int i = 0; i < matoms; ++i) {
    const int atomic_number = ToAtomicNumber(bt0.atoms(i));
    const int shielding = static_cast<int>(nmr[i]);
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
    int v = std::abs(outlier_values[i]);
    maybe_mol->set_isotope(outlier_atoms[i], v);
    if (v > max_outlier_value) {
      max_outlier_value = v;
      outlier_atomic_symbol = maybe_mol->atomic_symbol(outlier_atoms[i]);
    }
  }

  constexpr char sep = ' ';

  output << maybe_mol->smiles() << sep
         << conformer.conformer_id() << sep
         << (conformer.conformer_id() / 1000) << sep
         << max_outlier_value << sep
         << outlier_atomic_symbol << sep
         << conformer.properties().single_point_energy_pbe0d3_6_311gd().value()
         << '\n';

  return 1;
}

void
CopyData(const Conformer& source,
         Conformer& destination) {
  destination.set_conformer_id(source.conformer_id());
  destination.set_fate(source.fate());
  for (const auto & existing_bt : source.bond_topologies()) {
    BondTopology* bt = destination.add_bond_topologies();
    *bt = existing_bt;
  }

  destination.mutable_properties()->mutable_single_point_energy_pbe0d3_6_311gd()->set_value(source.properties().single_point_energy_pbe0d3_6_311gd().value());
  for (const double existing_nmr : source.properties().nmr_isotropic_shielding_pbe0_6_31ppgdp().values()) {
    destination.mutable_properties()->mutable_nmr_isotropic_shielding_pbe0_6_31ppgdp()->add_values(existing_nmr);
  }
  // single_point_energy_pbe0d3_6_311gd
  // nmr_isotropic_shielding_pbe0_6_31ppgdp
}

int
GetLowEnergyNmr(const Conformer& conformer,
            JobOptions& options,
            std::unordered_map<int,  Conformer>& low_energy) {
  if (conformer.fate() != Conformer::FATE_SUCCESS) {
    return 1;
  }
  if (! conformer.properties().has_nmr_isotropic_shielding_pbe0_6_31ppgdp()) {
    return 1;
  }
  if (! conformer.properties().has_single_point_energy_pbe0d3_6_311gd()) {
    return 1;
  }
  if (conformer.bond_topologies().size() == 0) {
    return 1;
  }

  if (options.btid_to_fetch.size() > 0) {
    int btid = conformer.conformer_id() / 1000;
    if (auto iter = options.btid_to_fetch.find(btid);
        iter == options.btid_to_fetch.end()) {
      return 1;
    }
  }
      

  options.results_obtained++;

  const int bond_topology_id = conformer.conformer_id() / 1000;
  const auto iter = low_energy.find(bond_topology_id);
  // If this is higher energy than anything previously found for this bond_topology_id,
  // not interested.
  if (iter != low_energy.end()) {
    if (conformer.properties().single_point_energy_pbe0d3_6_311gd().value() >
        iter->second.properties().single_point_energy_pbe0d3_6_311gd().value()) {
      return 1;
    }
    // Update info stored.
    CopyData(conformer, iter->second);
  } else {  // emplace new information.
    Conformer to_store;
    CopyData(conformer, to_store);
    low_energy.emplace(bond_topology_id, std::move(to_store));
  }

  return 1;
}
                
int
GetLowEnergyNmr(const const_IWSubstring& data,
            JobOptions& options,
            std::unordered_map<int,  Conformer>& low_energy) {
  const std::string as_string(data.data(), data.length());
  Conformer conformer;
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
       std::unordered_map<int,  Conformer>& low_energy) {
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
            std::unordered_map<int, Conformer>& low_energy) {
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
  Command_Line cl(argc, argv, "vr:O:S:G:");
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

  if (cl.option_present('G')) {
    const char * fname = cl.option_value('G');
    if (! ReadBtids(fname, options.btid_to_fetch)) {
      cerr << "Cannot read conformer ID's to process '" << fname << "'\n";
      return 1;
    }
  }

  std::unordered_map<int, Conformer> low_energy;

  for (const char * fname : cl) {
    if (! GetLowEnergyNmr(fname, options, low_energy)) {
      cerr << "Fatal error processing " << fname << '\n';
      return 1;
    }
  }

  constexpr char sep = ' ';
  IWString_and_File_Descriptor output(1);

  // Per atomic number accumulation of shielding values.
  std::vector<std::unordered_map<int, int>> shielding(10);

  for (const auto& [btid, conformer] : low_energy) {
    output << conformer.bond_topologies(0).smiles() << sep 
           << conformer.conformer_id() << sep
           << conformer.properties().single_point_energy_pbe0d3_6_311gd().value() << sep
           << conformer.fate() << sep
           << conformer.bond_topologies().size() << '\n';
    const BondTopology& bt0 = conformer.bond_topologies(0);
    const int natoms = bt0.atoms().size();
    for (int i = 0; i < natoms; ++i) {
      const int atomic_number = ToAtomicNumber(bt0.atoms(i));
      double s = conformer.properties().nmr_isotropic_shielding_pbe0_6_31ppgdp().values(i);
      int j = static_cast<int>(s);
      shielding[atomic_number][j] += 1;
    }
  }
  output.flush();

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
      WriteIfExtremeValue(conformer, quantiles, stream_for_outliers);
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
