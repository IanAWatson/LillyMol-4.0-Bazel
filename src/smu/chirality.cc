// Discern different chirality within conformers of the same BondTopology

#include <stdlib.h>
#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/is_actually_chiral.h"
#include "Molecule_Lib/molecule.h"

namespace smu_chirality {

struct JobOptions {
  int verbose = 0;

  int molecules_read = 0;

  int bond_topology_ids_read = 0;

  resizable_array<int> fates_processed;

  int no_chirality_discerned = 0;

  extending_resizable_array<int> number_conformers;
  Accumulator_Int<int> acc_number_conformers;
  extending_resizable_array<int> chiral_centre_count;
  extending_resizable_array<int> chiral_variants;

  public:
    JobOptions() {
    }

    int Build(Command_Line& cl);

    int OkToProcess(Molecule& m);

    int Report(std::ostream& output) const;
};

int
JobOptions::Build(Command_Line& cl) {
  verbose = cl.option_count('v');

  if (cl.option_present('F')) {
    int f;
    for (int i = 0; cl.value('F', f, i); ++i) {
      fates_processed << f;
    }
  }

  return 1;
}

int
JobOptions::Report(std::ostream& output) const {
  output << "Read " <<  molecules_read << " molecules\n";
  if (molecules_read == 0) {
    return output.good();
  }

  output << "Read " << bond_topology_ids_read << " bond topology ids\n";

  output << "btids have between " << acc_number_conformers.minval() << " and " <<
            acc_number_conformers.maxval() << " mean " << acc_number_conformers.average() << 
            " conformers\n";
  for (int i = 0; i < number_conformers.number_elements(); ++i) {
    if (number_conformers[i]) {
      cerr << number_conformers[i] << " bond topology's had " << i << " conformers\n";
    }
  }
  for (int i = 0; i < chiral_centre_count.number_elements(); ++i) {
    if (chiral_centre_count[i]) {
      cerr << chiral_centre_count[i] << " bond topology's had " << i << " chiral centres\n";
    }
  }
  for (int i = 0; i < chiral_variants.number_elements(); ++i) {
    if (chiral_variants[i]) {
      cerr << chiral_variants[i] << " bond topology's had " << i << " chiral variants\n";
    }
  }
  output << no_chirality_discerned << " btids had no chirality discerned\n";
  return output.good();
}

int
JobOptions::OkToProcess(Molecule& m) {
  if (fates_processed.empty()) {
    return 1;
  }

  const_IWSubstring token;
  if (! m.name().word(2, token)) {
    cerr << "Cannot extract fate from " << m.name() << '\n';
    return 0;
  }

  int fate;
  if (! token.numeric_value(fate) || fate < 1) {
    cerr << "Invalid fate " << m.name() << '\n';
    return 0;
  }

  return fates_processed.contains(fate);
}

void
Usage(int rc) {
  ::exit(rc);
}

IWString
GetBondTopologyId(const IWString& name) {
  IWString result;
  if (! name.word(0, result)) {
    cerr << "GetBondTopologyId:cannot extract first word from " << name << '\n';
    return result;
  }

#ifdef CONFORMIDS_IN_INPUT
  for (int i = 0; i < 3; ++i) {
    char c = result.back();
    if (! isdigit(c)) {
      cerr << "GetBondTopologyId:non numeric " << name << '\n';
      return result;
    }
    result.pop();
  }
#endif

  return result;
}

int
Chirality(resizable_array_p<Molecule>& btids,
          JobOptions& options,
          const IWString& btid,
          IWString_and_File_Descriptor& output) {
#ifdef DEBUG_CHIRALITY
  cerr << "Group of " << btids.number_elements() << "\n";
  for (const Molecule* m : btids) {
    cerr << m->name() << '\n';
  }
#endif
  options.bond_topology_ids_read++;

  const int n = btids.number_elements();
  int matoms = btids[0]->natoms();
  for (int i = 1; i < n; ++i) {
    if (btids[i]->natoms() != matoms) {
      cerr << "Chirality:atom count mismatch " << matoms << " vs " << btids[i]->natoms() << '\n';
      return 0;
    }
  }

  int * chiral = new_int(matoms); std::unique_ptr<int[]> free_chiral(chiral);

  int chiral_centres_found = 0;
  for (int i = 0; i < matoms; ++i) {
    if (is_actually_chiral(*btids[0], i)) {
      chiral[i] = 1;
      ++chiral_centres_found;
    }
  }

  options.chiral_centre_count[chiral_centres_found]++;

#ifdef DEBUG_CHIRALITY
  cerr << "Found " << chiral_centres_found << " chiral centres\n";
#endif
  if (chiral_centres_found == 0) {
    return 1;
  }

  IW_STL_Hash_Map_int smiles_count;

  int chirality_not_discerned = 0;
  for (Molecule* m : btids) {
    m->discern_chirality_from_3d_structure();
    if (m->chiral_centres() == 0) {
      ++chirality_not_discerned;
      continue;
    }
    const IWString& s = m->unique_smiles();
    auto iter = smiles_count.find(s);
    if (iter == smiles_count.end()) {
      smiles_count[s] = 1;
    } else {
      ++iter->second;
    }
  }

  options.chiral_variants[smiles_count.size()]++;

#ifdef DEBUG_CHIRALITY
  cerr << "After enumeration have " << smiles_count.size() << " different smiles\n";
#endif
  if (chirality_not_discerned == btids.number_elements()) {
    options.no_chirality_discerned++;
    // cerr << btids[0]->smiles() << ' ' << btid << " no chirality discerned\n";
    return 1;
  }

  if (smiles_count.size() == 1) {
    return 1;
  }

  constexpr char sep = ' ';

  for (const auto& [smiles, count] : smiles_count) {
    output << smiles << sep << btid << sep << count << sep << chiral_centres_found << sep << smiles_count.size() << '\n';
  }
  output << "*\n";
  output.write_if_buffer_holds_more_than(32768);

  return 1;
}

void
Preprocess(Molecule& m) {
  m.remove_all(1);
}

int
Chirality(data_source_and_type<Molecule>& input,
          JobOptions& options,
          IWString_and_File_Descriptor& output) {
  // groups of molecules having the same BondTopologyID
  resizable_array_p<Molecule> btids;
  IWString current_btid;

  Molecule * m;
  while (( m = input.next_molecule()) != nullptr) {
    Preprocess(*m);
    options.molecules_read++;
    if (! options.OkToProcess(*m)) {
      delete m;
      continue;
    }
    IWString btid = GetBondTopologyId(m->name());
    if (btid == current_btid) {
      btids << m;
      continue;
    }

    if (current_btid.empty()) {
      current_btid = btid;
      btids << m;
      continue;
    }

    options.number_conformers[btids.number_elements()]++;
    options.acc_number_conformers.extra(btids.number_elements());

    if (btids.number_elements() > 1) {
      Chirality(btids, options, current_btid, output);
    }

    btids.resize_keep_storage(0);
    current_btid = btid;
  }

  return 1;
}


int
Chirality(const char * fname,
          FileType input_type,
          JobOptions& options,
          IWString_and_File_Descriptor& output) {
  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "Chirality:cannot open '" << fname << "'\n";
    return 0;
  }

  return Chirality(input, options, output);
}

int
Chirality(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:E:i:F:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  JobOptions options;
  if (! options.Build(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  FileType input_type = FILE_TYPE_SMI;

  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      cerr << "Cannot process input type\n";
      return 1;
    }
  } else {
  }

  IWString_and_File_Descriptor output(1);

  for (const char * fname: cl) {
    if (! Chirality(fname, input_type, options, output)) {
      cerr << "Chirality:fatal error processing '" << fname << "'\n";
      return 1;
    }
  }
  output.flush();

  if (options.verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace smu_chirality

int
main(int argc, char** argv) {
  int rc = smu_chirality::Chirality(argc, argv);

  return rc;
}
