// Analyse duplicates from the SMU pipeline
// Typical input record
// 878962,2,OC1=CC=C2C=C12,31.071,C1=[C:4]2[CH:1]=[CH:3][C:2]([OH:6])=[C:5]12,F,11.813,C1=[C:5]2[C:2]([OH:6])=[CH:3][CH:1]=[C:4]12,T

#include <iostream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/smiles.h"

namespace smu_duplicates {

int lines_read = 0;
int one_topology_found = 0;
int one_topology_is_start = 0;
int one_topology_not_found = 0;
int multiple_bond_topologies = 0;
int same_smiles = 0;
int same_unique_smiles = 0;

int found1 = 0;
int found2 = 0;
int notfound = 0;

int verbose = 0;

char input_separator = ',';
char output_separator = ',';

using std::cerr;

void
Usage(int rc) {
  exit(rc);
}

class MoleculeAndSmiles {
  private:
    Molecule _mol;
    IWString _smiles;
  public:
    MoleculeAndSmiles(const const_IWSubstring& smi);

    int natoms() const { return _mol.natoms();}
    Molecule& mol() { return _mol;}

    const IWString& smiles() const { return _smiles;}
};

MoleculeAndSmiles::MoleculeAndSmiles(const const_IWSubstring& smi) : _smiles(smi) {
  if (! _mol.build_from_smiles(smi)) {
    cerr << "MoleculeAndSmiles:invalid smiles " << smi << '\n';
    return;
  }
}

std::optional<MoleculeAndSmiles>
MakeMoleculeAndSmiles(const const_IWSubstring& smi) {
  MoleculeAndSmiles result(smi);
  if (result.mol().natoms() == 0) {
    return std::nullopt;
  }

  return result;
}

int
SmuDuplicates(const const_IWSubstring& id,
              MoleculeAndSmiles& starting_molecule,
              int found,
              MoleculeAndSmiles& ms1,
              MoleculeAndSmiles& ms2,
              std::ostream& output) {
  Molecule& m1 = ms1.mol();
  Molecule& m2 = ms2.mol();

  if (m1.smiles() == m2.smiles()) {
    same_smiles++;
    output << starting_molecule.smiles() << output_separator << ms1.smiles() << output_separator << ms2.smiles() << output_separator << m1.unique_smiles() << output_separator << id << output_separator << found << output_separator << "SS\n";
    return 1;
  }

  if (m1.unique_smiles() == m2.unique_smiles()) {
    same_unique_smiles++;
    output << starting_molecule.smiles() << output_separator << ms1.smiles() << output_separator << ms2.smiles() << output_separator << m1.unique_smiles() << output_separator << id << output_separator << found << output_separator << "SUS\n";
    return 1;
  }

  output << starting_molecule.smiles() << output_separator << ms1.smiles() << output_separator << ms2.smiles() << output_separator << m1.unique_smiles() << output_separator << m2.unique_smiles() << output_separator << id << output_separator << found << output_separator << "DUS\n";
  return 1;
}

int
SmuDuplicates(const const_IWSubstring& id,
              const const_IWSubstring& starting_smiles,
              int found,
              const const_IWSubstring& smiles1,
              const const_IWSubstring& smiles2,
              std::ostream& output) {
  std::optional<MoleculeAndSmiles> starting_molecule = MakeMoleculeAndSmiles(starting_smiles);
  if (! starting_molecule) {
    cerr << "Invalid starting smiles " << starting_smiles << '\n';
    return 0;
  }

  std::optional<MoleculeAndSmiles> ms1 = MakeMoleculeAndSmiles(smiles1);
  std::optional<MoleculeAndSmiles> ms2 = MakeMoleculeAndSmiles(smiles2);
  if (! ms1 || ! ms2) {
    cerr << "Cannot interpret smiles\n";
    cerr << smiles1 << '\n';
    cerr << smiles2 << '\n';
    return 0;
  }

  set_include_atom_map_with_smiles(0);

  const int natoms = ms1->natoms();
  Molecule& mol1 = ms1->mol();
  Molecule& mol2 = ms2->mol();

  for (int i = 0; i < natoms; ++i) {
    mol1.unset_all_implicit_hydrogen_information(i);
    mol2.unset_all_implicit_hydrogen_information(i);
  }

  mol1.recompute_implicit_hydrogens();
  mol1.recompute_implicit_hydrogens();

  return SmuDuplicates(id, *starting_molecule, found, *ms1, *ms2, output);
}

int
SmuDuplicatesOneMatch(const const_IWSubstring& id,
                      MoleculeAndSmiles& starting_molecule,
                      MoleculeAndSmiles& found_molecule,
                      std::ostream& output) {
  output << starting_molecule.smiles() << output_separator << found_molecule.smiles() << output_separator << id << output_separator << "NF1\n";

  return 1;
}

int
SmuDuplicatesOneMatch(const const_IWSubstring& id,
                      const const_IWSubstring& starting_smiles,
                      const const_IWSubstring& found_smiles,
                      std::ostream& output) {
  std::optional<MoleculeAndSmiles> starting_molecule = MakeMoleculeAndSmiles(starting_smiles);
  if (! starting_molecule) {
    cerr << "Invalid starting molecule smiles " << starting_smiles << '\n';
    return 1;
  }

  std::optional<MoleculeAndSmiles> found_molecule = MakeMoleculeAndSmiles(found_smiles);
  if (! found_molecule) {
    cerr << "Invalid found molecule smiles " << found_smiles << '\n';
    return 1;
  }

  return SmuDuplicatesOneMatch(id, *starting_molecule, *found_molecule, output);
}

int
SmuDuplicatesOneMatch(const const_IWSubstring& buffer,
                      std::ostream & output) {
  one_topology_found++;
  const_IWSubstring id;
  int i = 0;
  const_IWSubstring token;
  const_IWSubstring starting_smiles;
  const_IWSubstring found_smiles;
  const_IWSubstring tf;
  for (int word = 0; buffer.nextword(token, i, input_separator); ++word) {
    if (word == 0) {
      id = token;
    } else if (word == 2) {
      starting_smiles = token;
    } else if (word == 4) {
      found_smiles = token;
    } else if (word == 5) {
      tf = token;
    }
  }

  if (tf == 'T') {
    one_topology_is_start++;
    return 1;
  }

  one_topology_not_found++;

  return SmuDuplicatesOneMatch(id, starting_smiles, found_smiles, output);
}

int
SmuDuplicatesRecord(const const_IWSubstring& buffer,
              std::ostream& output) {
  if (buffer.nwords(input_separator) == 6) {
    return SmuDuplicatesOneMatch(buffer, output);
  }

  if (buffer.nwords(input_separator) != 9) {
    return 1;
  }

  multiple_bond_topologies++;

  const_IWSubstring id;
  const_IWSubstring starting_smiles;
  const_IWSubstring smiles1, smiles2;
  const_IWSubstring tf1, tf2;
  int i = 0;
  const_IWSubstring token;
  for (int word = 0; buffer.nextword(token, i, input_separator); ++word) {
    if (word == 0) {
      id = token;
    } else if (word == 2) {
      starting_smiles = token;
    } else if (word == 4) {
      smiles1 = token;
    } else if (word == 5) {
      tf1 = token;
    } else if (word == 7) {
      smiles2 = token;
    } else if (word == 8) {
      tf2 = token;
    }
  }

  if (smiles1.empty() || smiles2.empty()) {
    cerr << "No data\n";
    return 0;
  }

  int found;
  if (tf1 == 'T') {
    found1++;
    found = 1;
  } else if (tf2 == 'T') {
    found2++;
    found = 2;
  } else {
    notfound++;
    found = 3;
  }

  return SmuDuplicates(id, starting_smiles, found, smiles1, smiles2, output);
}

int
SmuDuplicates(iwstring_data_source& input,
              std::ostream& output) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    lines_read++;
    if (! SmuDuplicatesRecord(buffer, output)) {
      cerr << "Cannot parse '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;

}
int
SmuDuplicates(const char * fname,
              std::ostream& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return SmuDuplicates(input, output);
}

int
SmuDuplicates(int argc, char ** argv) {
  Command_Line cl(argc, argv, "vi:A:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  verbose = cl.option_count('v');

  set_global_aromaticity_type(EVERYTHING_HAS_A_PI_ELECTRON);
  set_default_unique_smiles_aromaticity(EVERYTHING_HAS_A_PI_ELECTRON);
  set_min_aromatic_ring_size(3);


  if (cl.empty()) {
    cerr << "MUst specify input file(s)\n";
    Usage(1);
  }

  for (const char * fname : cl) {
    if (! SmuDuplicates(fname, std::cout)) {
      cerr << "Cannot process '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    cerr << "Read " << lines_read << " lines, " << multiple_bond_topologies << " had multiple BondTopology's\n";
    cerr << one_topology_found << " one topology\n";
    cerr << one_topology_is_start << " only found topology is start\n";
    cerr << same_smiles << " had the same smiles\n";
    cerr << same_unique_smiles << " had the same unique smiles\n";
    cerr << found1 << " found as first  choice\n";
    cerr << found2 << " found as second choice\n";
    cerr << notfound << " not found\n";
  }

  return 0;
}

}

int
main(int argc, char ** argv) {
  return smu_duplicates::SmuDuplicates(argc, argv);
}
