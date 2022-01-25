// Examine the output of SMU BondTopology matching and report matches that
// differ in ring membership.
// 176028002,10,-433.388,C([NH:3][N:2]=[NH:1])=[N:4][O:5][F:6],0.0569,0.0,T,C(=[N:3][N:2]=[NH:1])=[N:4][O:5][F:6],0.0556,0.0,F,C([NH:3][N:2]=[NH:1])[NH:4][O:5][F:6],0.0389,0.0,F,C(=[N:3][N:2]=[NH:1])[NH:4][O:5][F:6],0.0375,0.0,F,C([NH:3][NH:2][NH2:1])=[N:4][O:5][F:6],0.033,0.0,F,C(=[N:3][NH:2][NH2:1])=[N:4][O:5][F:6],0.0317,0.0,F,C([N:3]=[N:2][NH2:1])=[N:4][O:5][F:6],0.0314,0.0,F,C([NH:3][NH:2][NH2:1])[NH:4][O:5][F:6],0.015,0.0,F,C(=[N:3][NH:2][NH2:1])[NH:4][O:5][F:6],0.0136,0.0,F,C([N:3]=[N:2][NH2:1])[NH:4][O:5][F:6],0.0133,0.0,F

#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"
#include "Foundational/iwmisc/misc.h"
#include "Molecule_Lib/molecule.h"

namespace smu_ring_membership_change {

using std::cerr;

char input_separator = ',';

int molecules_read = 0;
int no_rings = 0;
int single_structure = 0;
int changed_ring_membership = 0;

constexpr int kElements = 10;   // Fluorine + 1

void
Usage(int rc) {
  exit(rc);
}

struct RingAtoms {
  int ring_atoms;
};

bool 
SameRingAtoms(const int * r1, const int * r2, int n) {
  for (int i = 0; i < n; ++i) {
    if (r1[i] != r2[i]) {
      return false;
    }
  }

  return true;
}

int
UpdateRingAtomCount(Molecule& m,
                    int * atoms_in_ring) {
  if (m.nrings() == 0) {
    return 0;
  }

  const int matoms = m.natoms();
  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (! m.is_ring_atom(i)) {
      continue;
    }
    atoms_in_ring[m.atomic_number(i)]++;
    rc++;
  }

  return rc;
}

int
SmuRingMembershipChange(const IWString& conformer_id,
                        Molecule* mols,
                        const int nmols,
                        const int * nrings,
                        std::ostream & output) {
  std::vector<int> elements {6, 7, 8, 9};
  int found_mismatch = 0;
  IWString output_buffer;
  for (int i = 0; i < nmols; ++i) {
    for (int j = i + 1; j < nmols; ++j) {
      if (SameRingAtoms(nrings + i * kElements, nrings + j * kElements, kElements)) {
        continue;
      }
      found_mismatch++;
      output_buffer << mols[i].unique_smiles() << '.' << mols[j].unique_smiles() << ' ' << conformer_id << '\n';
    }
  }
  if (found_mismatch) {
//  output << conformer_id << ' ' << found_mismatch << '\n';
    output << output_buffer;
    changed_ring_membership++;
  }
  return 1;
}

int
SmuRingMembershipChange(const IWString& conformer_id,
                        Molecule* mols,
                        const int nmols,
                        std::ostream & output) {
  std::unique_ptr<int[]> ring_atoms(new_int(nmols * kElements));
  int ring_atoms_detected = 0;
  for (int i = 0; i < nmols; ++i) {
    if (UpdateRingAtomCount(mols[i], ring_atoms.get() + i * kElements)) {
      ring_atoms_detected++;
    }
  }

  if (ring_atoms_detected == 0) {
    no_rings++;
    return 1;
  }

  return SmuRingMembershipChange(conformer_id, mols, nmols, ring_atoms.get(), output);
}

int
UniqueMolecules(Molecule* mols, int nmols) {
  int ndx = 0;
  IW_STL_Hash_Set seen;
  for (int i = 0; i < nmols; ++i) {
    if (seen.contains(mols[i].unique_smiles())) {
      continue;
    }
    seen.emplace(mols[i].unique_smiles());
    if (ndx < i) {
      mols[ndx] = std::move(mols[i]);
    }
    ndx++;
  }

  return ndx;
}

int
SmuRingMembershipChangeLine(const const_IWSubstring& buffer,
                        std::ostream & output) {
  molecules_read++;

  IWString conformer_id;
  int pos = 0;
  buffer.nextword(conformer_id, pos, input_separator);
  IWString token;
  buffer.nextword(token, pos, input_separator);
  int nstructures;
  if (! token.numeric_value(nstructures) || nstructures < 1) {
    cerr << "Invalid number of matches\n";
    return 0;
  }

  if (nstructures == 1) {
    single_structure++;
    return 1;
  }

  buffer.nextword(token, pos, input_separator);  // Skip over energy

  std::unique_ptr<Molecule[]> mols = std::make_unique<Molecule[]>(nstructures);
  for (int i = 0; i < nstructures; ++i) {
    IWString smiles;
    buffer.nextword(smiles, pos, input_separator);
    if (! mols[i].build_from_smiles(smiles)) {
      cerr << "Invalid smiles " << smiles << '\n';
      return 0;
    }
    buffer.nextword(token, pos, input_separator);  // Score
    buffer.nextword(token, pos, input_separator);  // Probability
    buffer.nextword(token, pos, input_separator);  // T/F
  }

  nstructures = UniqueMolecules(mols.get(), nstructures);
  if (nstructures == 1) {
    single_structure++;
    return 1;
  }

  return SmuRingMembershipChange(conformer_id, mols.get(), nstructures, output);
}

int
SmuRingMembershipChange(iwstring_data_source& input,
                        std::ostream & output) {
  const_IWSubstring buffer;
  input.next_record(buffer);  // skip header
  while (input.next_record(buffer)) {
    if (! SmuRingMembershipChangeLine(buffer, output)) {
      cerr << "Fatal error processing " << buffer << '\n';
      return 0;
    }
  }

  return 1;
}

int
SmuRingMembershipChange(const char * fname,
                        std::ostream & output) {
  iwstring_data_source input(fname);
  if (! input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return SmuRingMembershipChange(input, output);
}

int
SmuRingMembershipChange(int argc, char ** argv) {
  Command_Line cl(argc, argv, "vA:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (cl.number_elements() == 0) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  if (! SmuRingMembershipChange(cl[0], std::cout)) {
    return 1;
  }

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules, " << no_rings << " contained no rings\n";
    cerr << single_structure << " contained just a single structure\n";
    cerr << changed_ring_membership << " conformers with changed ring membership\n";
  }

  return 0;
}

}  // namespace smu_ring_membership_change

int
main(int argc, char ** argv) {
  return smu_ring_membership_change::SmuRingMembershipChange(argc, argv);
}
