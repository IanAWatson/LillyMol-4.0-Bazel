#ifndef MOLECULE_TOOLS_DICER_API_H_
#define MOLECULE_TOOLS_DICER_API_H_

#include "absl/container/flat_hash_map.h"

#include "Foundational/iwbits/fixed_bit_vector.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iwstring/iwhash.h"

#include "Molecule_Lib/molecule.h"

#include "Molecule_Tools/dicer_api.pb.h"

namespace dicer_api {


struct AtomPair {
  atom_number_t a1;
  atom_number_t a2;

  AtomPair(atom_number_t c1, atom_number_t c2) : a1(c1), a2(c2) {
  }
};

struct CurrentMoleculeData {
    int natoms;
    // A list of AtomPair that are the breakable bonds.
    std::vector<AtomPair> breakable;

    // At each level we need an array to hold fragment membership
    std::vector<std::unique_ptr<int[]>> fragment_membership;

    // And at each level, an array that is used to hold subset information.
    std::vector<std::unique_ptr<int[]>> subset;

    // We keep track of which subsets of atoms have already been found.
    resizable_array_p<fixed_bit_vector::FixedBitVector> found;

  public:
    CurrentMoleculeData(int n);

    // Once we know how many bonds will be broken, we know the
    // recursion depth, and can therefore allocate the either_side vector.
    int SetMaxBondsBroken(int nb, int natoms);

    bool IsNew(std::unique_ptr<fixed_bit_vector::FixedBitVector>& candidate);
};

class DicerApi {
  private:
    int _min_fragment_size;
    int _max_fragment_size;
    int _max_bonds_to_break;
    int _break_amide_bonds;
    int _break_acid_bonds;
    int _break_cc_bonds;
    isotope_t _isotope;

    int _skip_bad_valence;

    absl::flat_hash_map<IWString, uint32_t, IWStringHash> _frag_count;
    absl::flat_hash_map<IWString, uint32_t, IWStringHash> _frag_to_ndx;

  // private functions.
    bool OkSize(int matoms) const;

    int IdentifyBondsToBreak(Molecule& m,
                CurrentMoleculeData& current_molecule_data);

    int MaybeFragment(Molecule& m, CurrentMoleculeData& current_molecule_data);
    int Process(Molecule& m,
                CurrentMoleculeData& current_molecule_data,
                int bonds_already_broken,
                uint32_t bstart,
                Sparse_Fingerprint_Creator& sfc);
    void ProcessFragments(Molecule& m,
                           CurrentMoleculeData& current_molecule_data,
                           int ndx,
                           Sparse_Fingerprint_Creator& sfc);
    int FinalProcessing(Molecule& m, CurrentMoleculeData& current_molecule_data,
                Sparse_Fingerprint_Creator& sfc);
    int ProcessFragment(Molecule& m,
                          CurrentMoleculeData& current_molecule_data,
                          const int* subset,
                          Sparse_Fingerprint_Creator& sfc);

  public:
    DicerApi();

    int Initialise(Command_Line& cl);

    int Process(Molecule& m, Sparse_Fingerprint_Creator& sfc);

    int Report(std::ostream& output) const;
};

}  // namespace dicer_api

#endif  // MOLECULE_TOOLS_DICER_API_H_
