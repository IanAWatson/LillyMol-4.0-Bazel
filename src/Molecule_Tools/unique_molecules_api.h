#ifndef MOLECULE_TOOLS_UNIQUE_MOLECULES_API_H
#define MOLECULE_TOOLS_UNIQUE_MOLECULES_API_H

#include <string>
#include <unordered_set>

#ifdef UNIQUE_MOLECULES_USES_ABSL
#include "absl/container/flat_hash_set.h"

#include "Foundational/iwstring/absl_hash.h"
#else
#include "Foundational/iwstring/iw_stl_hash_set.h"
#endif  // UNIQUE_MOLECULES_USES_ABSL

#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/iwreaction.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rmele.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/standardise.h"

namespace unique_molecules {

// Identifies unique molecules based on a hash of unique smiles.
// Typical usage:
//   UniqueMoleculesImplementation um;
//   Command_Line cl(argc, argv, "vIzclhaj");
//   um.Initialise(cl)
//   for (Molecule& m: molecules) {
//     if (umn.IsUnique(m)) {
//       do something with new molecules
//     } else {
//       do something with duplicate molecules
//     }
//   }
//   um.Report(cerr);  // if interested in statistics.
//
class UniqueMoleculesImplementation {
#ifdef UNIQUE_MOLECULES_USES_ABSL
  using UMHash = absl::flat_hash_set<IWString>;
#else
  using UMHash = IW_STL_Hash_Set;
#endif
  private:
    int _verbose;

    // The number of molecules examined vial IsUnique().
    int _molecules_processed;

    // The number of duplicates found by IsUnique().
    int _duplicates_found;

    // the number of unique molecules found by IsUnique().
    int _unique_molecules;

    // The various standardization actiona available.

    // If set, input molecules are transformed to their largest fragment.
    int _reduce_to_largest_fragment;

    // If set, chirality information is excluded from unique smiles formation.
    int _exclude_chiral_info;

    // If set, cis/trans information is excluded from unique smiles formation.
    int _exclude_cis_trans_bonding_info;

    // We can optionally ignore isotopic information.
    int _ignore_isotopes;

    Chemical_Standardisation _chemical_standardisation;

    // Useful for things like compressing the heavy halogens.
    Element_Transformations _element_transformations;

    // Specific elements can be removed before unique smiles generation.
    Elements_to_Remove _elements_to_remove;

    resizable_array_p<IWReaction> _reaction;
    int _molecules_changed_by_reactions;

    // A unique smiles hash for every atom count.
    resizable_array_p<UMHash> _smiles_hash;

    // An array of atom hashes from quick_bond_hash();
    resizable_array_p<std::unordered_set<uint64_t>> _atom_hash;
    int _use_atom_hash;

    // One possible comparison is to compare molecules via their graph form.
    int _compare_as_graph;

    // Timing results show that generally the molecular formula hash is not
    // worth computing.
    resizable_array_p<UMHash> _formula_hash;
    int _perform_formula_check;

    // A variant on unique molecule finding is that molecules are unique only
    // if the name field also matches.
    int _only_same_if_structure_and_name_the_same;

    // private functions
    int PerformReactions(Molecule& m);
    int Preprocess(Molecule& m);
    int IsUniqueInner(Molecule& m);
    void FormUniqueSmiles(Molecule& m, IWString& usmi) const;
    void FormUniqueSmilesInner(Molecule& m, IWString& usmi) const;

  public:
    UniqueMoleculesImplementation();

    // Setter functions.
    void set_verbose(int s) {
      _verbose = s;
    }
    void set_reduce_to_largest_fragment(int s) {
      _reduce_to_largest_fragment = s;
    }
    void set_exclude_chiral_info(int s) {
      _exclude_chiral_info = s;
    }
    void set_exclude_cis_trans_bonding_info(int s) {
      _exclude_cis_trans_bonding_info = s;
    }
    void set_ignore_isotopes(int s) {
      _ignore_isotopes = s;
    }

    int unique_molecules() const {
      return _unique_molecules;
    }

    // Initialise state from command line options.
    int Initialise(Command_Line& cl);

    // To activate all chemical standardisations, use 'all' for `directive`.
    int ActivateChemicalStandardization(const std::string& directive);

    // Add the information for `m` to the internal hashes.
    // It adds the unique smiles to the hash, and returns whether or not the molecule
    // is unique, but does not update any counters.
    int IngestPreviousMolecule(Molecule& m);

    // Returns 1 if `m` is unique.
    int IsUnique(Molecule& m);

    // Interpret `smiles` as a molecule and call IsUnique(m).
    int SmilesIsUnique(const std::string& smiles);

    int Report(std::ostream& output) const;
};

}  // namespace unique_molecules

#endif  // MOLECULE_TOOLS_UNIQUE_MOLECULES_API_H
