#ifndef MOLECULE_TOOLS_SCAFFOLDS_H
#define MOLECULE_TOOLS_SCAFFOLDS_H

// Form all possible scaffold subsets for a given molecule.
// One of the criticisms of the Murcko scaffold is that the sequence
// methyl->ethyl->propyl->clclopropyl results in a different scaffold.
// This tool supports the idea of enumerating all possible scaffold
// subsets.

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/molecule.h"

#include "Molecule_Tools/scaffolds.pb.h"

namespace scaffolds {

// Avoid passing around a bunch of separate arguments.
class PerMoleculeData {
  private:
    // As returned by m.label_atoms_by_ring_system
    int* _ring_sys;

    // An indicator of whether or not an atom is in the current system
    // being assembled.
    int* _in_system;

    // The number of ring systems (returned by mlabel_atoms_by_ring_system).
    int _nsys;

    // For each non ring atom, a label which will be the same for all atoms
    // that are in the same inter-ring region.
    int* _region;

    // When building the subsets, we need an array for that.
    int* _atoms_in_subset;

    // Avoid duplicates
    IW_STL_Hash_Set _seen;

    // The path tracing functions need a temporary array to keep track
    // of which atoms have been visited.
    int* _tmp;

  public:
    PerMoleculeData(Molecule& m);
    ~PerMoleculeData();

    int nsys() const {
      return _nsys;
    }

    const int* ring_sys() const {
      return _ring_sys;
    }
    int* ring_sys() {
      return _ring_sys;
    }
    int ring_sys(int atom) const {
      return _ring_sys[atom];
    }
    const int* in_system() const {
      return _in_system;
    }

    int* in_system() {
      return _in_system;
    }
    int in_system(int atom) const {
      return _in_system[atom];
    }

    int* region() {
      return _region;
    }
    int region(int atom) const {
      return _region[atom];
    }

    int* atoms_in_subset() {
      return _atoms_in_subset;
    }

    // Return 1 if this molecule has been seen before.
    int Seen(Molecule& m);

    int* tmp() {
      return _tmp;
    }
};

class ScaffoldFinder {

  // All settable configuration options in the proto.
  Scaffolds::ScaffoldsOptions _config;

  // Keep track of the number of ring systems in the incoming molecules.
  extending_resizable_array<int> _nsys;
  // The number of scaffolds generated.
  extending_resizable_array<int> _generated;

  // We can specify the maximum number of ring systems lost in any particular scaffold
  // combination
  int _max_systems_lost;
  // We can also specify a minumim number of ring systems in a combination.
  int _min_systems_in_subset;

  private:

  // private functions.

    int OkSubsetSize(uint32_t in_molecule, uint32_t in_subset) const;

    int Process(Molecule& m,
                 Scaffolds::ScaffoldData& results);
    int Process(const Molecule& m,
                 const std::vector<int>& state,
                 PerMoleculeData& pmd,
                 const resizable_array<int>* systems_reachable,
                 Scaffolds::ScaffoldData& results);

    int AddToResultsIfUnique(Molecule& candidate,
                 PerMoleculeData& pmd,
                 int nsys,
                 Scaffolds::ScaffoldData& result);

  public:
    ScaffoldFinder();

    int Initialise(Command_Line& cl, char flag);
    int Initialise(const Scaffolds::ScaffoldsOptions& proto);

    // Place in `results` all scaffold subsets of `m`.
    int MakeScaffolds(Molecule& m, Scaffolds::ScaffoldData& results);

    int Report(std::ostream& output) const;
};

}  // namespace scaffolds

#endif // MOLECULE_TOOLS_SCAFFOLDS_H
