#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <unordered_set>

#include "Foundational/iwmisc/combinations.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"

#include "Molecule_Tools/scaffolds.h"
#include "Molecule_Tools/scaffolds.pb.h"

namespace scaffolds {

using std::cerr;

// `in_sys` identifies atons in the different ring systems in `m`.
// Identify =* groups and add them to `in_sys`.
int
ExtendToSinglyAttachedDoublyBonded(Molecule& m,
                                   int* in_sys) {
  const int matoms = m.natoms();
  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (in_sys[i]) {
      continue;;
    }

    const Atom& a = m.atom(i);
    if (a.ncon() != 1) {
      continue;
    }

    const Bond* b = a[0];
    if (! b->is_double_bond()) {
      continue;
    }

    atom_number_t j = b->other(i);
    if (in_sys[j] == 0) {
      continue;
    }

    in_sys[j] = in_sys[i];
    ++rc;
  }

  return rc;
}

PerMoleculeData::PerMoleculeData(Molecule& m) {
  const int matoms = m.natoms();
  _ring_sys = new int[matoms];
  _in_system = new int[matoms];
  _region = new int[matoms];

  _nsys = m.label_atoms_by_ring_system(_ring_sys);

  ExtendToSinglyAttachedDoublyBonded(m, _ring_sys);

  _atoms_in_subset = new int[matoms];

  _tmp = new int[matoms];
}

PerMoleculeData::~PerMoleculeData() {
  delete [] _ring_sys;
  delete [] _in_system;
  delete [] _tmp;
  delete [] _region;
  delete [] _atoms_in_subset;
}

int
PerMoleculeData::Seen(Molecule& m) {
  const auto iter = _seen.find(m.unique_smiles());
  // Has been seen before
  if (iter != _seen.end()) {
    return 1;
  }

  // Never seen before.
  _seen.emplace(m.unique_smiles());

  // Do not output unique smiles.
  m.invalidate_smiles();

  return 0;
}

ScaffoldFinder::ScaffoldFinder() {
  _max_systems_lost = 0;
  _min_systems_in_subset = 0;

  return;
}

int
ScaffoldFinder::Initialise(Command_Line& cl, char flag) {
  if (! cl.option_present(flag)) {
    return 1;
  }

  IWString fname = cl.string_value(flag);

  std::optional<Scaffolds::ScaffoldsOptions> maybe_options =
       iwmisc::ReadTextProto<Scaffolds::ScaffoldsOptions>(fname);
   if (! maybe_options) {
     cerr << "ScaffoldFinder::Initialise:cannot read proto '" << fname << "'\n";
     return 0;
   }

  return Initialise(*maybe_options);
}

int 
ScaffoldFinder::Initialise(const Scaffolds::ScaffoldsOptions& proto) {
  _config = proto;
  cerr << _config.ShortDebugString() << '\n';

  return 1;
}

// Return true if a molecule that starts with `in_molecule` ring systems,
// and a subset that contains `in_subset` ring systems is OK wrt any
// constraints on ring systems lost or a minimum count.
int
ScaffoldFinder::OkSubsetSize(uint32_t in_molecule, uint32_t in_subset) const {
  if (_config.max_systems_lost() > 0 &&
      (in_molecule - in_subset) > _config.max_systems_lost()) {
    return 0;
  }

  if (_config.min_systems_in_subset() > 0 && in_subset < _config.min_systems_in_subset()) {
    return 0;
  }

  return 1;
}

int
ScaffoldFinder::MakeScaffolds(Molecule& m,
                              Scaffolds::ScaffoldData& result) {
  const int nr = m.nrings();
  if (nr == 0) {
    return 1;
  }

  std::unique_ptr<int[]> in_sys = std::make_unique<int[]>(m.natoms());

  // Strip `m` to the scaffold.

  m.identify_spinach(in_sys.get());
#ifdef DEBUG_SCAFFOLD_FINDER
  for (int i = 0; i < m.natoms(); ++i) {
    cerr << "Atom " << i << " scafold " << in_sys[i] << '\n';
  }
#endif

  Molecule subset;
  m.create_subset(subset, in_sys.get(), 0);

  subset.set_name(m.name());
  // cerr << "Subset is " << subset.smiles() << '\n';

  int rc = Process(subset, result);

  ++_generated[result.subset_size()];

  return rc;
}

// Recursively visit atoms labeling `region` and populating
// systems_reachable.
int
LabelLinkerRegion(const Molecule& m,
             int* region,
             atom_number_t zatom,
             const int* ring_sys,
             int label,
             resizable_array<int>& systems_reachable) {
  region[zatom] = label;
  //cerr << "LabelLinkerRegion atom " << zatom << ' ' << m.smarts_equivalent_for_atom(zatom) << '\n';

  const Atom& a = m.atom(zatom);
  for (const Bond* b : a) {
    atom_number_t j = b->other(zatom);

    // Already assigned
    if (region[j] >= 0) {
      continue;
    }

    //cerr << "From atom " << zatom << " to atom " << j << " rs " << ring_sys[j] << '\n';
    if (ring_sys[j]) {
      systems_reachable.add_if_not_already_present(ring_sys[j]);
      continue;
    }

    LabelLinkerRegion(m, region, j, ring_sys, label, systems_reachable);
  }

  return 1;
}

// For each non-ring atom, assign a value to `region`. Atoms that are
// in the same inter-ring region will have the same region number.
// At the same time, fill `systems_reachable`. For each region, it
// is a list of the ring system identifiers that connect to a 
// region number.
int
LabelLinkerRegions(Molecule& m,
             PerMoleculeData& pmd,
             resizable_array<int>* systems_reachable) {
  const int matoms = m.natoms();

  int* region = pmd.region();
  std::fill_n(region, matoms, -1);

  // cerr << "LabelLinkerRegions processing " << m.smiles() << '\n';

  int label = 0;

  const int* ring_sys = pmd.ring_sys();

  for (int i = 0; i < matoms; ++i) {
    if (ring_sys[i]) {
      continue;
    }
    if (region[i] >= 0) {
      continue;
    }

    LabelLinkerRegion(m, region, i, ring_sys, label, systems_reachable[label]);

    ++label;
  }

  return 1;
}

// #define DEBUG_SCAFFOLD_FINDER

// `m` is the scaffold of an incoming molecule.
int
ScaffoldFinder::Process(Molecule& m,
                        Scaffolds::ScaffoldData& result) {
  PerMoleculeData pmd(m);

  uint32_t nsys = pmd.nsys();
#ifdef DEBUG_SCAFFOLD_FINDER
  cerr << "Scaffold molecule " << m.smiles() << " has " << m.nrings() << " ringa and " << nsys << " ring systems\n";
#endif

  result.set_ring_sys(nsys);

  ++_nsys[nsys];

  // If just 1 ring system in the molecule, process it directly.
  if (nsys == 1) {
    Scaffolds::ScaffoldSubset* s = result.add_subset();
    s->set_smi(m.smiles().AsString());
    s->set_ring_sys(1);
    return 1;
  }

  // Multiple ring systems present.

#ifdef DEBUG_SCAFFOLD_FINDER
  for (int i = 0; i < m.natoms(); ++i) {
    cerr << i << ' ' << m.smarts_equivalent_for_atom(i) << " sys " << pmd.ring_sys(i) << '\n';
  }
#endif

  // For each inter-ring region, a list of the ring system numbers that are adjacent.
  // There will always be < nsys such regions.
  std::unique_ptr<resizable_array<int>[]> systems_reachable = 
                std::make_unique<resizable_array<int>[]>(nsys);

  LabelLinkerRegions(m, pmd, systems_reachable.get());

#ifdef DEBUG_SCAFFOLD_FINDER
  cerr << "Region assignments\n";
  for (int q = 0; q < matoms; ++q) {
    cerr << " atom " << q << ' ' << m.smarts_equivalent_for_atom(q) << " region " << pmd.region()[q] << '\n';
  }
  for (int q = 0; q < nsys; ++q) {
    cerr << "region " << q << " attached to";
    for (int s : systems_reachable[q]) {
      cerr << ' ' << s;
    }
    cerr << '\n';
  }
#endif

  // Set up arrays needed for combinations.
  std::vector<int> count(nsys);
  std::fill(count.begin(), count.end(), 2);
  std::vector<int> state(nsys);

  combinations::Combinations comb(count);

  while (comb.Next(state)) {
    const uint32_t nset = std::count(state.begin(), state.end(), 1);
#ifdef DEBUG_SCAFFOLD_FINDER
    for (int a: state) {
      cerr << ' ' << a;
    }
    cerr << " nset " << nset << '\n';
#endif
    if (nset == 0) {
      continue;
    }

    if (! OkSubsetSize(pmd.nsys(), nset)) {
      continue;
    }

    if (nset < _config.min_systems_in_subset()) {
      continue;
    }
    // cerr << "Will be processed\n";

    Process(m, state, pmd, systems_reachable.get(), result);
  }

  return 1;
}

// Two ring systems, `sys1` and `sys2` are part of a subset. For each region
// that is joined to those two ring systems, add all the atoms in those regions
// to pmd.atoms_in_subset.
void
AddInterRingAtoms(const Molecule& m,
                  PerMoleculeData& pmd,
                  const resizable_array<int>* systems_reachable,
                  int sys1,
                  int sys2) {
  const int matoms = m.natoms();

  const int* region = pmd.region();
  int* atoms_in_subset = pmd.atoms_in_subset();
  // cerr << "Filling in between " << sys1 << " and " << sys2 << '\n';

  // For each region (which is < pmd.nsys) check to see if both of the ring
  // systems are accessible from that region.
  for (int i = 0; i < pmd.nsys(); ++i) {
    if (! systems_reachable[i].contains(sys1) ||
        ! systems_reachable[i].contains(sys2)) {
      continue;
    }
    // cerr << "Region " << i << " is OK\n";

    // Region `i` is joined to both ring systems. Include all those atoms.
    for (int j = 0; j < matoms; ++j) {
      if (region[j] == i) {
        atoms_in_subset[j] = 1;
      }
    }
  }
}

// Process a subset of the ring systems, as in `state`.
// `region` holds, for each non ring atom, a label that
// indicates a particular set of inter-ring atoms.
// `systems_reachable` holds, for each region, a list of
// the ring systems that can be reached from atoms in that region.
int
ScaffoldFinder::Process(const Molecule& m,
                 const std::vector<int>& state,
                 PerMoleculeData& pmd,
                 const resizable_array<int>* systems_reachable,
                 Scaffolds::ScaffoldData& result) {
  assert(state.size() > 1);

  const int matoms = m.natoms();

  const int* ring_sys = pmd.ring_sys();

  int* atoms_in_subset = pmd.atoms_in_subset();
  std::fill_n(atoms_in_subset, matoms, 0);

  // First add the atoms in the ring systems present.
  // `state` is numbered 0,1,2, but ring_sys starts at 1 for atoms that are in a ring.
#ifdef DEBUG_SCAFFOLD_FINDER
  for (int i = 0; i < m.natoms(); ++i) {
    cerr << " atom " << i << ' ' << m.smarts_equivalent_for_atom(i) << " ring_sys " << ring_sys[i] << " region " << pmd.region()[i] << '\n';
  }
#endif

  const int n = state.size();
  int nsys = 0;
  for (int i = 0; i < n; ++i) {
    if (state[i] == 0) {
      continue;
    }

    ++nsys;
    for (int j = 0; j < matoms; ++j) {
      if (ring_sys[j] == i + 1) {
        atoms_in_subset[j] = 1;
      }
    }
  }

#ifdef DEBUG_SCAFFOLD_FINDER
  cerr << "At end of ring system assignments\n";
  for (int q = 0; q < m.natoms(); ++q) {
    cerr << q << " atoms_in_subset " << atoms_in_subset[q] << '\n';
  }
#endif

  // Now each region that is between any pair of ring systems included.
  for (int i = 0; i < n; ++i) {
    if (state[i] == 0) {
      continue;
    }
    for (int j = i + 1; j < n; ++j) {
      if (state[j] == 0) {
        continue;
      }

      AddInterRingAtoms(m, pmd, systems_reachable, i + 1, j + 1);
    }
  }

#ifdef DEBUG_SCAFFOLD_FINDER
  cerr << "Subset atoms ";
  for (int i = 0; i < matoms; ++i) {
    if (atoms_in_subset[i]) {
      cerr << ' ' << i;
    }
  }
  cerr << '\n';
#endif

  Molecule subset;
  m.create_subset(subset, atoms_in_subset, 1);
  // cerr << "Subset " << subset->smiles() << '\n';
  if (subset.number_fragments() > 1) {
    return 0;
  }

  return AddToResultsIfUnique(subset, pmd, nsys, result);
}

int
ScaffoldFinder::AddToResultsIfUnique(Molecule& candidate,
                PerMoleculeData& pmd,
                int nsys,
                Scaffolds::ScaffoldData& result) {
  if (pmd.Seen(candidate)) {
    return 0;
  }

  Scaffolds::ScaffoldSubset* subset = result.add_subset();
  subset->set_smi(candidate.smiles().AsString());
  subset->set_ring_sys(nsys);

  // cerr << "Added unique " << candidate.smiles() << " result now contains " << result.subset_size() << " items\n";

  return 1;
}

int
ScaffoldFinder::Report(std::ostream& output) const {
  for (int i = 0; i < _nsys.number_elements(); ++i) {
    if (_nsys[i]) {
      output << _nsys[i] << " molecules had " << i << " ring systems\n";
    }
  }

  for (int i = 0; i < _generated.number_elements(); ++i) {
    if (_generated[i]) {
      output << _generated[i] << " molecules generated " << i << " scaffold subsets\n";
    }
  }

  return 1;
}

}  // namespace scaffolds
