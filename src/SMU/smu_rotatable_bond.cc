//Examine SMU 3d structures to figure out which bonds are rotatable.

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <tuple>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"

namespace smu_rotatable_bonds {

using std::cerr;

void
Usage(int rc) {
  cerr << "Identifies variability in torsion angles across SMU molecules\n";
  cerr << " -h            remove all explicit Hydrogen atoms\n";
  cerr << " -p <n>        support level for writing a torsion\n";
  cerr << " -g ...        chemical standardisation\n";
  cerr << " -A ...        aromaticity\n";
  cerr << " -i ...        input type specification\n";
  cerr << " -v            verbose output\n";
  exit(rc);
}

// A smarts, which will include all atoms attached to
// a bond, and the identify of the query atoms that are
// the canonical form.
struct SmartsAndDihedral {
  IWString smarts;
  Set_of_Atoms canonical_atoms;
};


struct Torsion {
  // For each example of this torsion, an accumulation of
  // the observed ranges.
  // For each example of this torsion, an accumulation of
  // the observed angles.
  Accumulator<double> acc_angle;
  // Ultimately we want to compute the variance, and that is complicated by
  // the fact that angles can be either side of 0,180 and be nearly the same.
  // Therefore use two accumulators, and the variance will be the minimum of
  // those.
  Accumulator<double> acc_angle_pi2;

  // The first molecule that exemplifies this torsion.
  IWString exemplar;

  // Smarts derived from that first molecule.
  // Note that subsequent molecules will not necessarily
  // match this query due to differences in ring membership, etc.
  IWString first_smarts;

  // Discretized to pi/20
  extending_resizable_array<int> all_angles;
};

struct JobParameters {
  // When comparing dihedral angles, if the difference across conformers
  // is less than `min_angle`, consider that non rotatable.

  float min_angle = 0.00;  // Radians

  int molecules_read = 0;

  extending_resizable_array<uint64_t> acc_conformers;

  IW_STL_Hash_Map<IWString, Torsion> acc_angle;

  // Do we process all molecules, or only those with multiple
  // conformers.
  bool process_as_groups = false;

  bool consider_all_attachments = true;

  // Should probably always be true.
  bool remove_hydrogens = false;

  Chemical_Standardisation chemical_standardisation;

  int verbose = 0;
};

IWString
MakeExemplar(Molecule& m,
             const Set_of_Atoms& atoms) {
  Molecule mcopy(m);
  mcopy.reset_all_atom_map_numbers(); 
  for (int i = 0; i < atoms.number_elements(); ++i) {
    mcopy.set_atom_map_number(atoms[i], i + 1);
  }

  return mcopy.smiles();
}

IWString
MakeExemplarIsotopic(Molecule& m,
                     const Set_of_Atoms& atoms) {
  Molecule mcopy(m);
  mcopy.transform_to_non_isotopic_form();
  mcopy.reset_all_atom_map_numbers();

  for (int i = 0; i < atoms.number_elements(); ++i) {
    mcopy.set_isotope(atoms[i], i + 1);
  }

  return mcopy.smiles();
}

// Return true if `b` looks like a possibly rotatable bond.
bool
BondOkToProcess(Molecule& m,
                const Bond & b) {
  if (! b.is_single_bond()) {
    return false;
  }

  const atom_number_t a2 = b.a1();
  const atom_number_t a3 = b.a2();

  // If terminal, not rotatable.
  if (m.ncon(a2) == 1 || m.ncon(a3) == 1) {
    return false;
  }

  if (m.ring_bond_count(a2) > 0 || m.ring_bond_count(a3) > 0) {
    return 0;
  }

  m.ring_membership();

  if (b.nrings()) {
    return false;
  }

  // disallow any unsaturation.
  if (m.ncon(a2) < m.nbonds(a2)) {
    return false;
  }
  if (m.ncon(a3) < m.nbonds(a3)) {
    return false;
  }

  return true;
}

void
AppendBondSymbol(const Bond& b,
                 IWString& destination) {
  if (b.is_single_bond()) {
    return;
  }

  if (b.is_double_bond()) {
    destination << '=';
  } else if (b.is_triple_bond()) {
    destination << '#';
  }
}

// Return a smarts for `atoms` in `m`.
// `atoms` must be a linear, bonded path.
IWString
MakeSmarts(Molecule& m,
           const Set_of_Atoms& atoms) {
  IWString result;

  atom_number_t previous_atom = INVALID_ATOM_NUMBER;
  for (auto current_atom : atoms) {
    if (previous_atom != INVALID_ATOM_NUMBER) {
      const Bond * b = m.bond_between_atoms(previous_atom, current_atom);
      AppendBondSymbol(*b, result);
    }
    result << m.smarts_equivalent_for_atom(current_atom);
    previous_atom = current_atom;
  }
  return result;
}

//#define DEBUG_ASSEMBLE_DIHEDRAL
// Generate a cononical set of atoms describing the bond `b`.
// unique_smiles will be set to the unique smiles of the atoms
// surrounding `b`.
// `dihedral` will be set to 4 canonical atoms that define the
// angle.
// THe only reason `m` is not const is that the smarts atom
// formation is non const.
int
AssembleAtoms(Molecule& m,
              const Bond& central_bond,
              IWString& unique_smiles,
              Set_of_Atoms& dihedral) {
  // Assemble all the atoms connected to `central_bond`, including the two central atoms.
  Set_of_Atoms subset_atoms;

  const Atom& a1 = m.atom(central_bond.a1());
  for (const Bond * b : a1) {
    if (b->is_triple_bond()) {
      return 0;
    }
    subset_atoms << b->other(central_bond.a1());
  }

  const Atom& a2 = m.atom(central_bond.a2());
  for (const auto * b : a2) {
    if (b->is_triple_bond()) {
      return 0;
    }
    const atom_number_t j = b->other(central_bond.a2());
    subset_atoms << j;
  }

#ifdef DEBUG_ASSEMBLE_DIHEDRAL
  cerr << "Around bond " << central_bond << " atoms " << subset_atoms << '\n';
  cerr << m.smiles() << '\n';
#endif

  const int matoms = m.natoms();

  std::unique_ptr<int[]> atom_in_subset(new_int(matoms));
  subset_atoms.set_vector(atom_in_subset.get(), 1);

  Molecule subset;
  std::unique_ptr<int[]> xref(new_int(matoms));
  if (! m.create_subset(subset, atom_in_subset.get(), 1, xref.get())) {
    return 0;
  }

  // xref is a cross reference from the original molecule to `subset`.
  // We need the inverse. Re-use the atom_in_subset array.

  std::unique_ptr<int[]> subset_to_parent(new_int(subset.natoms()));

  for (int i = 0; i < matoms; ++i) {
    if (xref[i] >= 0) {  // Atom is part of the subset.
      subset_to_parent[xref[i]] = i;
    }
  }

  const int * canonical_order = subset.canonical_ranks();

  // The unique 4 atoms from subset_atoms that will define the angle.
  // Remember, `dihedral` contains atom numbers in `m`.
  dihedral.resize_keep_storage(0);
  dihedral.extend(4, INVALID_ATOM_NUMBER);

  // The atom numbers of the central bond in the subset.
  const atom_number_t a1_in_subset = xref[central_bond.a1()];
  const atom_number_t a2_in_subset = xref[central_bond.a2()];

  int direction;
  if (canonical_order[a1_in_subset] < canonical_order[a2_in_subset]) {
    dihedral[1] = central_bond.a1();
    dihedral[2] = central_bond.a2();
    direction = 1;
  } else {
    dihedral[1] = central_bond.a2();
    dihedral[2] = central_bond.a1();
    direction = -1;
  }

  int highest_rank_bonded_to_a1 = -1;
  atom_number_t highest_ranked_atom_bonded_to_a1 = INVALID_ATOM_NUMBER;
  int highest_rank_bonded_to_a2 = -1;
  atom_number_t highest_ranked_atom_bonded_to_a2 = INVALID_ATOM_NUMBER;

#ifdef DEBUG_ASSEMBLE_DIHEDRAL
  cerr << "Subset atoms " << subset_atoms << " bond " << central_bond << '\n';
#endif
  const int atoms_in_subset = subset.natoms();
  for (int atom = 0; atom < atoms_in_subset; ++atom) {
    if (atom == a1_in_subset || atom == a2_in_subset) {
      continue;
    }

    const atom_number_t in_parent = subset_to_parent[atom];

    const int corder = canonical_order[atom];
    if (m.are_bonded(central_bond.a1(), in_parent)) {
      if (corder > highest_rank_bonded_to_a1) {
        highest_rank_bonded_to_a1 = corder;
        highest_ranked_atom_bonded_to_a1 = atom;
      }
    } else {  // in_parent is bonded to central_bond.a2()
      if (corder > highest_rank_bonded_to_a2) {
        highest_rank_bonded_to_a2 = corder;
        highest_ranked_atom_bonded_to_a2 = atom;
      }
    }
  }

  if (direction == 1) {
    dihedral[0] = subset_to_parent[highest_ranked_atom_bonded_to_a1];
    dihedral[3] = subset_to_parent[highest_ranked_atom_bonded_to_a2];
  } else {
    dihedral[0] = subset_to_parent[highest_ranked_atom_bonded_to_a2];
    dihedral[3] = subset_to_parent[highest_ranked_atom_bonded_to_a1];
  }

#ifdef DEBUG_ASSEMBLE_DIHEDRAL
  cerr << "dihedral " <<  dihedral << '\n';
#endif

  if (! m.are_bonded(dihedral[0], dihedral[1])) {
    cerr << "Atoms 0 1 not bonded " << dihedral << '\n';
  }

  if (! m.are_bonded(dihedral[1], dihedral[2])) {
    cerr << "Atoms 1 2 not bonded " << dihedral << '\n';
  }
  if (! m.are_bonded(dihedral[2], dihedral[3])) {
    cerr << "Atoms 2 3 not bonded " << dihedral << '\n';
  }

  unique_smiles = MakeSmarts(m, dihedral);

  return 1;
}

constexpr float float_pi = static_cast<float>(M_PI);

int
AngleToIndex(float angle) {
  constexpr float multiplier = 1.0f / static_cast<float>(M_PI) * 20.0f;

  return static_cast<int>(angle * multiplier);
}

double
MinVaranmce(const Torsion& torsion) {
  if (torsion.acc_angle.n() < 2) {
    return 0.0;
  }
  double v1 = torsion.acc_angle.variance();
  double v2 = torsion.acc_angle_pi2.variance();

  return std::min(v1, v2);
}

int
CheckVariability(Molecule& m,
                 const Bond& bond,
                 const Set_of_Atoms& atoms,
                 IWString& usmi,
                 const SmartsAndDihedral& smarts,
                 JobParameters & params) {
  auto angle = m.dihedral_angle(atoms[0], atoms[1], atoms[2], atoms[3]);

  auto iter = params.acc_angle.find(usmi);
  if (iter != params.acc_angle.end()) {
    iter->second.acc_angle.extra(angle);
    iter->second.all_angles[AngleToIndex(angle)]++;
    return 1;
  }

  Torsion torsion;
  torsion.acc_angle.extra(angle);
  torsion.all_angles[AngleToIndex(angle)]++;
  torsion.exemplar = MakeExemplarIsotopic(m, atoms);
  if (smarts.smarts.empty()) {
    torsion.first_smarts = std::move(MakeSmarts(m, atoms));
  } else {
    torsion.first_smarts = smarts.smarts;
  }
  params.acc_angle.emplace(std::move(usmi), torsion);  // OK to destroy usmi

  return 1;
}

int
CheckVariability(const resizable_array_p<Molecule>& mols,
                 const Bond& bond,
                 const Set_of_Atoms& atoms,
                 IWString& usmi,
                 JobParameters & params) {
  Accumulator<double> acc_angle;

  constexpr float float_pi = static_cast<float>(M_PI);

  for (const Molecule * m : mols) {
    auto angle = m->dihedral_angle(atoms[0], atoms[1], atoms[2], atoms[3]);
    if (angle < 0.0f)  {
      angle = - angle;
    }
    if (angle >= float_pi) {
      acc_angle.extra(angle - float_pi);
    } else {
      acc_angle.extra(angle);
    }
    if (angle < 0.0f)  {
      cerr << "NEgative angle " << angle << '\n';
      angle = 0.0f;
    }
  }

  auto iter = params.acc_angle.find(usmi);
  if (iter != params.acc_angle.end()) {
    iter->second.acc_angle.extra(acc_angle);
    return 1;
  }

  Torsion torsion;
  torsion.acc_angle.extra(acc_angle);
  torsion.exemplar = MakeExemplarIsotopic(*mols[0], atoms);
  torsion.first_smarts = std::move(MakeSmarts(*mols[0], atoms));
  params.acc_angle.emplace(std::move(usmi), torsion);  // OK to destroy usmi

  return 1;
}

int
ProcessConformerVariants(const resizable_array_p<Molecule>& mols,
                         JobParameters& params) {
  const int n = mols.number_elements();
  params.acc_conformers[n]++;
  if (n == 0) {
    return 1;
  }
  if (n == 1) {
    return 1;
  }

  for (int i = 0; i < n; ++i) {
    mols[i]->ring_membership();
  }

  const int matoms = mols[0]->natoms();

#ifdef CHECK_INTERASTOMIC_DISTANCES
  Accumulator<float> acc_range;
  for (int i = 0; i < matoms; ++i) {
    for (int j = i + 1; j < matoms; ++j) {
      Accumulator<float> acc_dist;
      for (int k = 0; k < n; ++k) {
        const auto dist = mols[k]->distance_between_atoms(i, j);
        acc_dist.extra(dist);
      }
      acc_range.extra(acc_dist.range());
    }
  }
  cerr << mols[0]->name() << " interatomic distances range " << acc_range.range() << '\n';
#endif

  for (int i = 1; i < n; ++i) {
    if (mols[i]->natoms() != matoms) {
      cerr << "Atom count mismatch " << mols[i]->natoms() << " vs " << matoms << '\n';
      return 0;
    }

    for (int j = 0; j < matoms; ++j) {
      if (mols[0]->atomic_number(j) != mols[i]->atomic_number(j)) {
        cerr << "Atomic number mismatch\n";
        return 0;
      }
    }
  }

  const int nedges = mols[0]->nedges();
  for (int i = 0; i < nedges; ++i) {
    const Bond * b = mols[0]->bondi(i);
    if (! BondOkToProcess(*mols[0], *b)) {
      continue;
    }

    IWString usmi;
    Set_of_Atoms dihedral;
    if (! AssembleAtoms(*mols[0], *b, usmi, dihedral)) {
      continue;
    }

    // usmi may be destroyed, ok.
    CheckVariability(mols, *b, dihedral, usmi, params);
  }

  return 1;
}

bool
SameId(const IWString& stripped,
       const IWString& unstripped) {
  if (! unstripped.starts_with(stripped)) {
    return 0;
  }

  return stripped.length() + 3 == unstripped.length();
}

IWString
IDFromConformerNumber(const IWString& id) {
  IWString result(id);
  result.chop(3);
  return result;
}

void
Preprocess(JobParameters& params,
           Molecule& m) {
  m.reset_all_atom_map_numbers();

  if (params.remove_hydrogens) {
    m.remove_all(1);
  }

  if (params.chemical_standardisation.active()) {
    params.chemical_standardisation.process(m);
  }
}

int
SmuRotatableBonds(Molecule& m,
                  JobParameters& params) {
  Preprocess(params, m);

  (void) m.ring_membership();

  SmartsAndDihedral smarts;  // Not used here.

  const int nedges = m.nedges();
  for (int i = 0; i < nedges; ++i) {
    const Bond * b = m.bondi(i);
    if (! BondOkToProcess(m, *b)) {
      continue;
    }

    IWString usmi;
    Set_of_Atoms dihedral;
    if (! AssembleAtoms(m, *b, usmi, dihedral)) {
      continue;
    }
    CheckVariability(m, *b, dihedral, usmi, smarts, params);
  }

  return 1;
}

// Instead of only processing molecules with multiple conformers,
// process all molecules singly.
int
SmuRotatableBondsEach(data_source_and_type<Molecule>& input,
                  JobParameters& params) {
  Molecule * m;
  while (NULL != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);
    params.molecules_read++;
    if (! SmuRotatableBonds(*m, params)) {
      cerr << "Fatal error processing " << m->name() << '\n';
      return 0;
    }
  }

  return 1;
}

int
SmuRotatableBondsGrouped(data_source_and_type<Molecule>& input,
                  JobParameters& params) {
  resizable_array_p<Molecule> mols;
  Molecule * m;
  IWString id;
  while (NULL != (m = input.next_molecule())) {
    Preprocess(params, *m);
    IWString name = m->name();
    name.truncate_at_first(' ');
    m->set_name(name);
//  cerr << "Checking SameId '" << id << "' and '" << m->name() << "' " << SameId(id, m->name()) << '\n';
    if (SameId(id, m->name())) {
      mols << m;
    } else {
      ProcessConformerVariants(mols, params);
      id = IDFromConformerNumber(m->name());
//    cerr << "From '" << m->name() << " id is '" << id << "'\n";
      mols.resize_keep_storage(0);
      mols << m;
      params.molecules_read++;
    }
  }

  ProcessConformerVariants(mols, params);
  params.molecules_read++;

  return 1;
}

int
NumberBonds(const Bond & bond) {
  if (bond.is_aromatic()) {
    return 4;
  }
  return bond.number_of_bonds();
}

uint32_t
ComputeHash(const Molecule& m,
            atom_number_t a1,
            const Bond & b12,
            const atom_number_t a2,
            const atom_number_t a3,
            const Bond & b34,
            const atom_number_t a4) {
  static uint32_t magic[] { 27, 46, 93 };
  return m.atomic_number(a1) + magic[0] * NumberBonds(b12) * m.atomic_number(a2) +
                     magic[1] * m.atomic_number(a3) +
                     magic[2] * NumberBonds(b34) * m.atomic_number(a4);
}

struct AtomsAndCount {
  atom_number_t a1;
  atom_number_t a2;
  atom_number_t a3;
  atom_number_t a4;
  int instances = 1;

  AtomsAndCount() {
    instances = 1;
  }
  AtomsAndCount(atom_number_t q1,
                atom_number_t q2,
                atom_number_t q3,
                atom_number_t q4) : a1(q1), a2(q2), a3(q3), a4(q4) {
    instances = 1;
  }
};

void
AddToHash(uint32_t score, AtomsAndCount&& atoms_and_count,
          std::unordered_map<uint32_t, AtomsAndCount>& scores) {
  auto iter = scores.find(score);
  if (iter != scores.end()) {
    iter->second.instances++;
  } else {
    scores.emplace(std::make_pair(score, std::move(atoms_and_count)));
  }
}

// The objective is to identify a canonical dihedral angle
// with the bond between a2-a3 in the middle.
// If it is possible to identify canonical atoms a1 and a4,
// those are returned. If no canonical choice can be made,
// nullopt is returned.
std::optional<AtomsAndCount>
CanonicalConnections(const Molecule& m,
                     const Bond & b) {
  const atom_number_t a2 = b.a1();
  const atom_number_t a3 = b.a2();
  const Atom& atom2 = m.atom(a2);
  const Atom& atom3 = m.atom(a3);
  // Mapping from score to dihedral atoms.
  std::unordered_map<uint32_t, AtomsAndCount> scores;

  for (const Bond * b2 : atom2) {
    atom_number_t j = b2->other(a2);
    if (j == a3) {
      continue;
    }
    for (const Bond * b3 : atom3) {
      const atom_number_t k = b3->other(a3);
      if (k == a2) {
        continue;
      }
      // Compute forward and reverse hashes.
      uint32_t v = ComputeHash(m, j, *b2, a2, a3, *b3, k);
      AddToHash(v, AtomsAndCount(j, a2, a3, k), scores);
      v = ComputeHash(m, k, *b3, a3, a2, *b2, j);
      AddToHash(v, AtomsAndCount(k, a3, a2, j), scores);
    }
  }

  // Looking for the highest score that has one instance.
  uint32_t highest_score = 0;
  AtomsAndCount canonical_atoms;

  for (const auto& [score, atoms] : scores) {
    if (score < highest_score) {
      continue;
    }
    if (atoms.instances > 1) {
      continue;
    }
    highest_score = score;
    canonical_atoms = atoms;
  }

  if (highest_score == 0) {
    return std::nullopt;
  }

  return canonical_atoms;
}

int
AssembleAtomsAllAttachments(Molecule& m,
                            const AtomsAndCount& ac,
                            IWString& usmi,
                            SmartsAndDihedral& smarts) {
  const atom_number_t a1 = ac.a1;
  const atom_number_t a2 = ac.a2;
  const atom_number_t a3 = ac.a3;
  const atom_number_t a4 = ac.a4;
  Set_of_Atoms connected2 = m.connections(a2);
  Set_of_Atoms connected3 = m.connections(a3);

  const int matoms = m.natoms();

  std::unique_ptr<int[]> atom_in_subset(new_int(matoms));
  connected2.set_vector(atom_in_subset.get(), 1);
  connected3.set_vector(atom_in_subset.get(), 1);

  Molecule subset;
  std::unique_ptr<int[]> xref(new_int(matoms));
  if (! m.create_subset(subset, atom_in_subset.get(), 1, xref.get())) {
    return 0;
  }

  usmi = subset.unique_smiles();

  // By convention, we write the canonical atoms as #atomic_number

  smarts.smarts << "[#" << m.atomic_number(a1) << ']';
  smarts.smarts << m.smarts_equivalent_for_atom(a2);
  for (atom_number_t atom : connected2) {
    if (atom == a3 || atom == a1) {
      continue;
    }
    smarts.smarts << "(-" << m.atomic_symbol(atom) << ')';
  }
  smarts.smarts << '-';
  smarts.smarts << m.smarts_equivalent_for_atom(a3);
  smarts.smarts << "(-[#" << m.atomic_number(a4) << "])";
  for (atom_number_t atom : connected3) {
    if (atom == a2 || atom == a4) {
      continue;
    }
    smarts.smarts << "(-" << m.atomic_symbol(atom) << ')';
  }
  return 1;
}

// Examine the atoms around `b` and if a canonical form can be
// identified, fill in `usmu`, `smarts`, and `dihedral`.
int
AssembleAtomsAllAttachments(Molecule& m,
                            const Bond & b,
                            IWString& usmi,
                            SmartsAndDihedral& smarts,
                            Set_of_Atoms& dihedral) {
  std::optional<AtomsAndCount> canonical_atoms = CanonicalConnections(m, b);
  if (canonical_atoms == std::nullopt) {
    return 0;
  }

  const AtomsAndCount& ac = canonical_atoms.value();
  dihedral << ac.a1 << ac.a2 << ac.a3 << ac.a4;
  return AssembleAtomsAllAttachments(m, ac, usmi, smarts);
}

int
ExtractTorsionsAllAttachments(Molecule& m,
                              JobParameters& params) {
//cerr << "Processing " << m.name() << '\n';
  const int nedges = m.nedges();
  for (int i = 0; i < nedges; ++i) {
    const Bond * b = m.bondi(i);
    if (! BondOkToProcess(m, *b)) {
//    cerr << "cannot process bond\n";
      continue;
    }

    IWString usmi;
    SmartsAndDihedral smarts;
    Set_of_Atoms dihedral;
    if (! AssembleAtomsAllAttachments(m, *b, usmi, smarts, dihedral)) {
//    cerr << "AssembleAtomsAllAttachments failed\n";
      continue;
    }
//  cerr << "CheckVariability " << usmi << '\n';
    CheckVariability(m, *b, dihedral, usmi, smarts, params);
  }

  return 1;
}

int
SmuRotatableBondsAllAttachments(data_source_and_type<Molecule>& input,
                        JobParameters& params)  {
  Molecule * m;
  while (NULL != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);
    Preprocess(params, *m);
    ExtractTorsionsAllAttachments(*m, params);
    params.molecules_read++;
  }

  return 1;
}

int
SmuRotatableBonds(const char * fname,
                  JobParameters& params) {
  data_source_and_type<Molecule> input(FILE_TYPE_SMI, fname);
  if (! input.good())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 1;
  }

  if (params.consider_all_attachments) {
    return SmuRotatableBondsAllAttachments(input, params);
  } else if (params.process_as_groups) {
    return SmuRotatableBondsGrouped(input, params);
  } else {
    return SmuRotatableBondsEach(input, params);
  }
}

int
SmuRotatableBonds(int argc, char ** argv) {
  Command_Line cl(argc, argv, "vi:A:hp:g:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  FileType input_type = FILE_TYPE_SMI;

  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      return 1;
    }
  }

  JobParameters params;

  params.verbose = cl.option_count('v');

  if (cl.option_present('h')) {
    params.remove_hydrogens = true;
    if (params.verbose) {
      cerr << "Will remove all explicit Hydrogen atoms\n";
    }
  }

  if (cl.option_present('g')) {
    if (! params.chemical_standardisation.construct_from_command_line(cl, params.verbose, 'g')) {
      cerr << "Invalid chemical standardisation specification (-g)\n";
      return 1;
    }
  }

  unsigned int min_count_for_output = 0;
  if (cl.option_present('p')) {
    if (! cl.value('p', min_count_for_output) || min_count_for_output < 0) {
      cerr << "Invalid -p options\n";
      return 1;
    }
    if (params.verbose) {
      cerr << "Will only display torsions with " << min_count_for_output << " exemplars\n";
    }
  }

  if (cl.empty()) {
    cerr << "INsufficient arguments\n";
    Usage(1);
  }

  set_make_smarts_embedding(0);

  for (const char * fname : cl) {
    if (! SmuRotatableBonds(fname, params)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  cerr << params.molecules_read << " molecules read\n";
  for (int i = 0; i < params.acc_conformers; ++i) {
    if (params.acc_conformers[i] > 0) {
      cerr << params.acc_conformers[i] << " molecules had " << i << " conformers\n";
    }
  }

  const char sep = ',';
  // Header record
  std::cout << "smiles" << sep << "num.smi" << sep << "smarts" << sep <<
            "N" << sep <<
            "min" << sep <<
            "ave" << sep <<
            "max" << sep <<
            "range" << sep <<
            "var" << sep <<
            "nbuckets" << sep <<
            "highest_occupancy" << sep <<
            "highest_occupancyf" << sep << '\n';

  for (const auto& [smiles, torsion] : params.acc_angle) {
    const auto& acc = torsion.acc_angle;
    if (static_cast<float>(acc.minval()) < 0.0) {
      cerr << "negative result " << acc.minval() << '\n';
    }
    if (torsion.acc_angle.n() < min_count_for_output) {
      continue;
    }
    int buckets_occupied = 0;
    int highest_occupancy = 0;
    for (int b : torsion.all_angles) {
      if (b > 0) {
        buckets_occupied++;
        if (b > highest_occupancy) {
          highest_occupancy = b;
        }
      }
    }
    std::cout << std::setprecision(3);
    std::cout << torsion.exemplar << sep << smiles << sep <<
              torsion.first_smarts << sep <<
              acc.n() << sep <<
              static_cast<float>(acc.minval()) << sep <<
              static_cast<float>(acc.average()) << sep <<
              static_cast<float>(acc.maxval()) << sep <<
              static_cast<float>(acc.range()) << sep << 
              static_cast<float>(sqrt(MinVaranmce(torsion))) << sep << 
              buckets_occupied << sep << highest_occupancy << sep <<
              static_cast<float>(highest_occupancy) / static_cast<float>(acc.n()) << '\n';

  }

  return 0;
}

}  // namespace smu_rotatable_bonds

int
main(int argc, char ** argv) {
  return smu_rotatable_bonds::SmuRotatableBonds(argc, argv);
}
