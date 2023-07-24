#include <algorithm>
#include <iostream>

#include "Molecule_Tools/dicer_api.h"

namespace dicer_api {

using std::cerr;

using fixed_bit_vector::FixedBitVector;

CurrentMoleculeData::CurrentMoleculeData(int n) : natoms(n) {
}

int
CurrentMoleculeData::SetMaxBondsBroken(int nb, int natoms) {
  subset.resize(nb);
  fragment_membership.resize(nb);
  for (std::unique_ptr<int[]>& e : subset) {
    e.reset(new int[natoms]);
  }
  for (std::unique_ptr<int[]>& f : fragment_membership) {
    f.reset(new int[natoms]);
  }

  return nb;
}

uint32_t
AtomsInSubset(const int * subset, int matoms) {
  return std::accumulate(subset, subset + matoms, 0);
}

DicerApi::DicerApi() {
  _min_fragment_size = 1;
  _max_fragment_size = 12;
  _max_bonds_to_break = 3;
  _break_amide_bonds = 0;
  _break_acid_bonds = 0;
  _break_cc_bonds = 0;
  _isotope = 1;
  _skip_bad_valence = 1;
}

bool
DicerApi::OkSize(int matoms) const {
  if (matoms < _min_fragment_size) {
    return false;
  }
  if (matoms > _max_fragment_size) {
    return false;
  }

  return true;
}

int
IsCf3(Molecule& m,
      atom_number_t fluorine,
      atom_number_t carbon) {
  atomic_number_t atomic_number_connected = m.atomic_number(fluorine);

  int singly_connected = 0;
  for (const Bond* b : m[carbon]) {
    atom_number_t j = b->other(carbon);
    if (j == fluorine) {
      continue;
    }
    const Atom& aj = m[j];
    if (aj.ncon() == 1 && aj.atomic_number() == atomic_number_connected) {
      ++singly_connected;
    }
  }

  return singly_connected == 2;
}

int
IsCf3(Molecule& m, const Bond& b) {
  const Atom& a1 = m[b.a1()];
  const Atom& a2 = m[b.a2()];

  if (a1.ncon() == 1 && a1.atomic_number() == 9 && a2.ncon() == 4) {
    return IsCf3(m, b.a1(), b.a2());
  }
  if (a2.ncon() == 1 && a2.atomic_number() == 6 && a1.ncon() == 4) {
    return IsCf3(m, b.a2(), b.a1());
  }

  return 0;
}

int
IsAmide(Molecule& m,
        atomic_number_t nitrogen,
        atomic_number_t carbon) {
  if (m.atomic_number(carbon) == 6) {
    if (m.ncon(carbon) != 3) {
      return 0;
    }
  } else if (m.atomic_number(carbon) == 16) {
    if (m.ncon(carbon) != 4) {
      return 0;
    }
  } else {
    return 0;
  }

  int douly_bonded = 0;
  for (const Bond * b : m[carbon]) {
    if (b->is_single_bond()) {
      continue;
    }
    if (b->nrings()) {
      continue;
    }
    atom_number_t j = b->other(carbon);

    atomic_number_t zj = m.atomic_number(j);

    if (zj == 7) {
      return 1;
    }

    if (m.ncon(j) != 1) {
      return 0;
    }
    if (zj == 8 || zj == 16) {
      ++douly_bonded;
    } else {
      return 0;
    }
  }

  return douly_bonded;
}

// `ohsh` is either [OH] or [SH] and `cs` is [C,S].
// Return true if there are =[O,S] attached to `cs`
// We are not terribly careful about this...
bool
IsAcid(Molecule& m,
       atom_number_t ohsh,
       atom_number_t cs) {
  int doubly_bonded_cs = 0;

  for (const Bond* b : m[cs]) {
    atom_number_t j = b->other(cs);
    if (j == ohsh) {
      continue;
    }

    const Atom& aj = m[j];
    if (aj.ncon() != 1) {
      continue;
    }

    if (aj.atomic_number() == 8 || aj.atomic_number() == 16 || aj.atomic_number() == 7) {
      ++doubly_bonded_cs;
    }
  }

  return doubly_bonded_cs > 0;
}

// Does this bond describe something that looks like the single
// bond part of an acid. O-C=O, S-C=O, O-S=O...
bool
IsAcid(Molecule& m,
       const Bond& b) {
  assert(b.is_single_bond());

  const Atom& a1 = m[b.a1()];
  const Atom& a2 = m[b.a2()];

  const atomic_number_t z1 = a1.atomic_number();
  const atomic_number_t z2 = a2.atomic_number();

  if (a1.ncon() == 1 && (z1 == 8 || z1 == 16)) {
    return IsAcid(m, b.a1(), b.a2());
  } else if (a2.ncon() == 1 && (z2 == 8 || z2 == 16)) {
    return IsAcid(m, b.a2(), b.a1());
  }

  return 0;
}

int
IsAmide(Molecule& m, const Bond& b) {
  atomic_number_t z1 = m.atomic_number(b.a1());
  if (z1 == 7) {
    return IsAmide(m, b.a1(), b.a2());
  }
  atomic_number_t z2 = m.atomic_number(b.a2());
  if (z2 == 7) {
    return IsAmide(m, b.a2(), b.a1());
  }

  return 0;
}

bool
CurrentMoleculeData::IsNew(std::unique_ptr<FixedBitVector>& candidate) {
  for (const FixedBitVector* b : found) {
    if (*b == *candidate) {
      return false;
    }
  }

  found << candidate.release();

  return true;
}

int
DicerApi::IdentifyBondsToBreak(Molecule& m,
                CurrentMoleculeData& current_molecule_data) {
  m.compute_aromaticity_if_needed();

  const int matoms = m.natoms();

  current_molecule_data.breakable.reserve(matoms);

  for (const Bond* b : m.bond_list()) {
    // cerr << "btw " << m.smarts_equivalent_for_atom(b->a1()) << " and " << m.smarts_equivalent_for_atom(b->a2()) << '\n';
    if (b->nrings()) {
      continue;
    }

    if (b->is_aromatic() || ! b->is_single_bond()) {
      continue;
    }

    // Allow ring-chain bonds to break.
    if (m.ring_bond_count(b->a1()) || m.ring_bond_count(b->a2())) {
      current_molecule_data.breakable.emplace_back(AtomPair(b->a1(), b->a2()));
      continue;
    }

    // cerr << "Checking IsAmide\n";
    if (! _break_amide_bonds && IsAmide(m, *b)) {
      continue;
    }

    // cerr << "Checking IsAcid\n";
    if (! _break_acid_bonds && IsAcid(m, *b)) {
      continue;
    }

    // cerr << "Checking cf3\n";
    if (IsCf3(m, *b)) {
      continue;
    }

    if (_break_cc_bonds) {
    } else if (m.atomic_number(b->a1()) == 6 && m.saturated(b->a1()) &&
               m.atomic_number(b->a2()) == 6 && m.saturated(b->a2())) {
      continue;
    }

    current_molecule_data.breakable.emplace_back(AtomPair(b->a1(), b->a2()));
  }

  return current_molecule_data.breakable.size();
}

template <typename T>
void
Invert(T* values, int n) {
  for (int i = 0; i < n; ++i) {
    values[i] = !values[i];
  }
}

void
SetBothIsotopes(Molecule& m, atom_number_t a1, atom_number_t a2,
                isotope_t iso) {
  m.set_isotope(a1, iso);
  m.set_isotope(a2, iso);
}

int
DicerApi::Process(Molecule& m,
          Sparse_Fingerprint_Creator& sfc) {
  const int matoms = m.natoms();

  CurrentMoleculeData current_molecule_data(matoms);
  int breakable_bonds = IdentifyBondsToBreak(m, current_molecule_data);
  //cerr << "Find " << breakable_bonds << " breakable_bonds bonds\n";
  if (breakable_bonds == 0) {
    return 0;
  }
  current_molecule_data.SetMaxBondsBroken(_max_bonds_to_break, matoms);

  for (uint32_t i = 0; i < current_molecule_data.breakable.size(); ++i) {
    atom_number_t a1 = current_molecule_data.breakable[i].a1;
    atom_number_t a2 = current_molecule_data.breakable[i].a2;

    SetBothIsotopes(m, a1, a2, _isotope);
    m.remove_bond_between_atoms(a1, a2);
    ProcessFragments(m, current_molecule_data, 0, sfc);
    if (_max_bonds_to_break > 1) {
      Process(m, current_molecule_data, 1, i + 1, sfc);
    }
    m.add_bond(a1, a2, SINGLE_BOND);
    SetBothIsotopes(m, a1, a2, 0);
  }

  return FinalProcessing(m, current_molecule_data, sfc);
}

struct IsoSave {
  Molecule& mol;
  atom_number_t a1;
  atom_number_t a2;
  isotope_t iso1;
  isotope_t iso2;

  IsoSave(Molecule& m, atom_number_t x1, atom_number_t x2,
          isotope_t iso) : mol(m), a1(x1), a2(x2) {
    iso1 = m.isotope(a1);
    iso2 = m.isotope(a2);
    m.set_isotope(a1, iso);
    m.set_isotope(a2, iso);
  }

  ~IsoSave() {
    mol.set_isotope(a1, iso1);
    mol.set_isotope(a2, iso2);
  }
};

int
DicerApi::Process(Molecule& m,
                  CurrentMoleculeData& current_molecule_data,
                  int bonds_already_broken,
                  uint32_t bstart,
                  Sparse_Fingerprint_Creator& sfc) {
  if (bstart >= current_molecule_data.breakable.size()) {
    return 0;
  }

  for (uint32_t i = bstart; i < current_molecule_data.breakable.size(); ++i) {
    atom_number_t a1 = current_molecule_data.breakable[i].a1;
    atom_number_t a2 = current_molecule_data.breakable[i].a2;

    IsoSave iso_save(m, a1, a2, _isotope);
    m.remove_bond_between_atoms(a1, a2);
    ProcessFragments(m, current_molecule_data, 0, sfc);
    if (bonds_already_broken + 1 < _max_bonds_to_break) {
      Process(m, current_molecule_data, bonds_already_broken + 1, i + 1, sfc);
    }
    m.add_bond(a1, a2, SINGLE_BOND);
  }

  return 1;
}

void
DicerApi::ProcessFragments(Molecule& m,
                           CurrentMoleculeData& current_molecule_data,
                           int ndx,
                           Sparse_Fingerprint_Creator& sfc) {
  const int matoms = m.natoms();

  int * fragment_membership = current_molecule_data.fragment_membership[ndx].get();
  int * subset = current_molecule_data.subset[ndx].get();

  int nfrag = m.fragment_membership(fragment_membership);

  for (int frag = 0; frag < nfrag; ++frag) {
    std::unique_ptr<FixedBitVector> bitvector = std::make_unique<FixedBitVector>(matoms);
    for (int j = 0; j < matoms; ++j) {
      if (fragment_membership[j] == frag) {
        subset[j] = 1;
        bitvector->set_bit(j);
      } else {
        subset[j] = 0;
      }
    }
    if (! OkSize(bitvector->nset())) {
      continue;
    }

    if (!current_molecule_data.IsNew(bitvector)) {
      continue;
    }

    ProcessFragment(m, current_molecule_data, subset, sfc);
  }
}

// We have broken a bond involving `zatom`. 
int
DicerApi::ProcessFragment(Molecule& m,
                          CurrentMoleculeData& current_molecule_data,
                          const int* subset,
                          Sparse_Fingerprint_Creator& sfc) {

  Smiles_Information smiles_info;
  if (_skip_bad_valence && ! m.valence_ok()) {
    return 0;
  }

  const IWString& usmi = m.unique_smiles(smiles_info, subset);

  uint32_t bit;

  auto iter_ndx = _frag_to_ndx.find(usmi);

  if (iter_ndx == _frag_to_ndx.end()) {
    bit = _frag_to_ndx.size();
    _frag_to_ndx[usmi] = bit;
    _atoms_in_frag[bit] = AtomsInSubset(subset, m.natoms());
    _ndx_to_frag[bit] = usmi;
    _frag_count[bit] = 1;
  } else {
    bit = iter_ndx->second;
    ++_frag_count[bit];
  }

  sfc.hit_bit(bit);
  if (sfc.IsSet(bit) == 1) {
    std::cout << usmi << ' ' << m.name() << '\n';
  }

  return 1;
}

int
DicerApi::Initialise(Command_Line& cl) {
  return 1;
}

int
DicerApi::Report(std::ostream& output) const {
  return 1;
}

int
DicerApi::FinalProcessing(Molecule& m,
                CurrentMoleculeData& current_molecule_data,
                Sparse_Fingerprint_Creator& sfc) {
  // cerr << "Done " << m.name() << " set " << sfc.nbits() << '\n';
  return 1;
}

const IWString&
DicerApi::Smiles(uint32_t bit) const {
  static IWString empty_string;

  auto iter = _ndx_to_frag.find(bit);
  if (iter == _ndx_to_frag.end()) {
    return empty_string;
  }

  return iter->second;
}

}  // namespace dicer_api
