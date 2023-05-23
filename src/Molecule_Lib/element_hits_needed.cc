#include <stdlib.h>

#include "misc2.h"
#include "substructure.h"
#include "target.h"

Elements_Needed::Elements_Needed()
{
  _z = INVALID_ATOMIC_NUMBER;

  return;
}

Elements_Needed::Elements_Needed(atomic_number_t s) : _z(s)
{
}

int
Elements_Needed::ok() const
{
  return 1;
}

int
Elements_Needed::debug_print(std::ostream & os, const IWString & ind) const
{
  os << ind << "Elements_Needed: z = " << _z << " hits needed ";
  return _hits_needed.debug_print(os);
}

int
Elements_Needed::matches(Query_Atoms_Matched & qam) const
{
  if (INVALID_ATOMIC_NUMBER == _z)
  {
    cerr << "Elements_Needed::matches: cannot match invalid atom number\n";
    debug_print(cerr, "");
    iwabort();
    return 0;
  }

  int esize = qam.number_elements();

  int hits = 0;
  for (int i = 0; i < esize; i++)
  {
    const Substructure_Atom * a = qam[i];

    if (! a->include_in_embedding())
      continue;

    if (_z == a->current_hold_atom()->atomic_number())
      hits++;
  }

#ifdef DEBUG_ELEMENTS_NEEDED_MATCHES
  cerr << "Searched through " << esize << " matched atoms, found " << hits << " matching atomic number, returning " << Min_Max_Specifier<int>::matches(hits) << endl;
#endif

  return _hits_needed.matches(hits);
}

int
Elements_Needed::matches(Molecule_to_Match & target_molecule) const
{
  int count = target_molecule.atoms_with_atomic_number(_z);

#ifdef DEBUG_ELEMENTS_NEEDED_MATCHES
  cerr << "Elements_Needed::matches: target contains " << count << " returning " << Min_Max_Specifier<int>::matches(count) << endl;
#endif

  return _hits_needed.matches(count);
}

RequiredBond::RequiredBond() {
  _atomic_number_1 = -1;
  _atomic_number_2 = -1;
  _btype = INVALID_BOND_TYPE;
  _min_count = 1;
}

int
RequiredBond::Matches(const Molecule& m) const {
  if (_btype == SINGLE_BOND) {
    return MatchesSingle(m);
  }
  if (_btype == DOUBLE_BOND) {
    return MatchesDouble(m);
  }
  if (_btype == TRIPLE_BOND) {
    return MatchesTriple(m);
  }
  cerr << "RequiredBonds::Matches1:What kind of bond is " << _btype << '\n';
  return 0;
}

int
RequiredBond::MatchesSingle(const Molecule& m) const {
  int rc = 0;
  for (const Bond* b : m.bond_list()) {
    if (! b->is_single_bond()) {
      continue;
    }
    const atomic_number_t z1 = m.atomic_number(b->a1());
    const atomic_number_t z2 = m.atomic_number(b->a2());
    if (z1 == _atomic_number_1 && z2 == _atomic_number_2) {
      ++rc;
    } else if (z2 == _atomic_number_1 && z1 == _atomic_number_2) {
      ++rc;
    } else {
      continue;
    }

    if (rc >= _min_count) {
      return 1;
    }
  }

  return 0;
}

int
RequiredBond::MatchesDouble(const Molecule& m) const {
  int rc = 0;
  for (const Bond* b : m.bond_list()) {
    if (! b->is_double_bond()) {
      continue;
    }
    const atomic_number_t z1 = m.atomic_number(b->a1());
    const atomic_number_t z2 = m.atomic_number(b->a2());
    if (z1 == _atomic_number_1 && z2 == _atomic_number_2) {
      ++rc;
    } else if (z2 == _atomic_number_1 && z1 == _atomic_number_2) {
      ++rc;
    } else {
      continue;
    }

    if (rc >= _min_count) {
      return 1;
    }
  }

  return 0;
}

int
RequiredBond::MatchesTriple(const Molecule& m) const {
  int rc = 0;
  for (const Bond* b : m.bond_list()) {
    if (! b->is_triple_bond()) {
      continue;
    }
    const atomic_number_t z1 = m.atomic_number(b->a1());
    const atomic_number_t z2 = m.atomic_number(b->a2());
    if (z1 == _atomic_number_1 && z2 == _atomic_number_2) {
      ++rc;
    } else if (z2 == _atomic_number_1 && z1 == _atomic_number_2) {
      ++rc;
    } else {
      continue;
    }

    if (rc >= _min_count) {
      return 1;
    }
  }

  return 0;
}
