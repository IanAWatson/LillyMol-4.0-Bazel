#ifndef MOLECULE_TOOLS_RING_CLOSURE_ENUMERATION_LIB_H
#define MOLECULE_TOOLS_RING_CLOSURE_ENUMERATION_LIB_H

#include "Foundational/iwaray/iwaray.h"

#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"

namespace ring_closure_enumeration {

class BondFormation {
  private:
    // When dealing with multiple reagent sets, each BondFormation needs
    // to know to which component is it to be applied.
    int _component;

    // Each isotope is translated into a bond type and a ring closure number.
    int _isotope;
    int _bond_type;
    int _ring_closure;

    // We can speed things up by storing the bond types as character forms.
    char _as_char[4];

  public:
    BondFormation();

    // Specification looks like 95=91
    // which means that atom with isotope 95 is to be transformed
    // to ring opening/closing 91 with a double bond.
    int Build(const const_IWSubstring& buffer);

    int component() const {
      return _component;
    }
    int isotope() const {
      return _isotope;
    }
    char btype() const {
      return _as_char[_bond_type];
    }
    int ring_closure() const {
      return _ring_closure;
    }
};

class Reagents {
  private:
    resizable_array_p<Molecule> _mols;
    resizable_array_p<IWString> _smiles;

  // private functions
    int Read(data_source_and_type<Molecule>& input);

  public:
    int Read(const char * fname);

    int CreateFragments(const BondFormation& bond_formation);
    int CreateFragments(const resizable_array<BondFormation*>& bond_formations);

    int number_reagents() const {
      return _mols.number_elements();
    }

    const IWString& MolName(int ndx) {
      return _mols[ndx]->name();
    }
    // Primarily used for testing. Note that we assume ownership of `m`.
    int Add(Molecule* m) {
      _mols << m;
      return 1;
    }

    const resizable_array_p<Molecule>& mols() const {
      return _mols;
    }

    const resizable_array_p<IWString>& smiles() const {
      return _smiles;
    }
};

}  // namespace ring_closure_enumeration

#endif  // MOLECULE_TOOLS_RING_CLOSURE_ENUMERATION_LIB_H
