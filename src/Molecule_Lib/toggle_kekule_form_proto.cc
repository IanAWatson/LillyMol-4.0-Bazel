// Separate file for the Toggle_Kekule_Form class that reads
// a configuration from a proto file. That keeps the protobuf
// dependency from minimal_lillymol

// ConstructFromProto is hidden in the class header file unless this
// variable is set.

#include <iostream>

#define COMPILING_TOGGLE_KEKULE_FORM_PROTO

#include "Molecule_Lib/toggle_kekule_form.h"

using std::cerr;

int
Toggle_Kekule_Form::ConstructFromProto(const ToggleKekuleForm::ToggleKekuleForm& proto)
{
  if (proto.bond().empty()) {
    cerr << "ToggleKekuleForm::ConstructFromProto:no bonds " << proto.ShortDebugString() << '\n';
    return 0;
  }

  for (const auto& bond : proto.bond())
  {
    if (! bond.has_a1() || ! bond.has_a2()) {
      cerr << "Toggle_Kekule_Form::ConstructFromProto:no a1/a2 " << proto.ShortDebugString() << '\n';
      return 0;
    }

    if (bond.a1() == bond.a2()) {
      cerr << "Toggle_Kekule_Form::ConstructFromProto:invalid a1/a2 " << proto.ShortDebugString() << '\n';
      return 0;
    }

    if (! bond.has_btype()) {
      cerr << "Toggle_Kekule_Form::ConstructFromProto:no bond type " << proto.ShortDebugString() << '\n';
      return 0;
    }

    bond_type_t bt;
    switch (bond.btype())
    {
      case SubstructureSearch::SS_SINGLE_BOND:
        bt = SINGLE_BOND;
        break;
      case SubstructureSearch::SS_DOUBLE_BOND:
        bt = DOUBLE_BOND;
        break;
      default:
        cerr << "Toggle_Kekule_Form::ConstructFromProto:unrecognized bond type " << proto.ShortDebugString() << '\n';
        return 0;
    }

    add_bond(new Bond(bond.a1(), bond.a2(), bt));
  }

  return _bond.number_elements();
}

