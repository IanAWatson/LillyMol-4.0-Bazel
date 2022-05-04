#include "re2/re2.h"

#include "Molecule_Tools/ring_closure_enumeration_lib.h"

namespace ring_closure_enumeration {

BondFormation::BondFormation() {
  _isotope = 0;
  _bond_type = 0;
  _ring_closure = 0;
  _as_char[1] = '-';
  _as_char[2] = '=';
  _as_char[3] = '#';
}

int
BondFormation::Build(const const_IWSubstring& buffer) {
  std::unique_ptr<RE2> rx = std::make_unique<RE2>("^(\\d+):(\\d+)(.)(\\d+)$");
  re2::StringPiece tmp(buffer.data(), buffer.length());
  char btype;
  if (! RE2::FullMatch(tmp, *rx, &_component, &_isotope, &btype,  &_ring_closure)) {
    cerr << "BondFormation::Build:cannot match rx '" << buffer << "'\n";
    return 0;
  }

  switch (btype) {
    case '-':
      _bond_type = 1;
      break;
    case '=':
      _bond_type = 2;
      break;
    case '#':
      _bond_type = 3;
      break;
    default:
      cerr << "BondFormation::Build:Unrecognised bond type '" << btype << "'\n";
      return 0;
  }

  if (_isotope == 0) {
    cerr << "BondFormation::Build:zero isotope\n";
    return 0;
  }

  if (_ring_closure == 0) {
    cerr << "BondFormation::Build:zero ring closure\n";
    return 0;
  }

  return 1;

  const_IWSubstring num1, num2;
  if (buffer.split(num1, '-', num2) && ! num1.empty() && ! num2.empty()) {
    _bond_type = 1;
  } else if (buffer.split(num1, '=', num2) && ! num1.empty() && ! num2.empty()) {
    _bond_type = 2;
  } else if (buffer.split(num1, '#', num2) && ! num1.empty() && ! num2.empty()) {
    _bond_type = 3;
  } else {
    cerr << "BondFormation::Build:cannot parse '" << buffer << "'\n";
    return 0;
  }

  if (! num1.numeric_value(_isotope) || _isotope < 1) {
    cerr << "BondFormation::buffer:invalid isotope '" << buffer << "'\n";
    return 0;
  }

  if (! num2.numeric_value(_ring_closure) || _ring_closure < 1) {
    cerr << "BondFormation::buffer:invalid ring closure '" << buffer << "'\n";
    return 0;
  }

  return 1;
}

int
Reagents::Read(const char * fname) {
  data_source_and_type<Molecule> input(FILE_TYPE_SMI, fname);
  if (! input.good()) {
    cerr << "Reagents::Read:cannot open '" << fname << "'\n";
    return 0;
  }

  return Read(input);
}

// no-op until I figure out how this should be handled.
int
Preprocess(Molecule& m) {
  return 1;
}

int
Reagents::Read(data_source_and_type<Molecule>& input) {
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    if (! Preprocess(*m)) {
      delete m;
      continue;
    }
    _mols << m;
  }

  return _mols.number_elements();
}

int
Reagents::CreateFragments(const BondFormation& bond_formation) {
  _smiles.resize(_mols.number_elements());

  IWString change_from;
  change_from << "\\[" << bond_formation.isotope() << "([[:alnum:]]+)\\]";
  IWString change_to;
  change_to << "[" << bond_formation.isotope() << "\\1]" << bond_formation.btype() << '%' << bond_formation.ring_closure();
  std::string change_to_string = change_to.AsString();

  std::unique_ptr<RE2> rx = std::make_unique<RE2>(change_from.AsString());
  for (Molecule* m : _mols) {
    std::string smiles = m->smiles().AsString();
    if (! RE2::Replace(&smiles, *rx.get(), change_to_string)) {
      cerr << "Cannot replace " << change_from << " with " << change_to << " in '" << smiles << "'\n";
      return 0;
    }
    IWString *smi = new IWString(smiles);
    _smiles << smi;
  }

  return _smiles.number_elements();
}

int
Reagents::CreateFragments(const resizable_array<BondFormation*>& bond_formations) {
  _smiles.resize(_mols.number_elements());

  resizable_array_p<RE2> transformations(bond_formations.number_elements());
  std::vector<std::string> change_to;
  change_to.reserve(bond_formations.number_elements());
  for (const BondFormation* bond_formation : bond_formations) {
    IWString change_from;
    change_from << "\\[" << bond_formation->isotope() << "([[:alnum:]]+)\\]";
    RE2* rx = new RE2(change_from.AsString());
    transformations << rx;

    IWString local_change_to;
    local_change_to << "[" << bond_formation->isotope() << "\\1]" << bond_formation->btype() << '%' << bond_formation->ring_closure();
    change_to.emplace_back(std::move(local_change_to.AsString()));
  }

  for (Molecule* m : _mols) {
    std::string smiles = m->smiles().AsString();
    for (int i = 0; i < transformations.number_elements(); ++i) {
      if (! RE2::Replace(&smiles, *transformations[i], change_to[i])) {
        cerr << "Cannot replace " << transformations[i]->pattern() << " with " << change_to[i] << " in '" << smiles << "'\n";
        return 0;
      }
    }
    IWString *smi = new IWString(smiles);
    _smiles << smi;
  }

  return _smiles.number_elements();
}

}  // namespace ring_closure_enumeration
