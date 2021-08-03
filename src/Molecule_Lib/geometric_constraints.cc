#include "geometric_constraints.h"

namespace geometric_constraints {

using std::cerr;

// Initialize with arbitrary values that will prevent matching anything.
AllowedRange::AllowedRange() {
  _min_value = 1.0;
  _max_value = 0.0;
}

std::ostream&
operator<< (std::ostream& output, const AllowedRange& arange) {
  output << "AllowedRange [" << arange._min_value << ',' << arange._max_value << ']';
  return output;
}

int
AllowedRange::BuildFromProto(const GeometricConstraints::Range& proto) {
  _min_value = proto.min();
  _max_value = proto.max();

  // If an upper bound is not specified, assume infinite.
  if (_min_value > 0.0 && _max_value == 0.0) {
    _max_value = std::numeric_limits<float>::max();
  }

  if (_min_value > _max_value) {
    cerr << "AllowedRange::BuildFromProto:invalid range " << *this << '\n';
    return 0;
  }

  if (proto.one_time_scaling_factor() != 0.0) {
    _min_value *= proto.one_time_scaling_factor();
    _max_value *= proto.one_time_scaling_factor();
  }

  return 1;
}

int
ConstraintBaseClass::IsValid() const {
  // Unset is a valid state.
  if (!_allowed_range.Active() && _atoms.empty()) {
    return 1;
  }

  const int natoms = _atoms.number_elements();

  if (natoms < 2)  // The DistanceConstraint has 2 atoms.
    return 0;
  if (natoms > 4)  // Torsion has 4 atoms.
    return 0;

  // We have atoms specified, but no _allowed_range.
  if (!_allowed_range.Active()) {
    return 0;
  }

  // Are all the atom numbers unique?
  for (int i = 0; i < natoms; ++i) {
    int a1 = _atoms[i];
    for (int j = i + 1; j < natoms; ++j) {
      if (_atoms[j] == a1)
        return 0;
    }
  }

  return 1;
}

int
ConstraintBaseClass::SetAtoms(int a1, int a2) {
  return (_atoms.add_if_not_already_present(a1)> 0) +
         (_atoms.add_if_not_already_present(a2)> 0);
}

int
ConstraintBaseClass::SetAtoms(int a1, int a2, int a3) {
  return (_atoms.add_if_not_already_present(a1)> 0) +
         (_atoms.add_if_not_already_present(a2)> 0) +
         (_atoms.add_if_not_already_present(a3)> 0);
}

int
ConstraintBaseClass::SetAtoms(int a1, int a2, int a3, int a4) {
  return (_atoms.add_if_not_already_present(a1)> 0) +
         (_atoms.add_if_not_already_present(a2)> 0) +
         (_atoms.add_if_not_already_present(a3)> 0) +
         (_atoms.add_if_not_already_present(a4)> 0);
}

int
ConstraintBaseClass::AtomNumbersPresent(resizable_array<int>& atom_numbers) const {
  int rc = 0;
  for (int a : _atoms) {
    rc += atom_numbers.add_if_not_already_present(a);
  }
  return rc;
}

std::ostream& operator<<(std::ostream& output, const ConstraintBaseClass& constraint) {
  for (int i : constraint._atoms) {
    output << ' ' << i;
  }
  output << ' ' << constraint._allowed_range;

  return output;
}

int
DistanceConstraint::BuildFromProto(const GeometricConstraints::Distance& proto) {
  if (! _allowed_range.BuildFromProto(proto.range())) {
    cerr << "DistanceConstraint::BuildFromProto:invalid range " << proto.ShortDebugString() << '\n';
    return 0;
  }
  if (SetAtoms(proto.a1(), proto.a2()) != 2) {
    cerr << "DistanceConstraint::BuildFromProto:invalid atom numbers " << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (! IsValid()) {
    cerr << "DistanceConstraint::BuildFromProto:invalid " << *this << '\n';
    return 0;
  }

  return 1;
}


std::ostream& operator<< (std::ostream& output, const DistanceConstraint& constraint) {
  const ConstraintBaseClass& me = constraint;
  output << "DistanceConstraint:atoms " << me;
  return output;
}

int
DistanceConstraint::Matches(const Molecule& m) {
  return _allowed_range.Matches(m.distance_between_atoms(_atoms[0], _atoms[1]));
}

int
DistanceConstraint::Matches(const Molecule& m, const Set_of_Atoms& embedding) {
  return _allowed_range.Matches(m.distance_between_atoms(embedding[_atoms[0]], embedding[_atoms[1]]));
}

int
BondAngleConstraint::BuildFromProto(const GeometricConstraints::BondAngle& proto) {
  if (! _allowed_range.BuildFromProto(proto.range())) {
    cerr << "DistanceConstraint::BuildFromProto:invalid range " << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (SetAtoms(proto.a1(), proto.a2(), proto.a3()) != 3) {
    cerr << "BondAngle::BuildFromProto:invalid atom numbers " << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (! IsValid()) {
    cerr << "BondAngleConstraint::BuildFromProto:invalid " << *this << '\n';
    return 0;
  }

  return 1;
}

std::ostream&
operator<< (std::ostream& output, const BondAngleConstraint& constraint) {
  const ConstraintBaseClass& me = constraint;
  output << "BondAngleConstraint:atoms " <<  me;
  return output;
}

int
BondAngleConstraint::Matches(const Molecule& m) {
//std::cerr << "Angle is " << m.bond_angle(_a1, _a2, _a3) << " radians " << _allowed_range << '\n';
  return _allowed_range.Matches(m.bond_angle(_atoms[0], _atoms[1], _atoms[2]));
}

int
BondAngleConstraint::Matches(const Molecule& m, const Set_of_Atoms& embedding) {
  return _allowed_range.Matches(m.bond_angle(embedding[_atoms[0]], embedding[_atoms[1]], embedding[_atoms[2]]));
}

int
TorsionAngleConstraint::BuildFromProto(const GeometricConstraints::TorsionAngle& proto) {
  if (! _allowed_range.BuildFromProto(proto.range())) {
    cerr << "DistanceConstraint::BuildFromProto:invalid range " << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (SetAtoms(proto.a1(), proto.a2(), proto.a3(), proto.a4()) != 4) {
    cerr << "TorsionAngle::BuildFromProto:invalid atom numbers " << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (! IsValid()) {
    cerr << "TorsionAngle::BuildFromProto:invalid " << *this << '\n';
    return 0;
  }

  return 1;
}

std::ostream&
operator<< (std::ostream& output, const TorsionAngleConstraint& constraint) {
  const ConstraintBaseClass& me = constraint;
  output << "TorsionAngleConstraint:atoms " << me;
  return output;
}
int
TorsionAngleConstraint::Matches(const Molecule& m) {
  return _allowed_range.Matches(m.dihedral_angle(_atoms[0], _atoms[1], _atoms[2], _atoms[3]));
}

int
TorsionAngleConstraint::Matches(const Molecule& m, const Set_of_Atoms& embedding) {
  return _allowed_range.Matches(m.dihedral_angle(embedding[_atoms[0]], embedding[_atoms[1]], embedding[_atoms[2]], embedding[_atoms[3]]));
}

SetOfGeometricConstraints::SetOfGeometricConstraints() {
  _active = 0;
  _number_to_match = 0;
}

// The number of constraints across all different kinds.
int
SetOfGeometricConstraints::_number_constraints() const {
  return _distances.number_elements() +
         _bond_angles.number_elements() +
         _torsion_angles.number_elements();
}

int
SetOfGeometricConstraints::IsValid() const {
  int nset = _number_constraints();
  if (_number_to_match > nset) {
    return 0;
  }

  if (nset == 0 &&_number_to_match > 0) {
    return 0;
  }

  return 1;
}

int
SetOfGeometricConstraints::BuildFromProto(const GeometricConstraints::SetOfConstraints& proto) {
  if (proto.number_to_match() > 0) {
    _number_to_match = proto.number_to_match();
  }
  for (const auto& dist : proto.distances()) {
    std::unique_ptr<DistanceConstraint> c = std::make_unique<DistanceConstraint>();
    if (! c->BuildFromProto(dist)) {
      cerr << "SetOfConstraints::BuildFromProto:bad distance " << dist.ShortDebugString() << '\n';
      return 0;
    }
    _distances << c.release();
  }

  for (const auto& angle : proto.bond_angles()) {
    std::unique_ptr<BondAngleConstraint> c = std::make_unique<BondAngleConstraint>();
    if (! c->BuildFromProto(angle)) {
      cerr << "SetOfConstraints::BuildFromProto:bad angle " << angle.ShortDebugString() << '\n';
      return 0;
    }
    _bond_angles << c.release();
  }

  for (const auto& torsion : proto.torsion_angles()) {
    std::unique_ptr<TorsionAngleConstraint> c = std::make_unique<TorsionAngleConstraint>();
    if (! c->BuildFromProto(torsion)) {
      cerr << "SetOfConstraints::BuildFromProto:bad torsion " << torsion.ShortDebugString() << '\n';
      return 0;
    }
    _torsion_angles << c.release();
  }

  _active = 1;
  if (_number_to_match == 0) {
    _number_to_match =  _number_constraints();
  }

  return IsValid();
}

int
SetOfGeometricConstraints::Matches(const Molecule& m) const {
  int matched_here = 0;
  for (auto c : _distances) {
    if (c->Matches(m)) {
      matched_here++;
      if (matched_here >= _number_to_match) {
        return matched_here;
      }
    }
  }
  for (auto c : _bond_angles) {
    if (c->Matches(m)) {
      matched_here++;
      if (matched_here >= _number_to_match) {
        return matched_here;
      }
    }
  }
  for (auto c : _torsion_angles) {
    if (c->Matches(m)) {
      matched_here++;
      if (matched_here >= _number_to_match) {
        return matched_here;
      }
    }
  }

  return 0;
}

int
SetOfGeometricConstraints::Matches(const Molecule& m, const Set_of_Atoms& embedding) const {
  int matched_here = 0;
  for (auto c : _distances) {
    if (c->Matches(m, embedding)) {
      matched_here++;
      if (matched_here >= _number_to_match) {
        return matched_here;
      }
    }
  }
  for (auto c : _bond_angles) {
    if (c->Matches(m, embedding)) {
      matched_here++;
      if (matched_here >= _number_to_match) {
        return matched_here;
      }
    }
  }
  for (auto c : _torsion_angles) {
    if (c->Matches(m, embedding)) {
      matched_here++;
      if (matched_here >= _number_to_match) {
        return matched_here;
      }
    }
  }

  return 0;
}

template <typename T>
void
write_constraints(const resizable_array_p<T>& constraints, std::ostream& output) {
  for (const T * constraint : constraints) {
    output << ' ' << *constraint;
  }
  output << '\n';
}
void
SetOfGeometricConstraints::DebugPrint(std::ostream& output) const {
  output << "SetOfGeometricConstraints\n";
  if (_distances.size()) {
    output << " distances";
    write_constraints(_distances, output);
  }
  if (_bond_angles.size()) {
    output << " bond angles";
    write_constraints(_bond_angles, output);
  }
  if (_torsion_angles.size()) {
    output << " torsion angles";
    write_constraints(_torsion_angles, output);
  }
  output << "must match " << _number_to_match << '\n';
}

resizable_array<int>
SetOfGeometricConstraints::AtomNumbersPresent() const {
  resizable_array<int> result;
  for (const auto * c : _distances) {
    c->AtomNumbersPresent(result);
  }
  for (const auto * c : _bond_angles) {
    c->AtomNumbersPresent(result);
  }
  for (const auto * c : _torsion_angles) {
    c->AtomNumbersPresent(result);
  }

  return result;
}

}  // namespace geometric_constraints
