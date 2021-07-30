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
  output << "AllowedRange:btw [" << arange._min_value << ',' << arange._max_value << ']';
  return output;
}

int
AllowedRange::BuildFromProto(const GeometricConstraints::Range& proto) {
  _min_value = proto.min();
  _max_value = proto.max();
  if (_min_value > _max_value) {
    cerr << "AllowedRange::BuildFromProto:invalid range " << *this << '\n';
    return 0;
  }
  return 1;
}

DistanceConstraint::DistanceConstraint() {
  _a1 = -1;
  _a2 = -1;
}

std::ostream& operator<< (std::ostream& output, const DistanceConstraint& constraint) {
  output << "DistanceConstraint:atoms " << constraint._a1 << ',' << constraint._a2 << ' ';
  output << constraint._allowed_range;
  return output;
}

int
DistanceConstraint::BuildFromProto(const GeometricConstraints::Distance& proto) {
  if (! _allowed_range.BuildFromProto(proto.range())) {
    cerr << "DistanceConstraint::BuildFromProto:invalid range " << proto.ShortDebugString() << '\n';
    return 0;
  }

  _a1 = proto.a1();
  _a2 = proto.a2();
  if (! IsValid()) {
    cerr << "DistanceConstraint::BuildFromProto:invalid " << *this << '\n';
    return 0;
  }

  return 1;
}

int
DistanceConstraint::IsValid() const {
  if (! _allowed_range.Active()) {
    return 0;
  }

  if (_a1 < 0 || _a2 < 0) {
    return 0;
  }

  return _a1 != _a2;
}

int
DistanceConstraint::Matches(const Molecule& m) {
  return _allowed_range.Matches(m.distance_between_atoms(_a1, _a2));
}

int
DistanceConstraint::Matches(const Molecule& m, const Set_of_Atoms& embedding) {
  return _allowed_range.Matches(m.distance_between_atoms(embedding[_a1], embedding[_a2]));
}

BondAngleConstraint::BondAngleConstraint() {
  _a1 = -1;
  _a2 = -1;
  _a3 = -1;
}

int
BondAngleConstraint::BuildFromProto(const GeometricConstraints::BondAngle& proto) {
  if (! _allowed_range.BuildFromProto(proto.range())) {
    cerr << "DistanceConstraint::BuildFromProto:invalid range " << proto.ShortDebugString() << '\n';
    return 0;
  }

  _a1 = proto.a1();
  _a2 = proto.a2();
  _a3 = proto.a3();
  if (! IsValid()) {
    cerr << "BondAngleConstraint::BuildFromProto:invalid " << *this << '\n';
    return 0;
  }

  return 1;
}

std::ostream&
operator<< (std::ostream& output, const BondAngleConstraint& constraint) {
  output << "BondAngleConstraint:atoms " << constraint._a1 << ' ' <<
            constraint._a2 << ' ' <<
            constraint._a3 << ' ';
  output << constraint._allowed_range;
  return output;
}

int
BondAngleConstraint::IsValid() const {
  if (! _allowed_range.Active()) {
    return 0;
  }

  if (_a1 < 0 || _a2 < 0 || _a3 < 0) {
    return 0;
  }

  if (_a1 == _a2 || _a1 == _a3) {
    return 0;
  }

  if (_a2 == _a3) {
    return 0;
  }

  return 1;
}

int
BondAngleConstraint::Matches(const Molecule& m) {
//std::cerr << "Angle is " << m.bond_angle(_a1, _a2, _a3) << " radians " << _allowed_range << '\n';
  return _allowed_range.Matches(m.bond_angle(_a1, _a2, _a3));
}

int
BondAngleConstraint::Matches(const Molecule& m, const Set_of_Atoms& embedding) {
  return _allowed_range.Matches(m.bond_angle(embedding[_a1], embedding[_a2], embedding[_a3]));
}

int
TorsionAngleConstraint::BuildFromProto(const GeometricConstraints::TorsionAngle& proto) {
  if (! _allowed_range.BuildFromProto(proto.range())) {
    cerr << "DistanceConstraint::BuildFromProto:invalid range " << proto.ShortDebugString() << '\n';
    return 0;
  }

  _a1 = proto.a1();
  _a2 = proto.a2();
  _a3 = proto.a3();
  _a4 = proto.a4();
  if (! IsValid()) {
    cerr << "TorsionAngle::BuildFromProto:invalid " << *this << '\n';
    return 0;
  }

  return 1;
}

std::ostream&
operator<< (std::ostream& output, const TorsionAngleConstraint& constraint) {
  output << "TorsionAngleConstraint:atoms " << constraint._a1 << ' ' <<
            constraint._a2 << ' ' <<
            constraint._a3 << ' ' <<
            constraint._a4 << ' ';
  output << constraint._allowed_range;
  return output;
}

TorsionAngleConstraint::TorsionAngleConstraint() {
  _a1 = -1;
  _a2 = -1;
  _a3 = -1;
  _a4 = -1;
}

int
TorsionAngleConstraint::IsValid() const {
  if (! _allowed_range.Active()) {
    return 0;
  }

  if (_a1 < 0 || _a2 < 0 || _a3 < 0 || _a4 < 0) {
    return 0;
  }

  if (_a1 == _a2 || _a1 == _a3 || _a1 == _a4) {
    return 0;
  }

  if (_a2 == _a3 || _a2 == _a4) {
    return 0;
  }

  if (_a3 == _a4) {
    return 0;
  }

  return 1;
}

int
TorsionAngleConstraint::Matches(const Molecule& m) {
  return _allowed_range.Matches(m.dihedral_angle(_a1, _a2, _a3, _a4));
}

int
TorsionAngleConstraint::Matches(const Molecule& m, const Set_of_Atoms& embedding) {
  return _allowed_range.Matches(m.dihedral_angle(embedding[_a1], embedding[_a2], embedding[_a3], embedding[_a4]));
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
      if (matched_here > _number_to_match) {
        return matched_here;
      }
    }
  }
  for (auto c : _bond_angles) {
    if (c->Matches(m)) {
      matched_here++;
      if (matched_here > _number_to_match) {
        return matched_here;
      }
    }
  }
  for (auto c : _torsion_angles) {
    if (c->Matches(m)) {
      matched_here++;
      if (matched_here > _number_to_match) {
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
      if (matched_here > _number_to_match) {
        return matched_here;
      }
    }
  }
  for (auto c : _bond_angles) {
    if (c->Matches(m, embedding)) {
      matched_here++;
      if (matched_here > _number_to_match) {
        return matched_here;
      }
    }
  }
  for (auto c : _torsion_angles) {
    if (c->Matches(m, embedding)) {
      matched_here++;
      if (matched_here > _number_to_match) {
        return matched_here;
      }
    }
  }

  return 0;
}

}  // namespace geometric_constraints
