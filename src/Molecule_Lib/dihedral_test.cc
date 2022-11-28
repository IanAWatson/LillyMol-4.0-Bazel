// Tests for the dihedral angle functions

#include <array>
#include <cmath>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "molecule.h"
#include "misc2.h"

namespace {

using testing::FloatEq;

struct CoordsAngle {
  // Coordinates for 4 atoms.
  std::array<float, 3> coords[4];
  // The expected dihedral angle
  float expected_angle;
  // A rotation axis. The dihedral should be invariant to
  // rotations about this axis.
  std::array<float, 3> rot;
};

class TestDihedral : public testing::TestWithParam<CoordsAngle> {
  protected:
    Molecule mol;
};

class TestSignedDihedral : public testing::TestWithParam<CoordsAngle> {
  protected:
    Molecule mol;
};

TEST_P(TestDihedral, TestAngles) {
  mol.build_from_smiles("CCCC");
  const auto params = GetParam();

  for (int i = 0; i < 4; ++i) {
    mol.setxyz(i, params.coords[i][0], params.coords[i][1], params.coords[i][2]);
  }

  const float angle = mol.dihedral_angle(0, 1, 2, 3) * RAD2DEG;
  //std::cerr << " Initial config angle " << angle << '\n';
  EXPECT_NEAR(angle, params.expected_angle, 0.0001);

  Coordinates rot(params.rot[0], params.rot[1], params.rot[2]);
  rot.normalise();

  // We will do a series of rotations by 30 degrees.
  constexpr float k30 = 30.0 * RAD2DEG;
  for (int i = 1; i < 12; ++i) {
    // Reset coordinates each time to avoid float inaccuracies.
    for (int j = 0; j < 4; ++j) {
      mol.setxyz(j, params.coords[j][0], params.coords[j][1], params.coords[j][2]);
    }
    mol.rotate_atoms(rot, i * k30);
    const float angle = mol.dihedral_angle(0, 1, 2, 3) * RAD2DEG;
    //std::cerr << " i " << i << " angle " << angle << '\n';
    EXPECT_NEAR(angle, params.expected_angle, 0.0001);
  }
}
INSTANTIATE_TEST_SUITE_P(TestDihedral, TestDihedral, testing::Values(
  CoordsAngle{{{-1.0, 1.0, 0.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}}, 0.0, {1,2,3}},
  CoordsAngle{{{-1.0, -1.0, 0.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}}, 180.0, {-1,2,3}},
  CoordsAngle{{{-1.0, -1.0, 0.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 0.0, 1.0}}, 90.0, {0,2,-1}},
  CoordsAngle{{{-1.0, -1.0, 0.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 0.0, -1.0}}, 90.0, {-3,-2,-1}}
));

TEST_P(TestSignedDihedral, TestAngles) {
  mol.build_from_smiles("CCCC");
  const auto params = GetParam();

  for (int i = 0; i < 4; ++i) {
    mol.setxyz(i, params.coords[i][0], params.coords[i][1], params.coords[i][2]);
  }

  const float angle = mol.signed_dihedral_angle(0, 1, 2, 3) * RAD2DEG;
  // std::cerr << " Initial config angle " << angle << " expected " << params.expected_angle << '\n';
  EXPECT_NEAR(angle, params.expected_angle, 0.0001);

  Coordinates rot(params.rot[0], params.rot[1], params.rot[2]);
  rot.normalise();

  // We will do a series of rotations by 30 degrees.
  constexpr float k30 = 30.0 * RAD2DEG;
  for (int i = 1; i < 12; ++i) {
    // Reset coordinates each time to avoid float inaccuracies.
    for (int j = 0; j < 4; ++j) {
      mol.setxyz(j, params.coords[j][0], params.coords[j][1], params.coords[j][2]);
    }
    mol.rotate_atoms(rot, i * k30);
    const float angle = mol.signed_dihedral_angle(0, 1, 2, 3) * RAD2DEG;
    // std::cerr << " i " << i << " angle " << angle << '\n';
    EXPECT_NEAR(angle, params.expected_angle, 0.0001);
  }
}
INSTANTIATE_TEST_SUITE_P(TestSignedDihedral, TestSignedDihedral, testing::Values(
  CoordsAngle{{{-1.0, 1.0, 0.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}}, 0.0, {1,2,3}},
  CoordsAngle{{{-1.0, 1.0, -0.01}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}}, -0.57293868064880371, {1,2,3}},
  CoordsAngle{{{-1.0, 1.0, 0.01}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}}, 0.57293689250946045, {1,2,3}},
  CoordsAngle{{{-1.0, -1.0, 0.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}}, 180.0, {-1,2,3}},
  CoordsAngle{{{-1.0, -1.0, 0.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 0.0, 1.0}}, 90.0, {0,2,-1}},
  CoordsAngle{{{-1.0, -1.0, 0.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 0.0, -1.0}}, -90.0, {-3,-2,-1}}
));

struct MolDihedralAngle {
  IWString smiles;
  std::array<int, 4> atoms;
  float angle;
  std::array<float, 3> rot;
};

class TestSetDihedral : public testing::TestWithParam<MolDihedralAngle> {
  protected:
    Molecule _m;
};

TEST_P(TestSetDihedral, TestAngles) {
  const auto params = GetParam();
  _m.build_from_smiles(params.smiles);
#ifdef ECHO_MOL
  _m.debug_print(std::cerr);
  for (int i = 0; i < 4; ++i) {
    std::cerr << _m.x(i) << ',' << _m.y(i) << ',' << _m.z(i) << '\n';
  }
#endif

  std::unique_ptr<float[]> initial_coords = _m.GetCoords();

  const int a0 = params.atoms[0];
  const int a1 = params.atoms[1];
  const int a2 = params.atoms[2];
  const int a3 = params.atoms[3];

  Coordinates rot(params.rot[0], params.rot[1], params.rot[2]);
  rot.normalise();

  constexpr float k30 = 30.0 * DEG2RAD;

  _m.set_dihedral(a0, a1, a2, a3, params.angle);
  EXPECT_NEAR(_m.signed_dihedral_angle(a0, a1, a2, a3), params.angle, 0.001);
  // Setting it again should result in same outcome.
  _m.set_dihedral(a0, a1, a2, a3, params.angle);
  EXPECT_NEAR(_m.signed_dihedral_angle(a0, a1, a2, a3), params.angle, 0.001);

  for (int i = 0; i < 12; ++i) {
    _m.SetXyz(initial_coords.get());
    _m.rotate_atoms(rot, i * k30);
    _m.set_dihedral(a0, a1, a2, a3, params.angle);
    EXPECT_NEAR(_m.signed_dihedral_angle(a0, a1, a2, a3), params.angle, 0.001);
  }
}

INSTANTIATE_TEST_SUITE_P(TestSetDihedral, TestSetDihedral, testing::Values(
  MolDihedralAngle{ "C{{0,1,0}}C{{0,0,0}}C{{1,0,0}}C{{1,1,0}}", {0, 1, 2, 3}, 0.0, {0.5, -0.2, 0.2}},
  MolDihedralAngle{ "C{{0,1,0}}C{{0,0,0}}C{{1,0,0}}C{{1,1,0}}", {0, 1, 2, 3}, 0.1, {-0.2, 0.4, -0.3}},
  MolDihedralAngle{ "C{{0,1,0}}C{{0,0,0}}C{{1,0,0}}C{{1,1,0}}", {0, 1, 2, 3}, 0.2, {0.9, 0.4, -0.3}},
  MolDihedralAngle{ "C{{0,1,0}}C{{0,0,0}}C{{1,0,0}}C{{1,1,0}}", {0, 1, 2, 3}, -0.2, {-0.9, 0.2, -0.1}},
  MolDihedralAngle{ "C{{0,1,0}}C{{0,0,0}}C{{1,0,0}}C{{1,1,0.2}}", {0, 1, 2, 3}, 0.4, {-0.9, 0.4, -0.2}}
));

TEST(TestRotate, TestX) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C{{0,0,0}}C{{1,0,0}}C{{1,1,0}}"));

  const Space_Vector<float> xaxis(1.0f, 0.0f, 0.0f);
  const Space_Vector<float> yaxis(0.0f, 1.0f, 0.0f);

  constexpr float kAbs = 0.0001;

  m.rotate_atoms(xaxis, 0.0f);
  EXPECT_NEAR(m.x(0), 0.0f, kAbs);
  EXPECT_NEAR(m.y(0), 0.0f, kAbs);
  EXPECT_NEAR(m.z(0), 0.0f, kAbs);
  EXPECT_NEAR(m.x(1), 1.0f, kAbs);
  EXPECT_NEAR(m.y(1), 0.0f, kAbs);
  EXPECT_NEAR(m.z(1), 0.0f, kAbs);
  EXPECT_NEAR(m.x(2), 1.0f, kAbs);
  EXPECT_NEAR(m.y(2), 1.0f, kAbs);
  EXPECT_NEAR(m.z(2), 0.0f, kAbs);

  std::unique_ptr<float[]> initial_coords = m.GetCoords();
  m.rotate_atoms(xaxis, 90.0 * DEG2RAD);
  EXPECT_NEAR(m.x(2), 1.0f, kAbs);
  EXPECT_NEAR(m.y(2), 0.0f, kAbs);
  EXPECT_NEAR(m.z(2), 1.0f, kAbs);
  Coordinates f = m.get_coords(2);
  f.normalise();
  EXPECT_NEAR(f.angle_between_unit_vectors(yaxis), iwmisc::Deg2Rad(90.0), kAbs);

  m.SetXyz(initial_coords.get());
  m.rotate_atoms(xaxis, 45.0f * DEG2RAD);
  EXPECT_NEAR(m.x(2), 1.0f, kAbs);
  EXPECT_NEAR(m.y(2), sqrt(2.0f) * 0.5f, kAbs);
  EXPECT_NEAR(m.z(2), sqrt(2.0f) * 0.5f, kAbs);
  EXPECT_NEAR(m.z(2) / m.y(2), std::tan(iwmisc::Deg2Rad(45.0f)), kAbs);
  f = m.get_coords(2);
  f.set_x(0.0f);
  f.normalise();
  EXPECT_NEAR(f.angle_between_unit_vectors(yaxis), iwmisc::Deg2Rad(45.0), kAbs);

  m.SetXyz(initial_coords.get());
  m.rotate_atoms(xaxis, -45.0f * DEG2RAD);
  EXPECT_NEAR(m.x(2), 1.0f, kAbs);
  EXPECT_NEAR(m.y(2), sqrt(2.0f) * 0.5f, kAbs);
  EXPECT_NEAR(m.z(2), -sqrt(2.0f) * 0.5f, kAbs);
  EXPECT_NEAR(m.z(2) / m.y(2), std::tan(iwmisc::Deg2Rad(-45.0f)), kAbs);
  f = m.get_coords(2);
  f.set_x(0.0f);
  f.normalise();
  EXPECT_NEAR(f.angle_between_unit_vectors(yaxis), iwmisc::Deg2Rad(45.0), kAbs);

  m.SetXyz(initial_coords.get());
  m.rotate_atoms(xaxis, iwmisc::Deg2Rad(1.0f));
  EXPECT_NEAR(m.x(2), 1.0f, kAbs);
  EXPECT_NEAR(m.y(2), 0.999847, kAbs);
  EXPECT_NEAR(m.z(2), 0.0174524, kAbs);
  EXPECT_NEAR(m.z(2) / m.y(2), std::tan(iwmisc::Deg2Rad(1.0f)), kAbs);
  f = m.get_coords(2);
  f.set_x(0.0f);
  f.normalise();
  EXPECT_NEAR(f.angle_between_unit_vectors(yaxis), iwmisc::Deg2Rad(1.0), kAbs);

  m.SetXyz(initial_coords.get());
  m.rotate_atoms(xaxis, iwmisc::Deg2Rad(5.0f));
  EXPECT_NEAR(m.x(2), 1.0f, kAbs);
  EXPECT_NEAR(m.y(2), 0.996194, kAbs);
  EXPECT_NEAR(m.z(2), 0.08715573, kAbs);
  EXPECT_NEAR(m.z(2) / m.y(2), std::tan(iwmisc::Deg2Rad(5.0f)), kAbs);
  f = m.get_coords(2);
  f.set_x(0.0f);
  f.normalise();
  EXPECT_NEAR(f.angle_between_unit_vectors(yaxis), iwmisc::Deg2Rad(5.0), kAbs);

  m.SetXyz(initial_coords.get());
  m.rotate_atoms(xaxis, iwmisc::Deg2Rad(30.0f));
  EXPECT_NEAR(m.x(2), 1.0f, kAbs);
  EXPECT_NEAR(m.y(2), 0.8660253, kAbs);
  EXPECT_NEAR(m.z(2), 0.5, kAbs);
  EXPECT_NEAR(m.z(2) / m.y(2), std::tan(iwmisc::Deg2Rad(30.0f)), kAbs);
  f = m.get_coords(2);
  f.set_x(0.0f);
  f.normalise();
  EXPECT_NEAR(f.angle_between_unit_vectors(yaxis), iwmisc::Deg2Rad(30.0), kAbs);

  m.SetXyz(initial_coords.get());
  m.rotate_atoms(xaxis, iwmisc::Deg2Rad(60.0f));
  EXPECT_NEAR(m.x(2), 1.0f, kAbs);
  EXPECT_NEAR(m.y(2), 0.5, kAbs);
  EXPECT_NEAR(m.z(2), 0.8660254, kAbs);
  EXPECT_NEAR(m.z(2) / m.y(2), std::tan(iwmisc::Deg2Rad(60.0f)), kAbs);
  f = m.get_coords(2);
  f.set_x(0.0f);
  f.normalise();
  EXPECT_NEAR(f.angle_between_unit_vectors(yaxis), iwmisc::Deg2Rad(60.0), kAbs);
}

}  // namespace
