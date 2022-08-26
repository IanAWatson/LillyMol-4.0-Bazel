// Given a spatial arrangement of atoms, discern plausible topologies

#include <algorithm>
#include <iostream>
#include <memory>
#include <unordered_map>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/combinations.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/bond_length_distribution.pb.h"

namespace topology_from_geom {

using std::cerr;

void
Usage(int rc) {
  cerr << "Discerns plausible bond topologies from atomic positions\n";

  cerr << " -v                verbose output\n";

  ::exit(rc);
}

// Data about a particular bond length distribution. Instantiated from a
// BondLengthDistribution.Distribution proto.
class EmpiricalBondLengthDistribution {
  private:
    int _atomic_number1;
    int _atomic_number2;
    int _btype;
    float _min_distance;
    float _dx;
    float _max_distance;
    float* _distance;
    uint32_t* _count;
    float* _score;

  public:
    EmpiricalBondLengthDistribution();
    ~EmpiricalBondLengthDistribution();

    int Build(const BondLengthDistribution::Distribution& bld);

    float Score(float distance) const;

    float max_distance() const {
      return _max_distance;
    }
};

EmpiricalBondLengthDistribution::EmpiricalBondLengthDistribution() {
  _atomic_number1 = 0;
  _atomic_number2 = 0;
  _btype = 0;
  _min_distance = 0.0;
  _dx = 0.0;
  _max_distance = 0.0;
  _distance = nullptr;
  _count = nullptr;
  _score= nullptr;
}

EmpiricalBondLengthDistribution::~EmpiricalBondLengthDistribution() {
  if (_distance != nullptr) {
    delete [] _distance;
  }
  if (_count != nullptr) {
    delete [] _count;
  }
  if (_score != nullptr) {
    delete [] _score;
  }
}

float
EmpiricalBondLengthDistribution::Score(float distance) const {
  if (distance <= _min_distance || distance >= _max_distance) {
    return 0.0;
  }

  const int bucket = static_cast<int>((distance - _min_distance) / _dx + 0.4999f);
  return _score[bucket];
}

// Contains multiple EmpiricalBondLengthDistribution. The important
// function is to be able to return the pdf for two atomic numbers.
// For each pair of atomic numbers, we store an EmpiricalBondLengthDistribution.
class EmpiricalBondLengthDistributions {
  private:
    std::unordered_map<int, EmpiricalBondLengthDistribution> _dist;

  // private functions
    int Index(int at1, int at2, int btype) const;

  public:
    int Build(const BondLengthDistribution::Distributions& proto);

    float LongestPlausibleBond() const;

    // The number of different bond types we have.
    uint32_t size() const {
      return _dist.size();
    }

    // Given two atomic numbers and a bond type, what is the score
    // associated with `distance`.
    float Score(int at1, int at2, int btype, float distance) const;
};

int
EmpiricalBondLengthDistributions::Index(int at1, int at2, int btype) const {
  if (at1 > at2) {
    std::swap(at1, at2);
  }

  return 100 * btype + 10 * at1 + at2;
}

int
EmpiricalBondLengthDistributions::Build(const BondLengthDistribution::Distributions& proto) {
  for (const auto& dist : proto.distribution()) {
    int at1 = dist.atomic_number1();
    int at2 = dist.atomic_number2();
    int btype = dist.btype();
    int ndx = Index(at1, at2, btype);
    //auto [iter, _] = _dist.emplace(std::make_pair(ndx, EmpiricalBondLengthDistribution()));
    auto [iter, _] = _dist.emplace(ndx, EmpiricalBondLengthDistribution());
    if (! std::get<1>(*iter).Build(dist)) {
      cerr << "EmpiricalBondLengthDistributions::Build:cannot process " << dist.ShortDebugString() << '\n';
      return 0;
    }
  }

  if (_dist.empty()) {
    cerr << "EmpiricalBondLengthDistributions::Build:no data\n";
    return 0;
  }

  cerr << "Read " << _dist.size() << " bond length distributions\n";

  return _dist.size();
}

float
EmpiricalBondLengthDistributions::LongestPlausibleBond() const {
  float result = 0.0;
  for (const auto& [k, v]: _dist) {
    float d = v.max_distance();
    if (d > result) {
      result = d;
    }
  }

  return result;
}

template <typename T>
int
FirstNonZero(const T* values, int nvalues) {
  if (values[0] > 0) {
    return 0;
  }
  for (int i = 1; i < nvalues; ++i) {
    if (values[i] > 0) {
      return i;
    }
  }

  return -1;
}

// Interpolate nonzero values in `values`.
// `prev_nonzero_index` was nonzero, and `to` is nonzer.
// Interpolate the values in between.
template <typename T>
void
Interpolate(T* values, int prev_nonzero_index, int to) {
  float lhs = values[prev_nonzero_index];
  float rhs = values[to];
  float dx = (rhs - lhs) / static_cast<float>(to - prev_nonzero_index);
  int ndx = prev_nonzero_index + 1;
  for (int i = prev_nonzero_index + 1; i < to ; ++i, ++ndx) {
    values[i] = static_cast<T>(ndx * dx);
  }
}

template <typename T>
void
InterpolateZeros(T* values, int nvalues) {
  int prev_nonzero_index = FirstNonZero(values, nvalues);
  // This should never happen.
  if (prev_nonzero_index < 0) {
    return;
  }
  bool prev_was_zero = false;
  for (int i = prev_nonzero_index; i < nvalues; ++i) {
    if (values[i] == 0) {
      prev_was_zero = true;
      continue;
    }
    if (prev_was_zero) {
      Interpolate(values, prev_nonzero_index, i);
    }
    prev_nonzero_index = i;
    prev_was_zero = false;
  }
}

int
EmpiricalBondLengthDistribution::Build(const BondLengthDistribution::Distribution& bld) {
  _atomic_number1 = bld.atomic_number1();
  _atomic_number2 = bld.atomic_number2();
  switch (bld.btype()) {
    case BondLengthDistribution::SS_SINGLE_BOND:
      _btype = 1;
      break;
    case BondLengthDistribution::SS_DOUBLE_BOND:
      _btype = 2;
      break;
    case BondLengthDistribution::SS_TRIPLE_BOND:
      _btype = 3;
      break;
    default:
      cerr << "EmpiricalBondLengthDistribution::Build:unrecognised bond type " << bld.btype() << '\n';
      return 0;
  }

  _min_distance = bld.min_distance();
  _dx = bld.dx();
  if (_dx <= 0.0) {
    cerr << "EmpiricalBondLengthDistribution::Build:invalid dx " << _dx << '\n';
    return 0;
  }

  uint32_t ndistances = bld.count_size();
  if (ndistances < 1) {
    cerr << "EmpiricalBondLengthDistribution::Build:no count data\n";
    return 0;
  }

  _distance = new float[ndistances];
  _count = new uint32_t[ndistances];
  uint32_t total_count = 0;
  for (uint32_t i = 0; i < ndistances; ++i) {
    _distance[i] = _min_distance + _dx * i;
    _count[i] = bld.count(i);
    total_count += _count[i];
  }

  _max_distance = _min_distance + ndistances * _dx;
  cerr << "min " << _min_distance << " ndistances " << ndistances << " dx " << _dx << " max " << _max_distance << '\n';

  if (ndistances > 1) {
    InterpolateZeros(_count, ndistances);
  }

  _score = new float[ndistances];
  // load _score with the normalised distribution.
  for (uint32_t i = 0; i < ndistances; ++i) {
    _score[i] = static_cast<float>(_count[i]) / static_cast<float>(total_count);
  }

  return 1;
}

float
EmpiricalBondLengthDistributions::Score(int at1, int at2, int btype, float distance) const {
  const int ndx = Index(at1, at2, btype);
  auto iter = _dist.find(ndx);
  if (iter == _dist.end()) {
    // cerr << "EmpiricalBondLengthDistributions::Score:no data for " << at1 << ' ' << at2 << " btype " << btype << '\n';
    return 0.0;
  }
  return iter->second.Score(distance);
}

struct Parameters {
  float _sp2_tolerance;
  float _max_possible_bond_length;

  int Initialise(Command_Line& cl);
};

struct BtypeScore {
  int btype;
  float score;

  BtypeScore(int bt, float s) {
    btype = bt;
    score = s;
  }
};

// For each bond that can be varied, the two atoms, and the score
// associated with the plausible bonds.
struct ChangeableBond {
  int a1;
  int a2;
  std::vector<BtypeScore> score;

  public:
    ChangeableBond();
};

ChangeableBond::ChangeableBond() {
  a1 = -1;
  a2 = -1;
}

// Each atom will have a vector of AtomAndScores, describing
// the atoms to which it might be bonded, and the scores associated
// with each possible bond type.
struct AtomAndScores {
  int atom;
  std::vector<BtypeScore> score;
}

// A structure used for sorting neighbours.
struct IDDist{
  int id;
  float dist;
};

// There are many working arrays associated with the current molecule.
class CurrentMoleculeData {
  private:
    Molecule& _mol;
    const int _matoms;

    float* _distance_matrix;

    // All the other atoms that are within range of each atom.
    std::vector<Set_of_Atoms> _nbrs;

    // Whether or not an atom is likely sp2.
    int * _sp2;

  public:
    CurrentMoleculeData(Molecule& m);
    ~CurrentMoleculeData();
};

CurrentMoleculeData::CurrentMoleculeData(Molecule& m) : _mol(m), _matoms(m.natoms()) {

  _distance_matrix = new float[_matoms * _matoms];
  for (int i = 0; i < _matoms; ++i) {
    for (int j = i + 1; j < _matoms; ++j) {
      float d = m.distance_between_atoms(i, j);
      _distance_matrix[i * _matoms + j] = d;
      _distance_matrix[j * _matoms + i] = d;
    }
  }
}

#ifdef MAYBEDOTHIS
int
CurrentMoleculeData::Initialise() {
  _sp2 = new int[_matoms];
  for (int i = 0; i < _matoms; ++i) {
  }
}
#endif

class TopologyFromGeometry {
  private:
    int _verbose;

    int _molecules_processed;

    int _reduce_to_largest_fragment;

    float _longest_plausible_bond;

    float _sp2_tolerance;

    EmpiricalBondLengthDistributions _bld;

    FileType _input_type;

  // private functions.

    int UnivalentToNearest(Molecule& m, const float* distance_matrix);
    std::unique_ptr<Set_of_Atoms[]> GatherPlausibleNeighbours(Molecule& m,
                const float* distance_matrix);
    int GatherPlausibleNeighbours(Molecule& m,
                atom_number_t zatom,
                const float* distance_matrix,
                IDDist* id_dist,
                Set_of_Atoms& result);
    std::unique_ptr<int[]> DiscernSp2(Molecule& m,
                                const Set_of_Atoms* nbrs) const;
    int CountSp2(const Molecule& m, atom_number_t zatom, const Set_of_Atoms& nbrs) const;
    std::unique_ptr<int[]> DetermineBondsNeeded(Molecule& m) const;
    std::vector<float[4]> GatherPlausibleBondTypes(Molecule& m,
                        const Set_of_Atoms* nbrs,
                        const int * bond,
                        const float* distance_matrix) const;
    std::vector<ChangeableBond> GatherChangeableBonds(Molecule& m,
                        const Set_of_Atoms* nbrs,
                        const int * bond,
                        const float* distance_matrix) const;

  public:
    TopologyFromGeometry();

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);

    int Process(Molecule& m);

    int Report(std::ostream& output) const;

    FileType input_type() const {
      return _input_type;
    }
};

TopologyFromGeometry::TopologyFromGeometry() {
  _verbose = 0;
  _molecules_processed = 0;
  _longest_plausible_bond = 0.0;
  _sp2_tolerance = 3.0;
  _input_type = FILE_TYPE_INVALID;
}

int
TopologyFromGeometry::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (! cl.option_present('B')) {
    cerr << "TopologyFromGeometry::Initialise:must specify the distributions via the -B option\n";
    return 0;
  }

  IWString fname = cl.option_value('B');
  std::optional<BondLengthDistribution::Distributions> proto = iwmisc::ReadTextProto<BondLengthDistribution::Distributions>(fname);
  if (proto == std::nullopt) {
    cerr << "Cannot read bond length distributions '" << fname << "'\n";
    return 0;
  }

  if (! _bld.Build(*proto)) {
    cerr << "TopologyFromGeometry::Initialise:cannot construct bond length distributions '" << fname << "'\n";
    return 0;
  }

  _longest_plausible_bond = _bld.LongestPlausibleBond();
  if (_verbose) {
    cerr << "Across " << _bld.size() << " bond length distributions, longest dist " << _longest_plausible_bond << '\n';
  }

  if (1 == cl.number_elements() && 0 == strcmp("-", cl[0])) { // reading a pipe, assume smiles
    _input_type = FILE_TYPE_SMI;
  } else if (!all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot discern all file types, use the -i option\n";
    return 0;
  } else if (!process_input_type(cl, _input_type)) {
    return 0;
  }

  return 1;
}

int
TopologyFromGeometry::Preprocess(Molecule& m) {

  if (m.natoms() == 0) {
    return 0;
  }

  return 1;
}

std::unique_ptr<float[]>
MakeDistanceMatrix(const Molecule& m) {
  const int matoms = m.natoms();
  std::unique_ptr<float[]> result = std::make_unique<float[]>(matoms * matoms);
  for (int i = 0; i < matoms; ++i) {
    for (int j = i + 1; j < matoms; ++j) {
      float d = m.distance_between_atoms(i, j);
      result[i * matoms + j] = d;
      result[j * matoms + i] = d;
    }
  }

  return result;
}

int
TopologyFromGeometry::GatherPlausibleNeighbours(Molecule& m,
                atom_number_t zatom,
                const float* distance_matrix,
                IDDist* id_dist,
                Set_of_Atoms& destination) {
  const int matoms = m.natoms();

  int nbrs = 0;
  for (int i = 0; i < matoms; ++i) {
    if (i == zatom) {
      continue;
    }
    // cerr << " atom " << zatom << " type " << m.atomic_number(zatom) << " to " << i << " z " << m.atomic_number(i) << " dist " << distance_matrix[zatom*matoms+i] << '\n';
    // cerr << " bonded " << m.are_bonded(zatom, i) << '\n';

    if (m.are_bonded(zatom, i)) {
      continue;
    }

    float d = distance_matrix[zatom * matoms + i];
    if (d > _longest_plausible_bond) {
      continue;
    }
    id_dist[nbrs].id = i;
    id_dist[nbrs].dist = d;
    ++nbrs;
  }

  if (nbrs == 0) {
    return 0;
  }

  if (nbrs == 1) {
    destination << id_dist[0].id;
    return 1;
  }

  if (nbrs == 2) {
    if (id_dist[0].dist < id_dist[1].dist) {
      destination << id_dist[0].id;
      destination << id_dist[1].id;
    } else {
      destination << id_dist[1].id;
      destination << id_dist[0].id;
    }
    return 1;
  }

  std::sort(id_dist, id_dist + nbrs,
                [](const IDDist& d1, const IDDist& d2) {
                  return d1.dist < d2.dist;
                }
  );
  for (int i = 0; i < nbrs; ++i) {
    destination << id_dist[i].id;
  }

  return nbrs;
}

// Examine distances to return the atoms that are within plausible
// bonding distance for each atom.
std::unique_ptr<Set_of_Atoms[]>
TopologyFromGeometry::GatherPlausibleNeighbours(Molecule& m,
                const float* distance_matrix) {

  const int matoms = m.natoms();
  std::unique_ptr<Set_of_Atoms[]> result = std::make_unique<Set_of_Atoms[]>(matoms);
  std::unique_ptr<IDDist[]> id_dist = std::make_unique<IDDist[]>(matoms);

  for (int i = 0; i < matoms; ++i) {
    GatherPlausibleNeighbours(m, i, distance_matrix, id_dist.get(), result[i]);
  }

  return result;
}

int
TopologyFromGeometry::CountSp2(const Molecule& m,
                atom_number_t zatom,
                const Set_of_Atoms& nbrs) const {
  constexpr float k120 = 120.0f;

  int rc = 0;

  uint32_t n = nbrs.size();
  for (uint32_t i = 0; i < n; ++i) {
    for (uint32_t j = i + 1; j < n; ++j) {
      const angle_t angle = m.bond_angle(nbrs[i], zatom, nbrs[j], 1);
      if (abs(angle - k120) < _sp2_tolerance) {
        ++rc;
      }
    }
  }

  return rc;
}

std::unique_ptr<int[]>
TopologyFromGeometry::DiscernSp2(Molecule& m,
                                const Set_of_Atoms* nbrs) const {
  const int matoms = m.natoms();
  std::unique_ptr<int[]> result = std::make_unique<int[]>(matoms);
  for (int i = 0; i < matoms; ++i) {
    result[i] = 0;
    uint32_t maybe_bonded = nbrs[i].size();
    if (maybe_bonded < 3) {
      continue;
    }
    // the number of bond angles that seem to be sp2.
    if (CountSp2(m, i, nbrs[i]) == 3) {
      result[i] = 1;
    }
  }

  return result;
}

std::unique_ptr<int[]>
TopologyFromGeometry::DetermineBondsNeeded(Molecule& m) const {
  const int matoms = m.natoms();

  std::unique_ptr<int[]> result = std::make_unique<int[]>(matoms);
  for (int i = 0; i < matoms; ++i) {
    result[i] = 0;
    switch (m.atomic_number(i)) {
      case 1:
        result[i] = 1;
        break;
      case 6:
        result[i] = 4;
        break;
      case 7:
        if (m.formal_charge(i) == 1) {
          result[i] = 4;
        } else {
          result[i] = 3;
        }
        break;
      case 8:
        if (m.formal_charge(i) == -1) {
          result[i] = 1;
        } else {
          result[i] = 2;
        }
        break;
      case 9:
        result[i] = 1;
        break;
      default:
        cerr << "Unhandled atom type " << m.smarts_equivalent_for_atom(i) << '\n';
    }
  }

  return result;
}

std::unique_ptr<int[]>
GatherBonds(const Molecule& m) {
  const int matoms = m.natoms();

  std::unique_ptr<int[]> result(new_int(matoms * matoms));

  for (const Bond* b : m.bond_list()) {
    if (b->is_single_bond()) {
      result[b->a1() * matoms + b->a2()] = 1;
      result[b->a2() * matoms + b->a1()] = 1;
    } else if (b->is_double_bond()) {
      result[b->a1() * matoms + b->a2()] = 2;
      result[b->a2() * matoms + b->a1()] = 2;
    } else if (b->is_triple_bond()) {
      result[b->a1() * matoms + b->a2()] = 3;
      result[b->a2() * matoms + b->a1()] = 3;
    }
  }

  return result;
}

int
TopologyFromGeometry::Process(Molecule& m) {
  ++_molecules_processed;

  const int matoms = m.natoms();

  std::unique_ptr<float[]> distance_matrix = MakeDistanceMatrix(m);

  UnivalentToNearest(m, distance_matrix.get());

  std::unique_ptr<Set_of_Atoms[]> nbrs = GatherPlausibleNeighbours(m, distance_matrix.get());

  std::unique_ptr<int[]> is_sp2 = DiscernSp2(m, nbrs.get());

  std::unique_ptr<int[]> bonds_needed = DetermineBondsNeeded(m);

  std::unique_ptr<int[]> initial_bonded = GatherBonds(m);

  std::vector<ChangeableBond> possible_bonds = GatherChangeableBonds(m, nbrs.get(), initial_bonded.get(), distance_matrix.get());

  if (_verbose) {
    cerr << m.name() << " has " << possible_bonds.size() << " possible bonds " << m.smiles() << '\n';
  }

  if (possible_bonds.empty()) {
    cerr << m.name() << " no changeable bonds\n";
  }

  const uint32_t npb = possible_bonds.size();

  std::vector<int> state(npb);
  for (uint32_t i = 0; i < npb; ++i) {
    state[i] = possible_bonds[i].score.size();
  }

  combinations::Combinations<int> combo(state);
  std::fill(state.begin(), state.end(), 0);
  int combinations_tried = 0;
  while (combo.Next(state)) {
    for (int s : state) {
      cerr << ' ' << s;
    }
    cerr << '\n';
    //std::copy_n(initial_bonded.get(), current_bonded.get(), npb);
    //if (! PlaceBonds(
    ++combinations_tried;
  }
  if (_verbose) {
    cerr << m.name() << " tried " << combinations_tried << " combinations\n";
  }

  return 1;
}

// Attach H and F atoms to their closest neighbour.
// Return a matrix of whether or not pairs of atoms are bonded.
int
TopologyFromGeometry::UnivalentToNearest(Molecule& m,
                        const float* distance_matrix) {
  const int matoms = m.natoms();

  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    const atomic_number_t z = m.atomic_number(i);
    if (z == 1 || z == 9) {
    } else {
      continue;
    }

    int closest = -1;
    float shortest_distance = std::numeric_limits<float>::max();
    for (int j = 0; j < matoms; ++j) {
      if (j == i || m.are_bonded(i, j)) {
        continue;
      }
      float d = distance_matrix[i * matoms + j];
      if (d < shortest_distance) {
        shortest_distance = d;
        closest = j;
      }
    }

    if (m.hcount(closest) == 0) {
      continue;
    }
    m.add_bond(i, closest, SINGLE_BOND);
    ++rc;
  }

  return rc;
}

// For each atom a list of possible bonds.
std::vector<AtomAndScores>
TopologyFromGeometry::GatherChangeableBonds(Molecule& m,
                        const Set_of_Atoms* nbrs,
                        const int * bond,
                        const float* distance_matrix) const {
  const int matoms = m.natoms();
  std::vector<AtomAndScores> result(matoms);
  for (int i = 0; i < matoms; ++i) {
    const atomic_number_t zi = m.atomic_number(i);
    for (int j : nbrs[i]) {
      const atomic_number_t zj = m.atomic_number(j);
      if (bond[i * matoms + j]) {
        continue;
      }
      for (int btype = 1; btype <= 3; ++btype) {
        float s = _bld.Score(zi, zj, btype, distance_matrix[i * matoms + j]);
        if (s > 0.0f) {
          result[i].score.emplace_back(ChangeableBond(i, j, btype, s));
        }
      }
    }
    if (result[i].score.size() > 0) {
      std::sort(result[i].score.begin(), result[i].score.end(),
                  [](const ChangeableBond& c1, const ChangeableBond& c2) {
                    return c1.score < c2.score;
                  }
      );
    }
  }
  return result;
}

#ifdef OLD_VWERSION
// For each pair of atoms, the score of each bond type.
std::vector<ChangeableBond>
TopologyFromGeometry::GatherChangeableBonds(Molecule& m,
                        const Set_of_Atoms* nbrs,
                        const int * bond,
                        const float* distance_matrix) const {
  const int matoms = m.natoms();
  // first count the number of changeable bonds.
  int changeable = 0;
  for (int i = 0; i < matoms; ++i) {
    // cerr << " atom " << i << " type " << m.atomic_number(i) << " has " << nbrs[i].size() << " nbrs\n";
    for (int j : nbrs[i]) {
      if (j < i || bond[i * matoms + j]) {
        continue;
      }
      ++changeable;
    }
  }

  std::vector<ChangeableBond> result(changeable);

  // Index into result.
  int ndx = 0;

  for (int i = 0; i < matoms; ++i) {
    const atomic_number_t zi = m.atomic_number(i);
    for (int j : nbrs[i]) {
      // SKip if already bonded.
      if (j < i || bond[i * matoms + j]) {
        continue;
      }
      const atomic_number_t zj = m.atomic_number(j);
      result[ndx].a1 = zi;
      result[ndx].a2 = m.atomic_number(j);
      const float dist = distance_matrix[i * matoms + j];
      for (int btype = 1; btype <= 3; ++btype) {
        const float score = _bld.Score(zi, zj, btype, dist);
        if (score <= 0.0f) {
          continue;
        }
        result[ndx].score.emplace_back(BtypeScore(btype, score));
      }
      ++ndx;
    }
  }

  assert(ndx == changeable);

  return result;
}
#endif

int
TopologyFromGeometry::Report(std::ostream& output) const {
  cerr << "TopologyFromGeometry:read " << _molecules_processed << " molecules\n";
  if (_molecules_processed == 0) {
    return 1;
  }

  return 1;
}

int
DoTopologyFromGeometry(TopologyFromGeometry& topology_from_geometry,
            Molecule& m,
            IWString_and_File_Descriptor& output) {
  if (! topology_from_geometry.Process(m)) {
    cerr << "Cannot process '" << m.name() << '\n';
    return 0;
  }

  return 1;
}


int
DoTopologyFromGeometry(TopologyFromGeometry& topology_from_geometry,
            data_source_and_type<Molecule>& input,
            IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    m->remove_all_bonds();

    if (! topology_from_geometry.Preprocess(*m)) {
      return 0;
    }

    if (! DoTopologyFromGeometry(topology_from_geometry, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
DoTopologyFromGeometry(TopologyFromGeometry& topology_from_geometry,
            const char * fname,
            IWString_and_File_Descriptor& output) {
  FileType input_type = topology_from_geometry.input_type();
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return DoTopologyFromGeometry(topology_from_geometry, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:i:B:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (! process_standard_aromaticity_options(cl, verbose, 'A')) {
    cerr << "Cannot process aromaticity options\n";
    return 1;
  }

  if (! process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process standard elements options (-E)\n";
    return 1;
  }

  TopologyFromGeometry topology_from_geometry;
  if (! topology_from_geometry.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);
  for (const char* fname : cl) {
    if (! DoTopologyFromGeometry(topology_from_geometry, fname, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    topology_from_geometry.Report(cerr);
  }

  return 0;
}

}  // namespace topology_from_geom

int
main(int argc, char** argv) {
  int rc = topology_from_geom::Main(argc, argv);

  return rc;
}
