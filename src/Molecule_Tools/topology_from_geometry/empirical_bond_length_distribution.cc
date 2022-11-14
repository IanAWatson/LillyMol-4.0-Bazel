#include <algorithm>
#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Tools/topology_from_geometry/empirical_bond_length_distribution.h"
#include "Molecule_Tools/topology_from_geometry/bond_length_distribution.pb.h"

namespace topology_from_geometry {

using std::cerr;

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
  // THis is asymmetric due to the downward truncation of the distribution.
  if (distance < _min_distance || distance >= _max_distance) {
    return 0.0;
  }

  const int bucket = static_cast<int>((distance - _min_distance) / _dx + 0.4999f);
  return _score[bucket];
}

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

  // The normalized counts.
  std::unique_ptr<double[]>values = std::make_unique<double[]>(ndistances);
  for (uint32_t i = 0; i < ndistances; ++i) {
    values[i] = iwmisc::Fraction<double>(_count[i], total_count);
  }

  // Cumulative sum.
  std::unique_ptr<uint32_t[]> csum = std::make_unique<uint32_t[]>(ndistances);
  // We need to know where the median is found.
  csum[0] = _count[0];
  const uint32_t c50 = total_count / 2;
  int ndx50 = -1;
  for (uint32_t i = 1; i < ndistances; ++i) {
    csum[i] = csum[i - 1] + _count[i];
    if (ndx50 < 0 && csum[i] >= c50) {
      ndx50 = i;
    }
  }

  for (uint32_t i = 0; i < ndistances; ++i) {
    cerr << "csum[" << i << "] " << csum[i] << '\n';
  }
  cerr << "ndx50 " << ndx50 << " c50 " << c50 << " total_count " << total_count << '\n';

  _score = new float[ndistances];
  for (uint32_t i = 0; i <= ndx50; ++i) {
    _score[i] = iwmisc::Fraction<float>(csum[i], c50);
  }

  // Re-use csum
  csum[ndistances - 1] = _count[ndistances - 1];
  for (uint32_t i = ndistances - 2; i > ndx50; --i) {
    csum[i] = csum[i + 1] + _count[i];
  }

  for (uint32_t i = ndx50 + 1; i < ndistances; ++i) {
    _score[i] = iwmisc::Fraction<float>(csum[i], c50);
  }

  cerr << "Before normalisation\n";
  for (uint32_t i = 0; i < ndistances; ++i) {
    cerr << " _score[" << i << "] " << _score[i] << '\n';
  }
  // Normalize the scores.

  float sum_scores = std::accumulate(_score, _score + ndistances, 0.0f);
  for (uint32_t i = 0; i < ndistances; ++i) {
    _score[i] /= sum_scores;
  }

  for (uint32_t i = 0; i < ndistances; ++i) {
    cerr << " _score[" << i << "] " << _score[i] << '\n';
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


}   // namespace topology_from_geometry
