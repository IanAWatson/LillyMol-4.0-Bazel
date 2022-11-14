#ifndef MOLECULE_TOOLS_TOPOLOGY_FROM_GEOMETRY_EMPIRICAL_BOND_LENGTH_DISTRIBUTION_H
#define MOLECULE_TOOLS_TOPOLOGY_FROM_GEOMETRY_EMPIRICAL_BOND_LENGTH_DISTRIBUTION_H

#include "Molecule_Tools/topology_from_geometry/bond_length_distribution.pb.h"

namespace topology_from_geometry {

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
    float LongestPlausibleBond(int z1, int z2) const;

    // The number of different bond types we have.
    uint32_t size() const {
      return _dist.size();
    }

    // Given two atomic numbers and a bond type, what is the score
    // associated with `distance`.
    float Score(int at1, int at2, int btype, float distance) const;

    uint32_t Mask(int at1, int at2, float distance) const;
};

}  // namespace topology_from_geometry

#endif // MOLECULE_TOOLS_TOPOLOGY_FROM_GEOMETRY_EMPIRICAL_BOND_LENGTH_DISTRIBUTION_H
