#ifndef MOLECULE_TOOLS_TOPOLOGY_FROM_GEOM_H
#define MOLECULE_TOOLS_TOPOLOGY_FROM_GEOM_H

#include <memory>
#include <unordered_map>
#include <vector>

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/molecule.h"

#include "Molecule_Tools/bond_length_distribution.pb.h"

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
    float LongestPlausibleBond(atomic_number_t z1, atomic_number_t z2) const;

    // The number of different bond types we have.
    uint32_t size() const {
      return _dist.size();
    }

    // Given two atomic numbers and a bond type, what is the score
    // associated with `distance`.
    float Score(int at1, int at2, int btype, float distance) const;

    uint32_t Mask(int at1, int at2, float distance) const;
};

struct BtypeScore {
  int btype;
  float score;

  BtypeScore(int bt, float s) {
    btype = bt;
    score = s;
  }
};

// Each atom will have a vector of AtomAndScores, describing
// the atoms to which it might be bonded, and the scores associated
// with each possible bond type.
struct AtomAndScores {
  int atom;
  int sp2;
  int atomic_number;
  std::vector<int> bonds_needed;
  int current_bonds;
  std::vector<BtypeScore> score;
};

// A structure used for sorting neighbours.
struct IDDist{
  int id;
  float dist;
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

}  // namespace topology_from_geometry

#endif // MOLECULE_TOOLS_TOPOLOGY_FROM_GEOM_H
