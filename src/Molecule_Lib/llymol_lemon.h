#ifndef MOLECULE_LIB_LLYMOL_LEMON_H
#define MOLECULE_LIB_LLYMOL_LEMON_H

#include <memory>

#include "lemon/list_graph.h"

#include "Molecule_Lib/molecule.h"

namespace llymol_lemon {


// A class to hold node and edge constants for lemon.
class NodesAndEdges {
  private:
    // The number of atoms in the molecule.
    int _natoms;

    // Node and edge constants.
    uint32_t* _node;
    uint32_t* _edge;

  public:
    NodesAndEdges();
    ~NodesAndEdges();

    int Build(Molecule& m);

    uint32_t* NodeLabel() const {
      return _node;
    }
    uint32_t* EdgeLabel() const {
      return _edge;
    }
};

struct NodeValue {
  uint32_t ordinal = 0;
  uint8_t atomic_number = 0;
  enum class Aromatic {
    kUnspecified,
    kAromatic,
    kAliphatic
  };
  Aromatic aromaticity;
};

struct EdgeValue {
  enum class Btype {
    kUnspecified,
    kSingle,
    kDouble,
    kTriple,
    kAromatic
  };
  Btype btype;
};

int
BuildLemonGraph(Molecule& m,
                NodesAndEdges& nodes_and_edges,
                lemon::ListGraph& result);

}  // namespace llymol_lemon

#endif // MOLECULE_LIB_LLYMOL_LEMON_H
