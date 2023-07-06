#ifndef MOLECULE_LIB_LLYMOL_LEMON_H
#define MOLECULE_LIB_LLYMOL_LEMON_H

#include <memory>

#include "lemon/list_graph.h"
#include "lemon/vf2pp.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/substructure.h"

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

class NodeValue {
  public:
    static constexpr uint32_t kAtomicNumber = 0;
    static constexpr uint32_t kAromatic = 1;
    static constexpr uint32_t kNcon = 2;
    static constexpr uint32_t kRingBondCount = 3;
    static constexpr uint32_t kHcount = 4;
    static constexpr uint32_t kUnsaturation = 5;
    static constexpr uint32_t kRingSize = 6;
    static constexpr uint32_t kNrings = 7;
    static constexpr uint32_t kAttachedHeteroatoms = 8;
    static constexpr uint32_t kIsotope = 9;

    static constexpr uint32_t kNproperties = kIsotope + 1;
  private:

    uint32_t _properties[kNproperties];

    int8_t _charge;

    int _atom_number;

  public:
    NodeValue();
    NodeValue(const NodeValue& rhs);
    NodeValue& operator=(const NodeValue& rhs);
    NodeValue& operator=(NodeValue&& rhs);
    bool operator==(const NodeValue& rhs) const;
    bool operator<(const NodeValue& rhs) const;

    void set_atom_number(int s) {
      _atom_number = s;
    }
    int atom_number() const {
      return _atom_number;
    }
    void set_property(uint32_t ndx, uint8_t value) {
      _properties[ndx] = value;
    }
    void set_charge(int8_t chg) {
      _charge = chg;
    }
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

class LemonMolecule {
  private:
    // The molecule being represented. Must remain in scope elsewhere.
    Molecule* _m;

    lemon::ListGraph _graph;

    //std::unique_ptr<lemon::ListGraph::NodeMap<NodeValue>> _node_map;
    //std::unique_ptr<lemon::ListGraph::EdgeMap<EdgeValue>> _edge_map;

    lemon::ListGraph::NodeMap<NodeValue> _node_map;
    lemon::ListGraph::EdgeMap<EdgeValue> _edge_map;

  public:
    LemonMolecule();

    int Build(Molecule& m);

    const lemon::ListGraph& graph() const {
      return _graph;
    }

    const lemon::ListGraph::NodeMap<NodeValue>& node_map() const {
      return _node_map;
    }
    const lemon::ListGraph::EdgeMap<EdgeValue>& edge_map() const {
      return _edge_map;
    }
};

class LemonQuery {
  private:
    Substructure_Query _query;

    lemon::ListGraph _graph;

    std::unique_ptr<lemon::ListGraph::NodeMap<NodeValue>> _node_map;
    std::unique_ptr<lemon::ListGraph::EdgeMap<EdgeValue>> _edge_map;


  // private functions.
    int AddChildren(const Substructure_Atom* parent_atom,
                    lemon::ListGraph::Node& parent_node,
                    std::unordered_map<int, lemon::ListGraph::Node>& id_to_node);

  public:
    int Build(Substructure_Query& query);
    int Build(Single_Substructure_Query& query);

    const lemon::ListGraph& graph() const {
      return _graph;
    }

    const lemon::ListGraph::NodeMap<NodeValue>& node_map() const {
      return *_node_map;
    }
    const lemon::ListGraph::EdgeMap<EdgeValue>& edge_map() const {
      return *_edge_map;
    }
};

int
BuildLemonGraph(Molecule& m,
                NodesAndEdges& nodes_and_edges,
                lemon::ListGraph& result);

int BuildLemonGraph(Molecule& m, lemon::ListGraph& result);

}  // namespace llymol_lemon

#endif // MOLECULE_LIB_LLYMOL_LEMON_H
