#include <stdlib.h>

#include "Foundational/iwbits/iwbits.h"
#include "Molecule_Lib/llymol_lemon.h"

namespace llymol_lemon {

NodesAndEdges::NodesAndEdges() {
  _natoms = 0;
  _node = nullptr;
  _edge = nullptr;
}

NodesAndEdges::~NodesAndEdges() {
  if (_node != nullptr) {
    delete [] _node;
  }
  if (_edge != nullptr) {
    delete [] _edge;
  }
}

// The only two properties that work for a substructure search
// are atomic number and aromaticity.
// But can we set up something so that a less than relationship
// will work for comparison.
// For example if ring_bond_count is 3 in the query, then
// anything with 2 in the molecule being searched cannot match,
// but the other way around is fine.
// Aromatic might be a little tricky...

// For an atom return a single 32 bit unsigned int that encodes
//   atomic number 8 bits.
//   ring bond count 3 bits : 11
//   aromatic 1 bit : 12
//   ncon 3 bits : 15
//   hcount 3 bits : 18
//   unsaturation 2 bits : 20
uint32_t
ComputeLabel(Molecule& m, atom_number_t zatom) {
  const Atom& a = m[zatom];

  uint32_t rc = a.atomic_number() << (32 - 8);

  uint32_t tmp = m.ring_bond_count(zatom);
  if (tmp) {
    rc |= tmp << (32 - 11);
    if (m.is_aromatic(zatom)) {
      rc |= one_bit_32[32 - 12];
    }
  }

  tmp = a.ncon();
  rc |= (tmp << one_bit_32[32 - 15]);

  return rc;
}

int
NodesAndEdges::Build(Molecule& m) {
  _natoms = m.natoms();

  if (_natoms == 0) {
    return 1;
  }

  _node = new uint32_t[_natoms];
  _edge = new uint32_t[m.nedges()];

  m.compute_aromaticity_if_needed();

  for (int i = 0; i < _natoms; ++i) {
    _node[i] = ComputeLabel(m, i);
  }

  return 1;
}

int
BuildLemonGraph(Molecule& m,
                lemon::ListGraph& result,
                lemon::ListGraph::Node* nodes) {
  const int matoms = m.natoms();
  result.reserveNode(matoms);
  result.reserveEdge(m.nedges());

  for (int i = 0; i < matoms; ++i) {
    nodes[i] = result.addNode();
  }

  for (const Bond* b : m.bond_list()) {
    result.addEdge(nodes[b->a1()], nodes[b->a2()]);
  }

  return 1;
}

NodeValue
GetLabel(Molecule& m,
         atom_number_t zatom) {
  NodeValue rc;
  rc.atomic_number = m.atomic_number(zatom);
  if (m.is_aromatic(zatom)) {
    rc.aromaticity = NodeValue::Aromatic::kAromatic;
  } else {
    rc.aromaticity = NodeValue::Aromatic::kAliphatic;
  }

  return rc;
}

EdgeValue
GetLabel(Molecule& m,
         const Bond* b) {
  EdgeValue rc;
  if (b->is_aromatic()) {
    rc.btype = EdgeValue::Btype::kAromatic;
  } else if (b->is_single_bond()) {
    rc.btype = EdgeValue::Btype::kSingle;
  } else if (b->is_double_bond()) {
    rc.btype = EdgeValue::Btype::kDouble;
  } else if (b->is_triple_bond()) {
    rc.btype = EdgeValue::Btype::kTriple;
  }

  return rc;
}

int
AttachMap(Molecule& m, NodeValue* node_value, EdgeValue* edge_value,
                const lemon::ListGraph::Node* nodes,
                lemon::ListGraph& result) {
  lemon::ListGraph::NodeMap<NodeValue> atom_label(result);
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    atom_label[nodes[i]] = GetLabel(m, i);
  }

  lemon::ListGraph::ArcMap<EdgeValue> bond_label(result);
  for (int ndx = 0; lemon::ListGraph::ArcIt e(result); e != INVALID; ++e) {
    bond_label->set(e, GetLabel(m, m.bondi(ndx));
    ++ndx;
  }

  return 1;
}

};  // namespace llymol_lemon
