#include <stdlib.h>

#include <algorithm>
#include <limits>
#include <optional>
#include <unordered_map>

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

int
BuildLemonGraph(Molecule& m,
                lemon::ListGraph& result) {
  const int matoms = m.natoms();
  std::unique_ptr<lemon::ListGraph::Node[]> nodes = std::make_unique<lemon::ListGraph::Node[]>(matoms);

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

constexpr uint8_t kUndefined = std::numeric_limits<uint8_t>::max();
constexpr uint8_t kChargeUndefined = std::numeric_limits<int8_t>::max();

constexpr int kInvaldAtomNumber = -1;

NodeValue::NodeValue() {
  std::fill_n(_properties, kNproperties, kUndefined);

  _charge =  kChargeUndefined;
  
  _atom_number = kInvaldAtomNumber;
}

NodeValue::NodeValue(const NodeValue& rhs) {
  std::copy_n(rhs._properties, kNproperties, _properties);
  _charge = rhs._charge;
  _atom_number = rhs._atom_number;
}

NodeValue&
NodeValue::operator=(const NodeValue& rhs) {
  std::copy_n(rhs._properties, kNproperties, _properties);
  _charge = rhs._charge;
  _atom_number = rhs._atom_number;
  return *this;
}

NodeValue&
NodeValue::operator=(NodeValue&& rhs) {
  std::copy_n(rhs._properties, kNproperties, _properties);
  _charge = rhs._charge;
  _atom_number = rhs._atom_number;
  return *this;
}

bool
NodeValue::operator==(const NodeValue& rhs) const {
  cerr << "Comparison invoked\n";
  if (_properties[kAtomicNumber] != rhs._properties[kAtomicNumber]) {
    return false;
  }

  return true;
}

std::optional<NodeValue>
GetLabel(const Substructure_Atom_Specifier& a) {

  const resizable_array<const Element*>& element = a.element();
  if (element.size() > 1) {
    cerr << "GetLabel::cannot process multiple elements\n";
    return std::nullopt;
  }

  NodeValue rc;

  rc.set_property(NodeValue::kAtomicNumber, element[0]->atomic_number());

  if (! a.ncon().is_set()) {
  } else if (const std::optional<int> ncon = a.ncon().IsSingleValue()) {
    rc.set_property(NodeValue::kNcon, *ncon);
  } else {
    return std::nullopt;
  }

  if (! a.ring_bond_count().is_set()) {
  } else if (const std::optional<int> rbc = a.ring_bond_count().IsSingleValue()) {
    rc.set_property(NodeValue::kRingBondCount, *rbc);
  } else {
    return std::nullopt;
  }

  if (! a.hcount().is_set()) {
  } else if (const std::optional<int> hcount = a.hcount().IsSingleValue()) {
    rc.set_property(NodeValue::kHcount, *hcount);
  } else {
    return std::nullopt;
  }

  if (! a.unsaturation().is_set()) {
  } else if (const std::optional<int> unsaturation = a.unsaturation().IsSingleValue()) {
    rc.set_property(NodeValue::kUnsaturation, *unsaturation);
  } else {
    return std::nullopt;
  }

  if (! a.ring_size().is_set()) {
  } else if (const std::optional<int> ring_size = a.ring_size().IsSingleValue()) {
    rc.set_property(NodeValue::kRingSize, *ring_size);
  } else {
    return std::nullopt;
  }

  if (! a.get_nrings().is_set()) {
  } else if (const std::optional<int> nrings = a.get_nrings().IsSingleValue()) {
    rc.set_property(NodeValue::kNrings, *nrings);
  } else {
    return std::nullopt;
  }

  if (! a.attached_heteroatom_count().is_set()) {
  } else if (const std::optional<int> ahc = a.attached_heteroatom_count().IsSingleValue()) {
    rc.set_property(NodeValue::kAttachedHeteroatoms, *ahc);
  } else {
    return std::nullopt;
  }

#ifdef IMPLEMENT_ISOTOPES_SOMETIME
  TODO:ianwatson implement isotope handling. Maybe change to matcher in Substructure_Atom_Specifier
  if (! a.isotope().is_set()) {
  } else if (const std::optional<int> iso = a.isotope().IsSingleValue()) {
    rc.set_property(NodeValue::kIsotope, *iso);
  } else {
    return std::nullopt;
  }
#endif

  return rc;
}

std::optional<NodeValue>
GetLabel(const Substructure_Atom& a) {
  if (a.ncomponents() == 0) {
    const Substructure_Atom_Specifier& spec = a;
    return GetLabel(spec);
  }

  if (a.ncomponents() > 1) {
    cerr << "GetLabel::cannot process Substructure_Atom with multiple components\n";
    return std::nullopt;
  }

  const Substructure_Atom_Specifier* spec = a.component(0);

  return GetLabel(*spec);
}

NodeValue
GetLabel(Molecule& m,
         atom_number_t zatom) {
  NodeValue rc;
  rc.set_property(NodeValue::kAtomicNumber, m.atomic_number(zatom));

  if (m.is_aromatic(zatom)) {
    rc.set_property(NodeValue::kAromatic, 1);
  } else {
    rc.set_property(NodeValue::kAromatic, 0);
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
  } else {
    rc.btype = EdgeValue::Btype::kUnspecified;
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

  lemon::ListGraph::EdgeMap<EdgeValue> bond_label(result);
  int ndx = 0;
  for (lemon::ListGraph::EdgeIt e(result); e != lemon::INVALID; ++e) {
    bond_label[e] = GetLabel(m, m.bondi(ndx));
    ++ndx;
  }

  return 1;
}

int
LemonQuery::Build(Substructure_Query& query) {
  if (query.empty()) {
    cerr << "LillyMolQuery::Build:empty query\n";
    return 0;
  }

  if (query.size() > 1) {
    cerr << "LillyMolQuery::Build:cannot process composite query\n";
    return 0;
  }

  return Build(*query[0]);
}

int
LemonQuery::Build(Single_Substructure_Query& query) {
  query.assign_unique_numbers();

  std::unordered_map<int, lemon::ListGraph::Node> id_to_node;

  // First build the connectivity of the graph.
  for (const Substructure_Atom* r : query.RootAtoms()) {
    lemon::ListGraph::Node a = _graph.addNode();
    id_to_node[r->unique_id()] = a;
    AddChildren(r, a, id_to_node);
  }

  // Now tha the graph is build, form the node and edge maps.

  _node_map = std::make_unique<lemon::ListGraph::NodeMap<NodeValue>>(_graph);

  // Now go back and fill in map properties.

  for (const Substructure_Atom* r : query.RootAtoms()) {

    std::optional<NodeValue> maybe_label = GetLabel(*r);
    if (! maybe_label) {
      return 0;
    }

    const lemon::ListGraph::Node& n = id_to_node[r->unique_id()];

    (*_node_map)[n] = *maybe_label;
  }

  // And the edges.
  _edge_map = std::make_unique<lemon::ListGraph::EdgeMap<EdgeValue>>(_graph);
  for (lemon::ListGraph::EdgeIt n(_graph); n != lemon::INVALID; ++n) {
  }

  return 1;
}

int
LemonQuery::AddChildren(const Substructure_Atom* parent_atom,
                           lemon::ListGraph::Node& parent_node,
                           std::unordered_map<int, lemon::ListGraph::Node>& id_to_node) {
  for (const Substructure_Atom* child : parent_atom->Children()) {
    lemon::ListGraph::Node n = _graph.addNode();
    id_to_node[child->unique_id()] = n;
    _graph.addEdge(parent_node, n);
    AddChildren(child, n, id_to_node);
  }

  return 1;
}

LemonMolecule::LemonMolecule() : _node_map(_graph), _edge_map(_graph) {
  _m = nullptr;
}

int
LemonMolecule::Build(Molecule& m) {
  if (m.empty()) {
    return 0;
  }

  const int matoms = m.natoms();

  std::unique_ptr<lemon::ListGraph::Node[]> nodes = std::make_unique<lemon::ListGraph::Node[]>(matoms);

  BuildLemonGraph(m, _graph, nodes.get());

  for (int i = 0; i < matoms; ++i) {
    _node_map[nodes[i]] = GetLabel(m, i);
  }

  // We assume that the edges are ordered the same way we added them.
  int ndx = 0;
  for (lemon::ListGraph::EdgeIt e(_graph); e != lemon::INVALID; ++e) {
    _edge_map[e] = GetLabel(m, m.bondi(ndx));
    ++ndx;
  }

  _m = &m;

  return 1;
}

};  // namespace llymol_lemon
