#include <string>
#include <vector>

#include "jlcxx/jlcxx.hpp"

#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"

namespace lillymol_julia {

std::string
greet() {
  return "hello world";
}

struct World
{
  World(const std::string& message = "default hello") : msg(message){}
  void set(const std::string& msg) { this->msg = msg; }
  std::string greet() { return msg; }
  std::string msg;
  ~World() { std::cout << "Destroying World with message " << msg << std::endl; }
};

JLCXX_MODULE define_types_module(jlcxx::Module& types)
{
  types.add_bits<FileType>("FileType", jlcxx::julia_type("CppEnum"));
  types.set_const("SMI", FILE_TYPE_SMI);
  types.set_const("SDF", FILE_TYPE_SDF);

  types.add_type<World>("World")
    .constructor<const std::string&>()
    .method("set", &World::set)
    .method("greet", &World::greet)
  ;
}

enum BondType {
  kInvalidBond = 0,
  kSingleBond = 1,
  kDoubleBond = 2,
  kTripleBond = 3,
  kAromaticBond = 4
};

BondType
ToBondType(const Bond& b) {
  if (b.is_aromatic()) {
    return kAromaticBond;
  }
  if (b.is_single_bond()) {
    return kSingleBond;
  }
  if (b.is_double_bond()) {
    return kDoubleBond;
  }
  if (b.is_triple_bond()) {
    return kTripleBond;
  }

  return kInvalidBond;
}

bond_type_t
BtypeEnumToBtype(const BondType btenum) {
  switch (btenum) {
    case kSingleBond:
      return SINGLE_BOND;
    case kDoubleBond:
      return DOUBLE_BOND;
    case kTripleBond:
      return TRIPLE_BOND;
    case kAromaticBond:
      return AROMATIC_BOND;
    default:
      return NOT_A_BOND;
  }
}

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
  mod.method("greet", &greet);

  mod.add_type<World>("World")
    .constructor<const std::string&>()
    .method("set", &World::set)
    .method("greet", &World::greet)
  ;

  mod.add_bits<BondType>("BondType", jlcxx::julia_type("CppEnum"));
  mod.set_const("SINGLE", kSingleBond);
  mod.set_const("DOUBLE", kDoubleBond);
  mod.set_const("TRIPLE", kTripleBond);
  mod.set_const("AROMATIC", kAromaticBond);
    
  mod.add_bits<FileType>("FileType", jlcxx::julia_type("CppEnum"));
  mod.set_const("SMI", FILE_TYPE_SMI);
  mod.set_const("SDF", FILE_TYPE_SDF);

  mod.add_type<Set_of_Atoms>("SetOfAtoms")
    .constructor<>()
    .method("contains", &Set_of_Atoms::contains)
  ;

  mod.add_type<Ring>("Ring")
    .constructor<>()
    .method("ring_number", &Ring::ring_number)
    .method("fragment_membership", &Ring::fragment_membership)
    .method("fused_system_identifier", &Ring::fused_system_identifier)
    .method("fused_ring_neighbours", &Ring::fused_ring_neighbours)
    .method("fused_neighbour", &Ring::fused_neighbour)
    .method("largest_number_of_bonds_shared_with_another_ring", &Ring::largest_number_of_bonds_shared_with_another_ring)
    .method("strongly_fused_ring_neighbours", &Ring::strongly_fused_ring_neighbours)
    .method("contains_bond", &Ring::contains_bond)
    .method("contains_both", &Ring::contains_both)

    .method("is_fused",
      [](const Ring& r)->bool{
        return r.is_fused();
      }
    )
    .method("is_aromatic",
      [](const Ring& r)->bool{
        return r.is_aromatic();
      }
    )
  ;

  mod.add_type<Chemical_Standardisation>("ChemicalStandardisation")
    .method("activate_all", &Chemical_Standardisation::activate_all)
    .method("process", 
      [](Chemical_Standardisation& s, jlcxx::BoxedValue<Molecule>& boxed_mol)->int {
        Molecule& m = jlcxx::unbox<Molecule&>(boxed_mol);
        return s.process(m);
      }
    )
  ;
    
  mod.add_type<Bond>("Bond")
    .method("a1", &Bond::a1)
    .method("a2", &Bond::a2)
    .method("btype",
      [](const Bond& b)->BondType{
        return ToBondType(b);
      }
    )
    .method("other", &Bond::other)
    .method("is_single_bond",
      [](const Bond& b)->bool{
        return b.is_single_bond();
      }
    )
    .method("is_double_bond",
      [](const Bond& b)->bool{
        return b.is_double_bond();
      }
    )
    .method("is_triple_bond",
      [](const Bond& b)->bool{
        return b.is_triple_bond();
      }
    )
    .method("is_aromatic",
      [](const Bond& b)->bool{
        return b.is_aromatic();
      }
    )
    .method("is_aromatic_bond",
      [](const Bond& b)->bool{
        return b.is_aromatic();
      }
    )
    .method("nrings",
      [](const Bond& b) {
        return b.nrings();
      }
    )
    .method("bond_number", &Bond::bond_number)
    .method("either_atom_set_in_array", &Bond::either_atom_set_in_array)
    .method("atoms",
      [](const Bond& b)->std::tuple<atom_number_t, atom_number_t>{
        return std::make_tuple(b.a1(), b.a2());
      }
    )
    .method("involves",
      [](const Bond& a, atom_number_t o)->bool{
        return a.involves(o);
      }
    )
    .method("involves",
      [](const Bond& a, atom_number_t o1, atom_number_t o2)->bool{
        return a.involves(o1, o2);
      }
    )
    .method("joins",
      [](const Bond& b1, const Bond& b2)->bool{
        return b1.joins(&b2);
      }
    )
  ;

  mod.add_type<Atom>("Atom")
    // .constructor<jlcxx::cxxint_t>(true) // with finalizer
    .method("valence_ok", &Atom::valence_ok)
    .method("ncon",
      [](const Atom& a) {
        return a.ncon();
      }
    )
    .method("isotope", &Atom::isotope)
    .method("lone_pair_count",
      [](Atom&a){
        int lp;
        if (a.lone_pair_count(lp)) {
          return lp;
        }
        return 0;
      }
    )
    .method("nbonds",
      [](const Atom& a) {
        return a.nbonds();
      }
    )
    .method("nbonds",
      [](Atom& a) {
        return a.nbonds();
      }
    )
    .method("bond_to_atom",
      [](const Atom& a, atom_number_t o)->const Bond*{
        return a.bond_to_atom(o);
      }
    )
    .method("is_bonded_to",
      [](const Atom& a, atom_number_t o)->bool{
        return a.is_bonded_to(o);
      }
    )
    .method("is_halogen", &Atom::is_halogen)
    .method("atomic_number", &Atom::atomic_number)
    .method("atomic_symbol",
      [](const Atom& a)->std::string{
        return a.atomic_symbol().AsString();
      }
    )
    .method("implicit_hydrogens", &Atom::implicit_hydrogens)
    .method("formal_charge", &Atom::formal_charge)
    .method("other", &Atom::other)
    .method("connections",
      [](const Atom& a, atom_number_t me)->Set_of_Atoms{
        return a.connections(me);
      }
    )
    .method("atomic_weight", &Atom::atomic_weight)
    .method("saturated", 
      [](const Atom& a)->bool{
        return a.fully_saturated();
      }
    );

    mod.set_override_module(jl_base_module);
    mod.method("getindex",
      [](const Atom& a, int i)->const Bond&{
        return *a[i-1];
      }
    );
    mod.unset_override_module();
  ;

  mod.add_type<Molecule>("Molecule")
    .constructor<>()
    //.constructor<jlcxx::cxxint_t>(false) // no finalizer
    .method("ok",
      [](const Molecule& m)->bool{
        return m.ok();
      }
    )
    .method("debug_string", &Molecule::debug_string)
    .method("empty",
      [](const Molecule& m)->bool{
        return m.empty();
      }
    )
    .method("x", &Molecule::x)
    .method("y", &Molecule::y)
    .method("z", &Molecule::z)
    .method("setx", &Molecule::setx)
    .method("sety", &Molecule::sety)
    .method("setz", &Molecule::setz)

    .method("add_bond",
      [](Molecule& m, atom_number_t a1, atom_number_t a2, BondType bt)->bool{
        return m.add_bond(a1, a2, BtypeEnumToBtype(bt));
      }
    )
    .method("are_bonded",
      [](const Molecule& m, atom_number_t a1, atom_number_t a2)->bool{
        return m.are_bonded(a1, a2);
      }
    )
    .method("are_adjacent",
      [](const Molecule& m, atom_number_t a1, atom_number_t a2)->bool{
        return m.are_adjacent(a1, a2);
      }
    )

    .method("has_formal_charges",
      [](const Molecule& m)->bool {
        return m.has_formal_charges();
      }
    )
    .method("number_formally_charged_atoms", &Molecule::number_formally_charged_atoms)
    .method("net_formal_charge", &Molecule::net_formal_charge)


    .method("set_name",
      [](Molecule& m, const std::string& s) {
        return m.set_name(s);
      }
    )
    .method("name",
      [](const Molecule& m)->std::string{
        return m.name().AsString();
      }
    )
    .method("natoms",
      [](const Molecule& m)->int{
        return m.natoms();
      }
    )
    .method("natoms",
      [](const Molecule& m, atomic_number_t z)->int {
        return m.natoms(z);
      }
    )
    .method("natoms",
      [](const Molecule& m, std::string& asymbol)->int{
        return m.natoms(asymbol.c_str());
      }
    )
    .method("nedges", &Molecule::nedges)
    .method("atomic_number", &Molecule::atomic_number)
    .method("nrings",
      [](Molecule& m) {
        return m.nrings();
      }
    )
    .method("nrings",
      [](Molecule& m, atom_number_t a) {
        return m.nrings(a);
      }
    )
    .method("nrings",
      [](Molecule& m, atom_number_t a, int rsize) {
        return m.nrings(a, rsize);
      }
    )
    .method("is_ring_atom",
      [](Molecule& m, atom_number_t a)->bool{
        return m.is_ring_atom(a);
      }
    )
    .method("ring_bond_count", 
      [](Molecule& m, atom_number_t a) {
        return m.ring_bond_count(a);
      }
    )
    .method("fused_system_size", &Molecule::fused_system_size)
    .method("rings_with_fused_system_identifier", &Molecule::rings_with_fused_system_identifier)
    .method("fused_system_identifier", &Molecule::fused_system_identifier)
    .method("in_same_ring", &Molecule::in_same_ring)
    .method("in_same_ring_system", &Molecule::in_same_ring_system)
    .method("ring_membership",
      [](Molecule& m)->std::vector<int>{
        const int* r = m.ring_membership();
        return std::vector<int>(r, r + m.natoms());
      }
    )
    .method("rings_containing_both",
      [](Molecule& m, atom_number_t a1, atom_number_t a2) {
        return m.in_same_rings(a1, a2);
      }
    )
    .method("is_part_of_fused_ring_system",
      [](Molecule& m, atom_number_t a)->bool{
        return m.is_part_of_fused_ring_system(a);
      }
    )
    .method("ring",
      [](Molecule& m, int rnum)->const Ring*{
        return m.ringi(rnum);
      }
    )
    .method("ring_containing_atom", &Molecule::ring_containing_atom)
    .method("sssr_rings",
      [](Molecule& m)->std::vector<const Ring*>{
        std::vector<const Ring*> result;
        result.reserve(m.nrings());
        for (const Ring* r : m.sssr_rings()) {
          result.push_back(r);
        }
        return result;
      }
    )
    .method("label_atoms_by_ring_system",
      [](Molecule& m)->std::vector<int>{
        const int matoms = m.natoms();
        std::unique_ptr<int[]> tmp = std::make_unique<int[]>(matoms);
        m.label_atoms_by_ring_system(tmp.get());
        return std::vector<int>(tmp.get(), tmp.get() + matoms);
      }
    )
    .method("label_atoms_by_ring_system_including_spiro_fused",
      [](Molecule& m)->std::vector<int>{
        const int matoms = m.natoms();
        std::unique_ptr<int[]> tmp = std::make_unique<int[]>(matoms);
        m.label_atoms_by_ring_system_including_spiro_fused(tmp.get());
        return std::vector<int>(tmp.get(), tmp.get() + matoms);
      }
    )
    .method("nrings_including_non_sssr_rings", &Molecule::nrings_including_non_sssr_rings)
    .method("non_sssr_rings", &Molecule::non_sssr_rings)
    .method("non_sssr_ring", &Molecule::non_sssr_ring)
    .method("is_spiro_fused",
      [](Molecule& m, atom_number_t a)->bool{
        return m.is_spiro_fused(a);
      }
    )
    .method("is_halogen", &Molecule::is_halogen)
    .method("ncon",
      [](const Molecule& m, atom_number_t a) {
        return m.ncon(a);
      }
    )
    .method("nbonds",
      [](const Molecule& m, atom_number_t a) {
        return m.nbonds(a);
      }
    )
    .method("maximum_connectivity", &Molecule::maximum_connectivity)
    .method("other", &Molecule::other)
    .method("connections",
      [](const Molecule& m, atom_number_t a)->std::vector<atom_number_t>{
        std::vector<atom_number_t> result;
        for (atom_number_t o : m.connections(a)) {
          result.push_back(o);
        }
        return result;
      }
    )

    .method("build_from_smiles", static_cast<int (Molecule::*)(const std::string& smiles)>(&Molecule::build_from_smiles))
    .method("smiles", 
      [](Molecule& m)->std::string{
        return m.smiles().AsString();
      }
    )
    .method("unique_smiles", 
      [](Molecule& m)->std::string{
        return m.unique_smiles().AsString();
      }
    )
    .method("random_smiles", 
      [](Molecule& m)->std::string{
        return m.random_smiles().AsString();
      }
    )
    .method("unique_kekule_smiles", 
      [](Molecule& m)->std::string{
        return m.UniqueKekuleSmiles().AsString();
      }
    )
    .method("isotopically_labelled_smiles",
      [](Molecule& m)->std::string{
        return m.isotopically_labelled_smiles().AsString();
      }
    )
    .method("is_aromatic", &Molecule::is_aromatic)
    .method("atom", &Molecule::atom)
  ;

  mod.set_override_module(jl_base_module);
  mod.method("getindex",
    [](const Molecule& m, int i)->const Atom&{
      return m[i-1];
    }
  );
  mod.unset_override_module();
}

}  // namespace lillymol_julia
