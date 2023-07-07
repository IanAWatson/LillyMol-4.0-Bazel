#include <string>

#include "jlcxx/jlcxx.hpp"

#include "Molecule_Lib/molecule.h"

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

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
  mod.method("greet", &greet);

  mod.add_type<World>("World")
    .constructor<const std::string&>()
    .method("set", &World::set)
    .method("greet", &World::greet)
  ;
    
  mod.add_bits<FileType>("FileType", jlcxx::julia_type("CppEnum"));
  mod.set_const("SMI", FILE_TYPE_SMI);
  mod.set_const("SDF", FILE_TYPE_SDF);

  mod.add_type<Molecule>("Molecule")
    .constructor<>()
    .constructor<jlcxx::cxxint_t>(false) // no finalizer
    .method("ok",
      [](const Molecule& m)->bool{
        return m.ok();
      }
    )
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

    .method("build_from_smiles", static_cast<int (Molecule::*)(const std::string& smiles)>(&Molecule::build_from_smiles))
    .method("smiles", 
      [](Molecule& m)->std::string{
        return m.smiles().AsString();
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
  ;
}

}  // namespace lillymol_julia
