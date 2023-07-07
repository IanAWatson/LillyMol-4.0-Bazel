module LillyMol
  using CxxWrap
  @wrapmodule(joinpath("bazel-bin/julia/","lillymol_julia.so"))

  function __init__()
    @initcxx
  end
  export World

  export FileType
  export BondType
  export greet
  export Molecule
  import Base: getindex, iterate, in
  getindex(m::Molecule, a::Int)=atom(m, a)
  iterate(m::Molecule, state=1) = (state >= natoms(m) ? nothing : (m[state], state + 1))
  iterate(a::Atom) = (a, 1)
  iterate(a::Atom, state) = (state == ncon(a) ? nothing : (a, state + 1))
  in(z::Int, m::Molecule) = (natoms(m, z) > 0)
  in(atom::Int, a::Atom) = involves(a, atom)
  length(m::Molecule) = natoms(m)
  size(m::Molecule) = natoms(m)
  export getindex
  export iterate
  export length
  export in
  export natoms, smiles, unique_smiles, nrings, build_from_smiles, is_aromatic
  export atomic_number
  export ncon, nbonds, involves
  export atomic_number
end
#getindex(m::LillyMol.Molecule, a::Int64)=LillyMol.atom(m, a)

using .LillyMol

# Call greet and show the result
@show LillyMol.greet()

w = LillyMol.World()
println(LillyMol.greet(w));
LillyMol.set(w, "hello")
println(LillyMol.greet(w))

m = Molecule()
build_from_smiles(m, "c1ccccc1OC foo")
build_from_smiles(m, "CNOF foo")
println("molecule has $(natoms(m)) atoms")
print(LillyMol.smiles(m))
println(LillyMol.name(m))

for (index, atom) in enumerate(m)
  println("$(index) is $(atom)")
  println("Atom $(index) is $(atomic_number(atom)) with $(ncon(atom)) connections, aromatic $(is_aromatic(m, index))")
end

for z in 1:8
  println("$(z) in $(z in m)")
end
