module LillyMol
  using CxxWrap

  @wrapmodule(joinpath("bazel-bin/julia/","lillymol_julia.so"))

  function __init__()
    @initcxx
  end
  export World
  export greet

  export FileType
  export BondType

  export Molecule, SetOfAtoms, Atom, Bond, ChemicalStandardisation

  import Base: getindex, iterate, in
  # Now done in C++
  # getindex(m::Molecule, a::Int)=atom(m, a)
  # getindex(a::Atom, b::Int)=atom.item(b)
  iterate(m::Molecule, state=1) = (println("foo state $(state) $(atomic_number(m[state]))"); state >= natoms(m) ? nothing : (m[state], state + 1))
  iterate(a::Atom, state=1) = (state >= ncon(a) ? nothing : (a[state], state + 1))
  in(z::Int, m::Molecule) = (natoms(m, z) > 0)
  in(atom::Int, a::Atom) = involves(a, atom)
  in(s::SetOfAtoms, a::Int) = contains(s, a)
  length(m::Molecule) = natoms(m)
  size(m::Molecule) = natoms(m)
  export getindex
  export iterate
  export length
  export in
  export natoms, smiles, unique_smiles, nrings, build_from_smiles, is_aromatic
  export atomic_number
  export ncon, nbonds, involves, other
  export atomic_number

  export activate_all, process
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

standardise = ChemicalStandardisation()
activate_all(standardise)
build_from_smiles(m, "CC[N+](=O)[O-]")
process(standardise, m)
println("After std $(smiles(m))")

for (index,atom) in enumerate(m)
  print("atom $(index) type $(atomic_number(atom)) connected to")
  for bond in atom
    nbr = other(bond, index)
    print(" $(nbr)")
  end
  println("")
end
