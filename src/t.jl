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

  export Molecule, SetOfAtoms, Atom, Bond, ChemicalStandardisation, BondList, Mol2Graph

  import Base: getindex, iterate, in, length
  # Now done in C++
  # getindex(m::Molecule, a::Int)=atom(m, a)
  # getindex(a::Atom, b::Int)=atom.item(b)
  iterate(m::Molecule, state=1) = (state >= natoms(m) ? nothing : (m[state], state + 1))
  iterate(a::Atom, state=1) = (state >= ncon(a) ? nothing : (a[state], state + 1))
  iterate(b::Bond, state=1) = (state == 1 ? (b.a1(), 2) : state == 2 ? (b.a2(), 2) : nothing)
  iterate(b::BondList, state=1) = (state >= size(b) ? nothing : (b[state], state + 1))
  iterate(r::Ring, state=1) = (state > length(r) ? nothing : (r[state], state + 1))
  iterate(s::SetOfAtoms, state=1) = (state >= s.size() ? nothing : (s[state], state + 1))
  in(z::Int, m::Molecule) = (natoms(m, z) > 0)
  in(atom::Int, a::Atom) = involves(a, atom)
  in(s::SetOfAtoms, a::Int) = contains(s, a)
  length(m::Molecule) = natoms(m)
  length(r::Ring) = atoms_in_ring(r)
  # length(r::Ring) = size(r)
  size(m::Molecule) = natoms(m)
  size(r::Ring) = (atoms_in_ring(r),)
  export getindex
  export iterate
  export length
  export in
  export atoms_in_ring
  export natoms, smiles, unique_smiles, nrings, build_from_smiles, is_aromatic, set_name, name
  export atomic_number, molecular_formula, nedges, is_ring_atom, fused_system_size, fused_system_identifier
  export rings_with_fused_system_identifier, in_same_ring, in_same_aromatic_ring, in_same_ring_system
  export ring_membership, rings_containing_both, is_part_of_fused_ring_system, ring, ring_containing_atom
  export ncon, nbonds, involves, other
  export atomic_number

  export activate_all, process
end
#getindex(m::LillyMol.Molecule, a::Int64)=LillyMol.atom(m, a)

using .LillyMol
using Test

# Call greet and show the result
@show LillyMol.greet()

w = LillyMol.World()
println(LillyMol.greet(w));
LillyMol.set(w, "hello")
println(LillyMol.greet(w))

function test_build_from_smiles()::Bool
  m = Molecule()
  build_from_smiles(m, "c1ccccc1OC foo")
end

function test_set_name()::Bool
  m = Molecule()
  build_from_smiles(m, "C")
  set_name(m, "foo")
  name(m) == "foo"
end

function test_natoms()::Bool
  m = Molecule()
  for i in 1:20
    smiles = "C" ^ i
    build_from_smiles(m, smiles) || return false
    natoms(m) == i || return false
  end
  true
end

function test_natoms_atomic_number()::Bool
  m = Molecule()
  for i in 1:10
    smiles = "CN" ^ i
    build_from_smiles(m, smiles) || return false
    natoms(m, 6) == i || return false
    natoms(m, 7) == i || return false
  end
  true
end

function test_natoms_element()::Bool
  m = Molecule()
  for i in 1:10
    smiles = "CN" ^ i
    build_from_smiles(m, smiles) || return false
    natoms(m, "C") == i || return false
    natoms(m, "N") == i || return false
  end
  true
end

function test_length_molecule()::Bool
  m = Molecule()
  for i = 1:10
    smiles = "CN"^i
    build_from_smiles(m, smiles) || return false
    length(m) == (2 * i) || return false
  end
  true
end

function test_atomic_number()::Bool
  atomic_numbers = [9, 6, 7, 8, 92]
  m = Molecule()
  build_from_smiles(m, "FCNO.[U]")
  for (ndx, atom) in enumerate(m)
    atomic_number(m, ndx) == atomic_numbers[ndx] || return false
    atomic_number(atom) == atomic_numbers[ndx] || return false
  end
  true
end

function test_nedges()::Bool
  mol_and_edges = Dict{String, Int}(""=>0, "[Fe]"=>0, "II"=>1, "III"=>2, "I1II1"=>3)
  for (smi,edges) in mol_and_edges
    m = Molecule()
    build_from_smiles(m, smi) || return false
    nedges(m) == edges || "$(edges) failed" == ""
  end
  true
end

function test_atomic_number_in():Bool
  m = Molecule()
  for i in 1:10
    smiles = "CNOC(F)(F)" ^ i
    build_from_smiles(m, smiles) || return false
    for z in [6, 7, 8, 9]
      z in m || return false
    end
  end
  true
end

function test_molecular_formula()::Bool
  smi_to_mf = Dict{String, String}("C"=>"CH4", "CO"=>"CH4O",
        "C.C"=>"CH4.CH4", "c1ccccc1"=>"C6H6")
  for (smi, mf) in smi_to_mf
    m = Molecule()
    build_from_smiles(m, smi) || return false
    # println("Expct $(mf) get $(molecular_formula(m))")
    mf == molecular_formula(m) || return false
  end
  true
end

function test_nrings()::Bool
  smiles_and_nrings = Dict{String, Int}("C"=>0, "CC"=>0, "CCC"=>0,
        "C1CC1"=>1, "C12CC1C2"=>2, "C1CC1C1CC1"=>2)
  m = Molecule()
  for (smi, nr) in smiles_and_nrings
    build_from_smiles(m, smi) || return false
    nrings(m) == nr || return false
  end
  true
end

function test_nrings_atom()::Bool
  smiles_and_nrings = Dict{String, Array{Int}}("C"=>[0], "CC"=>[0,0], "CCC"=>[0,0,0],
        "C1CC1"=>[1,1,1], "C12CC1C2"=>[2,1,2,1], "C1CC1C1CC1"=>[1,1,1,1,1,1])
  m = Molecule()
  for (smi, rings) in smiles_and_nrings
    build_from_smiles(m, smi) || return false
    for (atom, nr) in enumerate(rings)
      nrings(m, atom) == nr || return false
    end
  end
  true
end

function test_nrings_size()::Bool
  m = Molecule()
  for rsize in 3:10
    smiles = "C1" * "C"^(rsize-1) * "1"
    build_from_smiles(m, smiles) || return false
    for i in 1:natoms(m)
      nrings(m, i, rsize) == 1 || return false
    end
  end
  build_from_smiles(m, "C12CC1C2") || return false
  nrings(m, 1, 3) == 2 || return false
  nrings(m, 2, 3) == 1 || return false
  nrings(m, 3, 3) == 2 || return false
  nrings(m, 4, 3) == 1 || return false
  true
end

function test_is_ring_atom()::Bool
  m = Molecule()
  build_from_smiles(m, "C")
  is_ring_atom(m, 1) == false || return false
  build_from_smiles(m, "C1CC1C")
  is_ring_atom(m, 1) == true || return false
  is_ring_atom(m, 2) == true || return false
  is_ring_atom(m, 3) == true || return false
  is_ring_atom(m, 4) == false || return false

  true
end

function test_ring_bond_count()::Bool
  m = Molecule()
  build_from_smiles(m, "C")
  ring_bond_count(m, 1) == 0 || return false
  build_from_smiles(m, "C1CC1C")
  ring_bond_count(m, 1) == 1 || return false
  ring_bond_count(m, 2) == 1 || return false
  ring_bond_count(m, 3) == 1 || return false
  ring_bond_count(m, 4) == 0 || return false

  build_from_smiles(m, "C12CC2C1") || return false
  ring_bond_count(m, 1) == 2 || return false
  ring_bond_count(m, 2) == 1 || return false
  ring_bond_count(m, 3) == 2 || return false
  ring_bond_count(m, 4) == 1 || return false

  true
end

function test_fused_system_size()::Bool
  m = Molecule();
  build_from_smiles(m, "C") || return false
  fused_system_size(m, 1) == 0 || return false
  build_from_smiles(m, "C1CC1C") || return false
  fused_system_size(m, 1) == 1 || return false
  fused_system_size(m, 4) == 0 || return false

  build_from_smiles(m, "C12CC2C1") || return false
  for i in 1:natoms(m)
    fused_system_size(m, i) == 2 || return false
  end
  true
end

function test_fused_system_identifier()::Bool
  m = Molecule()
  build_from_smiles(m, "CC") || return false
  fused_system_identifier(m, 1) == fused_system_identifier(m, 2) || return false
  build_from_smiles(m, "C1CC1") || return false
  fused_system_identifier(m, 1) == fused_system_identifier(m, 2) || return false
  fused_system_identifier(m, 2) == fused_system_identifier(m, 3) || return false
  build_from_smiles(m, "C12CC2C1") || return false
  for i in 2:natoms(m)
    fused_system_identifier(m, 1) == fused_system_identifier(m, i) || return false
  end
  true
end

function test_rings_with_fused_system_identifier()::Bool
  m = Molecule()
  build_from_smiles(m, "C1CC1")
  fsid = fused_system_identifier(m, 1);
  rings_with_fused_system_identifier(m, fsid) == 1 || return false
  build_from_smiles(m, "C1CC1C1CC1")
  fsid = fused_system_identifier(m, 1);
  rings_with_fused_system_identifier(m, fsid) == 1 || return false
  fsid = fused_system_identifier(m, 4);
  rings_with_fused_system_identifier(m, fsid) == 1 || return false
  build_from_smiles(m, "C12CC2C1") || return false
  fsid = fused_system_identifier(m, 4);
  rings_with_fused_system_identifier(m, fsid) == 2 || return false
  true
end

function test_in_same_ring()::Bool
  m = Molecule();
  build_from_smiles(m, "CC") || return false
  in_same_ring(m, 1, 2) && return false
  build_from_smiles(m, "C1CC1C") || return false
  in_same_ring(m, 1, 2) || return false
  in_same_ring(m, 2, 3) || return false
  in_same_ring(m, 1, 3) || return false
  ! in_same_ring(m, 1, 4) || return false
  ! in_same_ring(m, 2, 4) || return false
  ! in_same_ring(m, 3, 4) || return false

  build_from_smiles(m, "C12CC2C1") || return false
  in_same_ring(m, 1, 2) || return false
  in_same_ring(m, 1, 3) || return false
  in_same_ring(m, 1, 4) || return false
  in_same_ring(m, 2, 3) || return false
  ! in_same_ring(m, 2, 4) || return false
  in_same_ring(m, 3, 4) || return false
end

function test_in_same_aromatic_ring()::Bool
  m = Molecule()
  build_from_smiles(m, "Oc1cccnc1") || return false
  for i in 2:7
    ! in_same_aromatic_ring(m, 1, i) || return false
  end
  for i in 2:7
    for j in i:7
      in_same_aromatic_ring(m, i, j) || return false
    end
  end
  build_from_smiles(m, "c12cncc2cccc1") || return false
  for i in 2:9
    in_same_aromatic_ring(m, 1, i) || return false
  end
  in_same_aromatic_ring(m, 2, 3) || return false
  in_same_aromatic_ring(m, 2, 4) || return false
  in_same_aromatic_ring(m, 2, 5) || return false
  ! in_same_aromatic_ring(m, 2, 6) || return false
  ! in_same_aromatic_ring(m, 6, 2) || return false
  ! in_same_aromatic_ring(m, 2, 7) || return false
  ! in_same_aromatic_ring(m, 7, 2) || return false
  ! in_same_aromatic_ring(m, 2, 8) || return false
  ! in_same_aromatic_ring(m, 2, 9) || return false
end

function test_in_same_ring_system()::Bool
  m = Molecule()
  build_from_smiles(m, "C1CC1CC1CC1") || return false
  in_same_ring_system(m, 1, 2) || return false
  in_same_ring_system(m, 1, 3) || return false
  in_same_ring_system(m, 2, 3) || return false
  in_same_ring_system(m, 5, 6) || return false
  in_same_ring_system(m, 5, 7) || return false
  in_same_ring_system(m, 6, 7) || return false
  ! in_same_ring_system(m, 3, 4) || return false
  ! in_same_ring_system(m, 4, 5) || return false
  ! in_same_ring_system(m, 3, 5) || return false
  ! in_same_ring_system(m, 2, 6) || return false
  ! in_same_ring_system(m, 1, 7) || return false
end

function test_ring_membership()::Bool
  m = Molecule();
  build_from_smiles(m, "CCC") || return false
  r = ring_membership(m)
  sum(r) == 0 || return false
  build_from_smiles(m, "C1CC1") || return false
  r = ring_membership(m)
  # https://stackoverflow.com/questions/47564825/check-if-all-the-elements-of-a-julia-array-are-equal#:~:text=The%20shortest%20way%20I%20can,%3D%3D%20arr)%20.
  # all(y->y==r[1], r) || return false
  all(y->y==1, r) || return false
  build_from_smiles(m, "C12CC2C1")
  r = ring_membership(m)
  [2, 1, 2, 1] == r || return false
  true
end

function test_rings_containing_both()::Bool
  m = Molecule();
  build_from_smiles(m, "C1CC1C") || return false
  rings_containing_both(m, 1, 2) == 1 || return false
  rings_containing_both(m, 1, 3) == 1 || return false
  rings_containing_both(m, 2, 3) == 1 || return false
  rings_containing_both(m, 1, 4) == 0 || return false
  rings_containing_both(m, 2, 4) == 0 || return false
  rings_containing_both(m, 3, 4) == 0 || return false
  build_from_smiles(m, "C12CC2C1C") || return false
  rings_containing_both(m, 1, 2) == 1 || return false
  rings_containing_both(m, 1, 3) == 2 || return false
  rings_containing_both(m, 1, 4) == 1 || return false
  rings_containing_both(m, 2, 3) == 1 || return false
  rings_containing_both(m, 3, 4) == 1 || return false
  rings_containing_both(m, 4, 5) == 0 || return false
  true
end

function test_is_part_of_fused_ring_system()::Bool
  m = Molecule()
  build_from_smiles(m, "CCC") || return false
  ! is_part_of_fused_ring_system(m, 1) || return false
  build_from_smiles(m, "C1CC1") || return false
  ! is_part_of_fused_ring_system(m, 1) || return false
  ! is_part_of_fused_ring_system(m, 2) || return false
  build_from_smiles(m, "C12CC2C1CC") || return false
  is_part_of_fused_ring_system(m, 1) || return false
  is_part_of_fused_ring_system(m, 2) || return false
  is_part_of_fused_ring_system(m, 3) || return false
  is_part_of_fused_ring_system(m, 4) || return false
  ! is_part_of_fused_ring_system(m, 5) || return false
end

function test_ring()::Bool
  m = Molecule();
  build_from_smiles(m, "C1CC1C") || return false
  r = ring(m, 1)
  l = length(r)
  length(r) == 3 || return false
  for i in 1:3
    i in r || return false
  end
  atoms = collect(r)
  for i in 1:3
    i in atoms || return false
  end
  4 in atoms && return false

  build_from_smiles(m, "C1CC1CC1CC1") || return false
  for i in 1:nrings(m)
    r = ring(m, i)
    length(r) == 3 || return false
  end
  all([x in ring(m, 1) for x in 1:3]) || return false
  all([x in ring(m, 2) for x in 5:7]) || return false
  true
end

# Broken test
function test_rings()::Bool
  m = Molecule()
  build_from_smiles("C1CC1C2CCC2C3CCCC3") || return false
  for (i, r) in enumerate(sssr_rings(m))
  end
  true
end

function test_ring_containing_atom()::Bool
  m = Molecule();
  build_from_smiles(m, "C1CC1C2CCC2C3CCCC3") || return false
  for i in 1:3
    r = ring_containing_atom(m, i)
    length(r) == 3 || return false
  end
  for i in 4:7
    r = ring_containing_atom(m, i)
    length(r) == 4 || return false
  end
  for i in 8:natoms(m)
    r = ring_containing_atom(m, i)
    length(r) == 5 || return false
  end
  true
end

function test_standardise()
  standardise = ChemicalStandardisation()
  activate_all(standardise)
  m = Molecule()
  build_from_smiles(m, "CC[N+](=O)[O-]")
  process(standardise, m)
  return smiles(m) == "CCN(=O)=O"
end

#for (index,atom) in enumerate(m)
#  print("atom $(index) type $(atomic_number(atom)) connected to")
#  for bond in atom
#    nbr = other(bond, index)
#    print(" $(nbr)")
#  end
#  println("")
#end

@test test_build_from_smiles()
@test test_set_name()
@test test_natoms()
@test test_natoms_atomic_number()
@test test_natoms_element()
@test test_length_molecule()
@test test_atomic_number()
@test test_atomic_number_in()
@test test_nedges()
@test test_molecular_formula()
@test test_standardise()
@test test_nrings()
@test test_nrings_atom()
@test test_nrings_size()
@test test_is_ring_atom()
@test test_fused_system_size()
@test test_fused_system_identifier()
@test test_rings_with_fused_system_identifier()
@test test_in_same_ring()
@test test_in_same_aromatic_ring()
@test test_in_same_ring_system()
@test test_ring_membership()
@test test_rings_containing_both()
@test test_is_part_of_fused_ring_system()
@test test_ring()
@test test_rings() broken=true
@test test_ring_containing_atom()
