module LillyMol
  using CxxWrap

  @wrapmodule(joinpath("bazel-bin/julia/","lillymol_julia.so"))

  function __init__()
    @initcxx
  end
  export World
  export greet

  export FileType
  export BondType, SINGLE, DOUBLE, TRIPLE ,AROMATIC

  export Molecule, SetOfAtoms, Atom, Bond, ChemicalStandardisation, BondList, Mol2Graph

  import Base: getindex, iterate, in, length
  # Now done in C++
  # getindex(m::Molecule, a::Int)=atom(m, a)
  # getindex(a::Atom, b::Int)=atom.item(b)
  iterate(m::Molecule, state=0) = (state >= natoms(m) ? nothing : (m[state], state + 1))
  iterate(a::Atom, state=0) = (state >= ncon(a) ? nothing : (a[state], state + 1))
  iterate(b::Bond, state=1) = (state == 1 ? (b.a1(), 2) : state == 2 ? (b.a2(), 2) : nothing)
  iterate(b::BondList, state=0) = (state >= size(b) ? nothing : (b[state], state + 1))
  iterate(r::Ring, state=0) = (state >= length(r) ? nothing : (r[state], state + 1))
  iterate(s::SetOfAtoms, state=0) = (state >= s.size() ? nothing : (s[state], state + 1))
  in(z::Int, m::Molecule) = (natoms(m, z) > 0)
  in(atom::Int, a::Atom) = involves(a, atom)
  length(m::Molecule) = natoms(m)
  length(r::Ring) = atoms_in_ring(r)
  # length(r::Ring) = size(r)
  size(m::Molecule) = natoms(m)
  size(r::Ring) = (atoms_in_ring(r),)
  size(s::SetOfAtoms) = (length(s),)
  size(s::SetOfAtomsAllocated) = (length(s),)
  export getindex
  export iterate
  export length
  export in
  export atoms_in_ring, contains
  export natoms, smiles, unique_smiles, nrings, build_from_smiles, is_aromatic, set_name, name
  export atomic_number, molecular_formula, nedges, is_ring_atom, fused_system_size, fused_system_identifier
  export rings_with_fused_system_identifier, in_same_ring, in_same_aromatic_ring, in_same_ring_system
  export ring_membership, rings_containing_both, is_part_of_fused_ring_system, ring, ring_containing_atom
  export label_atoms_by_ring_system, label_atoms_by_ring_system_including_spiro_fused
  export nrings_including_non_sssr_rings, non_sssr_rings, non_sssr_ring, is_spiro_fused, is_halogen
  export maximum_connectivity, connections, isotopically_labelled_smiles, is_aromatic, atom, formal_charge
  export isotope, number_formally_charged_atoms, net_formal_charge, bond, bond_between_atoms
  export compute_aromaticity_if_needed, bond_between_atoms, number_symmetry_classes, symmetry_class, symmetry_equivalents
  export symmetry_classes, attached_heteroatom_count, add_bond, are_bonded, bond_between_atoms
  export bond_length, bond_angle, dihedral_angle
  export ncon, nbonds, involves, other
  export is_single_bond, is_double_bond, is_triple_bond, is_aromatic
  export atomic_number
  export a1, a2

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
    atomic_number(m, ndx - 1) == atomic_numbers[ndx] || return false
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
      atom = atom - 1
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
      nrings(m, i - 1, rsize) == 1 || return false
    end
  end
  build_from_smiles(m, "C12CC1C2") || return false
  nrings(m, 0, 3) == 2 || return false
  nrings(m, 1, 3) == 1 || return false
  nrings(m, 2, 3) == 2 || return false
  nrings(m, 3, 3) == 1 || return false
  true
end

function test_is_ring_atom()::Bool
  m = Molecule()
  build_from_smiles(m, "C")
  is_ring_atom(m, 0) == false || return false
  build_from_smiles(m, "C1CC1C")
  is_ring_atom(m, 0) == true || return false
  is_ring_atom(m, 1) == true || return false
  is_ring_atom(m, 2) == true || return false
  is_ring_atom(m, 3) == false || return false

  true
end

function test_ring_bond_count()::Bool
  m = Molecule()
  build_from_smiles(m, "C")
  ring_bond_count(m, 1) == 0 || return false
  build_from_smiles(m, "C1CC1C")
  ring_bond_count(m, 0) == 1 || return false
  ring_bond_count(m, 1) == 1 || return false
  ring_bond_count(m, 2) == 1 || return false
  ring_bond_count(m, 3) == 0 || return false

  build_from_smiles(m, "C12CC2C1") || return false
  ring_bond_count(m, 0) == 2 || return false
  ring_bond_count(m, 1) == 1 || return false
  ring_bond_count(m, 2) == 2 || return false
  ring_bond_count(m, 3) == 1 || return false

  true
end

function test_fused_system_size()::Bool
  m = Molecule();
  build_from_smiles(m, "C") || return false
  fused_system_size(m, 0) == 0 || return false
  build_from_smiles(m, "C1CC1C") || return false
  fused_system_size(m, 0) == 1 || return false
  fused_system_size(m, 3) == 0 || return false

  build_from_smiles(m, "C12CC2C1") || return false
  for i in 1:natoms(m)
    fused_system_size(m, i - 1) == 2 || return false
  end
  true
end

function test_fused_system_identifier()::Bool
  m = Molecule()
  build_from_smiles(m, "CC") || return false
  fused_system_identifier(m, 0) == fused_system_identifier(m, 1) || return false
  build_from_smiles(m, "C1CC1") || return false
  fused_system_identifier(m, 0) == fused_system_identifier(m, 1) || return false
  fused_system_identifier(m, 1) == fused_system_identifier(m, 2) || return false
  build_from_smiles(m, "C12CC2C1") || return false
  for i in 2:natoms(m)
    fused_system_identifier(m, 0) == fused_system_identifier(m, i - 1) || return false
  end
  true
end

function test_rings_with_fused_system_identifier()::Bool
  m = Molecule()
  build_from_smiles(m, "C1CC1")
  fsid = fused_system_identifier(m, 0);
  rings_with_fused_system_identifier(m, fsid) == 1 || return false
  build_from_smiles(m, "C1CC1C1CC1")
  fsid = fused_system_identifier(m, 0);
  rings_with_fused_system_identifier(m, fsid) == 1 || return false
  fsid = fused_system_identifier(m, 3);
  rings_with_fused_system_identifier(m, fsid) == 1 || return false
  build_from_smiles(m, "C12CC2C1") || return false
  fsid = fused_system_identifier(m, 3);
  rings_with_fused_system_identifier(m, fsid) == 2 || return false
  true
end

function test_in_same_ring()::Bool
  m = Molecule();
  build_from_smiles(m, "CC") || return false
  in_same_ring(m, 0, 1) && return false
  build_from_smiles(m, "C1CC1C") || return false
  in_same_ring(m, 0, 1) || return false
  in_same_ring(m, 1, 2) || return false
  in_same_ring(m, 0, 2) || return false
  ! in_same_ring(m, 0, 3) || return false
  ! in_same_ring(m, 1, 3) || return false
  ! in_same_ring(m, 2, 3) || return false

  build_from_smiles(m, "C12CC2C1") || return false
  in_same_ring(m, 0, 1) || return false
  in_same_ring(m, 0, 2) || return false
  in_same_ring(m, 0, 3) || return false
  in_same_ring(m, 1, 2) || return false
  ! in_same_ring(m, 1, 3) || return false
  in_same_ring(m, 2, 3) || return false
end

function test_in_same_aromatic_ring()::Bool
  m = Molecule()
  build_from_smiles(m, "Oc1cccnc1") || return false
  for i in 1:6
    in_same_aromatic_ring(m, 0, i - 1) && return false
  end
  for i in 1:6
    for j in i:6
      in_same_aromatic_ring(m, i, j) || return false
    end
  end
  build_from_smiles(m, "c12cncc2cccc1") || return false
  in_same_aromatic_ring(m, 1, 2) || return false
  in_same_aromatic_ring(m, 1, 3) || return false
  in_same_aromatic_ring(m, 1, 4) || return false
  ! in_same_aromatic_ring(m, 1, 5) || return false
  ! in_same_aromatic_ring(m, 5, 1) || return false
  ! in_same_aromatic_ring(m, 1, 6) || return false
  ! in_same_aromatic_ring(m, 6, 1) || return false
  ! in_same_aromatic_ring(m, 1, 7) || return false
  ! in_same_aromatic_ring(m, 1, 8) || return false
end

function test_in_same_ring_system()::Bool
  m = Molecule()
  build_from_smiles(m, "C1CC1CC1CC1") || return false
  in_same_ring_system(m, 0, 1) || return false
  in_same_ring_system(m, 0, 2) || return false
  in_same_ring_system(m, 1, 2) || return false
  in_same_ring_system(m, 4, 5) || return false
  in_same_ring_system(m, 4, 6) || return false
  in_same_ring_system(m, 5, 6) || return false
  ! in_same_ring_system(m, 2, 3) || return false
  ! in_same_ring_system(m, 3, 4) || return false
  ! in_same_ring_system(m, 2, 4) || return false
  ! in_same_ring_system(m, 1, 5) || return false
  ! in_same_ring_system(m, 0, 6) || return false
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
  rings_containing_both(m, 0, 1) == 1 || return false
  rings_containing_both(m, 0, 2) == 1 || return false
  rings_containing_both(m, 1, 2) == 1 || return false
  rings_containing_both(m, 0, 3) == 0 || return false
  rings_containing_both(m, 1, 3) == 0 || return false
  rings_containing_both(m, 2, 3) == 0 || return false
  build_from_smiles(m, "C12CC2C1C") || return false
  rings_containing_both(m, 0, 1) == 1 || return false
  rings_containing_both(m, 0, 2) == 2 || return false
  rings_containing_both(m, 0, 3) == 1 || return false
  rings_containing_both(m, 1, 2) == 1 || return false
  rings_containing_both(m, 2, 3) == 1 || return false
  rings_containing_both(m, 3, 4) == 0 || return false
  true
end

function test_is_part_of_fused_ring_system()::Bool
  m = Molecule()
  build_from_smiles(m, "CCC") || return false
  ! is_part_of_fused_ring_system(m, 0) || return false
  build_from_smiles(m, "C1CC1") || return false
  ! is_part_of_fused_ring_system(m, 0) || return false
  ! is_part_of_fused_ring_system(m, 1) || return false
  build_from_smiles(m, "C12CC2C1CC") || return false
  is_part_of_fused_ring_system(m, 0) || return false
  is_part_of_fused_ring_system(m, 1) || return false
  is_part_of_fused_ring_system(m, 2) || return false
  is_part_of_fused_ring_system(m, 3) || return false
  ! is_part_of_fused_ring_system(m, 4) || return false
end

function test_ring()::Bool
  m = Molecule();
  build_from_smiles(m, "C1CC1C") || return false
  r = ring(m, 0)
  l = length(r)
  length(r) == 3 || return false
  for i in 0:2
    i in r || return false
  end
  atoms = collect(r)
  for i in 0:2
    i in atoms || return false
  end
  3 in atoms && return false

  build_from_smiles(m, "C1CC1CC1CC1") || return false
  for i in 1:nrings(m)
    r = ring(m, i - 1)
    length(r) == 3 || return false
  end
  all([x in ring(m, 0) for x in 0:2]) || return false
  all([x in ring(m, 1) for x in 4:6]) || return false
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
  for i in 0:2
    r = ring_containing_atom(m, i)
    length(r) == 3 || return false
  end
  for i in 3:6
    r = ring_containing_atom(m, i)
    length(r) == 4 || return false
  end
  for i in 7:(natoms(m) - 1)
    r = ring_containing_atom(m, i)
    length(r) == 5 || return false
  end
  true
end

# Fix for zero arrays
function test_label_atoms_by_ring_system()::Bool
  m = LillyMol.MolFromSmiles("C1CC1C2CCC2C3CCCCC3C")
  rsys = label_atoms_by_ring_system(m)
  for i in 1:3
    rsys[i] == 1 || return false
  end
  for i in 4:7
    rsys[i] == 2 || return false
  end
  for i in 8:13
    rsys[i] == 3 || return false
  end
  rsys[14] == 0 || return false

  true
end

# Fix for zero arrays
function test_label_atoms_by_ring_system_including_spiro_fused()::Bool
  m = LillyMol.MolFromSmiles("C1CC12CC2")
  rsys = label_atoms_by_ring_system_including_spiro_fused(m)
  for i in 1:natoms(m)
    rsys[i] == 1 || return false
  end
  true
end

function test_nrings_including_non_sssr_rings()::Bool
  m = LillyMol.MolFromSmiles("C12C3C4C1C5C2C3C45")
  for i in 0:(natoms(m) - 1)
    nrings_including_non_sssr_rings(m, i) == 3 || return false
  end
  true
end

function test_non_sssr_rings()::Bool
  m = LillyMol.MolFromSmiles("C12C3C4C1C5C2C3C45")
  non_sssr_rings(m) == 1 || return false
  true
end

function test_non_sssr_ring()::Bool
  m = LillyMol.MolFromSmiles("C12C3C4C1C5C2C3C45")
  # 1 == non_sssr_rings(m) || return false
  r = non_sssr_ring(m, 0)
  length(r) == 4 || return false
  true
end

function test_is_spiro_fused()::Bool
  m = LillyMol.MolFromSmiles("C1CC12CC2")
  is_spiro_fused(m, 0) && return false
  is_spiro_fused(m, 1) && return false
  is_spiro_fused(m, 2) || return false
  is_spiro_fused(m, 3) && return false
  is_spiro_fused(m, 4) && return false
  true
end

function test_is_halogen()::Bool
  m = LillyMol.MolFromSmiles("FC(Cl)(Br)CI")
  is_halogen(m, 0) || return false
  is_halogen(m, 1) && return false
  is_halogen(m, 2) || return false
  is_halogen(m, 3) || return false
  is_halogen(m, 4) && return false
  is_halogen(m, 5) || return false
  true
end

function test_ncon_molecule()::Bool
  m = LillyMol.MolFromSmiles("CCC(C)C(C)(C)C")
  ncon(m, 0) == 1 || return false
  ncon(m, 1) == 2 || return false
  ncon(m, 2) == 3 || return false
  ncon(m, 4) == 4 || return false
  true
end

function test_nbonds_molecule()::Bool
  m = LillyMol.MolFromSmiles("CCC=CC#C")
  nbonds(m, 0) == 1 || return false
  nbonds(m, 1) == 2 || return false
  nbonds(m, 2) == 3 || return false
  nbonds(m, 4) == 4 || return false
  true
end

function test_maximum_connectivity()::Bool
  m = Molecule()
  maximum_connectivity(m) == 0 || return false
  build_from_smiles(m, "C") || return false
  maximum_connectivity(m) == 0 || return false
  build_from_smiles(m, "CC") || return false
  maximum_connectivity(m) == 1 || return false
  build_from_smiles(m, "CCC") || return false
  maximum_connectivity(m) == 2 || return false
  build_from_smiles(m, "CC(C)C") || return false
  maximum_connectivity(m) == 3 || return false
  build_from_smiles(m, "CC(C)(C)C") || return false
  maximum_connectivity(m) == 4 || return false
  return true
end

function test_connections_molecule()::Bool
  m = Molecule();
  build_from_smiles(m, "C") || return false
  c = connections(m, 0)
  length(c) == 0 || return false
  build_from_smiles(m, "CC") || return false
  c = connections(m, 0)
  length(c) == 1 || return false
  c[1] == 1 || return false
  c = connections(m, 1)
  length(c) == 1 || return false
  c[1] == 0 || return false
  build_from_smiles(m, "CCC") || return false
  c = connections(m, 1)
  (0 in c && 2 in c) || return false
  build_from_smiles(m, "CC(C)(C)C") || return false
  c = connections(m, 1)
  length(c) == 4 || return false
  (0 in c && 2 in c && 3 in c && 4 in c) || return false
  true
end

function test_isotopically_labelled_smiles()::Bool
  m = Molecule()
  build_from_smiles(m, "CNO")
  s = isotopically_labelled_smiles(m)
  s == "C[1NH][2OH]" || return false
  true
end

function test_is_aromatic()::Bool
  m = LillyMol.MolFromSmiles("C1CCC1")
  for i in natoms(m)
    is_aromatic(m, i - 1) && return false
  end

  build_from_smiles(m, "c1ccccc1") || return false
  for i in natoms(m)
    is_aromatic(m, i - 1) || return false
  end
  true
end

function test_getindex_molecule()::Bool
  m = LillyMol.MolFromSmiles("[1CH3]NOF")
  atomic_number(m[0]) == 6 || return false
  atomic_number(m, 0) == 6 || return false
  atomic_number(m[1]) == 7 || return false
  atomic_number(m, 1) == 7 || return false
  atomic_number(m[2]) == 8 || return false
  atomic_number(m, 2) == 8 || return false
  atomic_number(m[3]) == 9 || return false
  atomic_number(m, 3) == 9 || return false
  ncon(m[0]) == 1 || return false
  ncon(m, 0) == 1 || return false
  ncon(m[1]) == 2 || return false
  ncon(m, 1) == 2 || return false
  ncon(m[2]) == 2 || return false
  ncon(m, 2) == 2 || return false
  ncon(m[3]) == 1 || return false
  ncon(m, 3) == 1 || return false
  formal_charge(m[0]) == 0 || return false
  formal_charge(m, 0) == 0 || return false
  isotope(m[0]) == 1 || return false
  isotope(m, 0) == 1 || return false
  isotope(m[1]) == 0 || return false
  isotope(m, 1) == 0 || return false
  true
end

function test_number_formally_charged_atoms()::Bool
  m = LillyMol.MolFromSmiles("NC")
  number_formally_charged_atoms(m) == 0 || return false
  build_from_smiles(m, "[NH3+]C")
  number_formally_charged_atoms(m) == 1 || return false
  build_from_smiles(m, "C(=O)[O-]")
  number_formally_charged_atoms(m) == 1 || return false
  true
end

function test_net_formal_charge()::Bool
  m = LillyMol.MolFromSmiles("NC")
  net_formal_charge(m) == 0 || return false
  build_from_smiles(m, "[NH3+]C")
  net_formal_charge(m) == 1 || return false
  build_from_smiles(m, "C(=O)[O-]")
  net_formal_charge(m) == -1 || return false
  true
end

function test_bond_molecule()::Bool
  m = LillyMol.MolFromSmiles("CC=CC#C")
  b = bond(m, 0)
  is_single_bond(b) || return false
  is_double_bond(bond(m, 1)) || return false
  is_triple_bond(bond(m, 3)) || return false
  build_from_smiles(m, "c1ccccc1")
  compute_aromaticity_if_needed(m)
  for i in 1:nedges(m)
    is_aromatic(bond(m, i - 1)) || return false
  end
  true
end

function test_bond_between_atoms()::Bool
  m = LillyMol.MolFromSmiles("CC=C")
  b = bond_between_atoms(m, 0, 1)
  is_single_bond(b) || return false
  b = bond_between_atoms(m, 1, 2)
  is_double_bond(b) || return false
  true
end

function test_number_symmetry_classes()::Bool
  m = LillyMol.MolFromSmiles("C")
  number_symmetry_classes(m) == 1 || return false
  build_from_smiles(m, "CC") || return false
  number_symmetry_classes(m) == 1 || return false
  build_from_smiles(m, "FC(F)(F)C(C)C") || return false
  number_symmetry_classes(m) == 4 || return false
  build_from_smiles(m, "C1CC1") || return false
  number_symmetry_classes(m) == 1 || return false
  true
end

function test_symmetry_class()::Bool
  m = Molecule()
  build_from_smiles(m, "FC(F)(F)C(C)C") || return false
  symmetry_class(m, 0) == symmetry_class(m, 2) || return false
  symmetry_class(m, 0) == symmetry_class(m, 3) || return false
  symmetry_class(m, 5) == symmetry_class(m, 6) || return false
  build_from_smiles(m, "CN")
  symmetry_class(m, 0) == symmetry_class(m, 1) && return false
  true
end

# in SetOfAtomsAllocated is not working
function test_symmetry_equivalents()::Bool
  m = Molecule()
  build_from_smiles(m, "FC(F)(F)C(C)C") || return false
  s = symmetry_equivalents(m, 0)
  length(s) == 2 || return false
  (2 in s && 3 in s) || return false
  s = symmetry_equivalents(m, 6)
  length(s) == 1 || return false
  (6 in s) || return false
  true
end

# Fix for zero arrays
function test_symmetry_classes()::Bool
  m = Molecule()
  build_from_smiles(m, "FC(F)(F)C(C)C") || return false
  c = symmetry_classes(m)
  c[1] == c[3] || return false
  c[1] == c[4] || return false
  c[6] == c[7] || return false
  true
end

function test_attached_heteroatom_count()::Bool
  m = Molecule()
  build_from_smiles(m, "FC(F)(F)C(C)N") || return false
  attached_heteroatom_count(m, 0) == 0 || return false
  attached_heteroatom_count(m, 1) == 3 || return false
  attached_heteroatom_count(m, 2) == 0 || return false
  attached_heteroatom_count(m, 4) == 1 || return false
  attached_heteroatom_count(m, 5) == 0 || return false
  attached_heteroatom_count(m, 6) == 0 || return false
  true
end

function test_bond_length()::Bool
  m = Molecule()
  build_from_smiles(m, "C{{0,0,0}}C{{1,1,1}}") || return false
  @test bond_length(m, 0, 1) ≈ sqrt(3.0) atol=0.001
  true
end

function test_bond_angle()::Bool
  m = Molecule()
  build_from_smiles(m, "C{{0,0,0}}C{{1,0,0}}C{{1,1,1}}") || return false
  @test bond_angle(m, 0, 1, 2) ≈ (π / 2.0) atol=0.001
  true
end

function test_dihedral_angle()::Bool
  m = Molecule()
  build_from_smiles(m, "C{{0,0,0}}C{{1,0,0}}C{{1,1,1}}") || return false
  @test bond_angle(m, 0, 1, 2) ≈ (π / 2.0) atol=0.001
  true
end

function test_add_bond()::Bool
  m = Molecule()
  build_from_smiles(m, "C.C") || return false
  add_bond(m, 0, 1, SINGLE)
  smiles(m) == "CC" || return false
  build_from_smiles(m, "C.C") || return false
  add_bond(m, 0, 1, DOUBLE)
  smiles(m) == "C=C" || return false
  build_from_smiles(m, "C.C") || return false
  add_bond(m, 0, 1, TRIPLE)
  smiles(m) == "C#C" || return false
  true
end

function test_are_bonded()::Bool
  m = Molecule()
  build_from_smiles(m, "C.CC")
  are_bonded(m, 0, 1) && return false
  are_bonded(m, 1, 2) || return false
  true
end

function test_bond_between_atoms()::Bool
  m = Molecule()
  build_from_smiles(m, "CCC") || return false
  b = bond_between_atoms(m, 0, 1)
  ((a1(b) == 0 && a2(b) == 1) || (a1(b) == 1 && a2(b) == 0)) || return false
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
@test test_label_atoms_by_ring_system()
@test test_label_atoms_by_ring_system_including_spiro_fused()
@test test_nrings_including_non_sssr_rings()
@test test_non_sssr_rings()
@test test_non_sssr_ring()
@test test_is_spiro_fused()
@test test_is_halogen()
@test test_ncon_molecule()
@test test_nbonds_molecule()
@test test_maximum_connectivity()
@test test_connections_molecule()
@test test_isotopically_labelled_smiles()
@test test_is_aromatic()
@test test_getindex_molecule()
@test test_number_formally_charged_atoms()
@test test_net_formal_charge()
@test test_bond_molecule()
@test test_bond_between_atoms()
@test test_number_symmetry_classes()
@test test_symmetry_class()
@test test_symmetry_equivalents() broken=true
@test test_symmetry_classes()
@test test_attached_heteroatom_count()
@test test_bond_length()
@test test_bond_angle()
#@test test_dihedral_angle()
#@test test_signed_dihedral_angle()
#@test test_set_bond_length()
#@test test_set_bond_angle()
#@test test_set_dihedral_angle()
#@test test_xyx_molecule()
@test test_add_bond()
@test test_are_bonded()
@test test_bond_between_atoms()
