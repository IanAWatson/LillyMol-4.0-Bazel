#!/usr/bin/env ruby

# frozen_string_literal: true

# convert a SubstructureQuery::Query proto to C++ code

require_relative 'lib/iwcmdline'

require 'Molecule_Lib/substructure_pb'

require 'protobuf'

def add_property(qry, ndx, prop, fname, src)
  pmin = eval("qry.min_#{prop}")
  pmax = eval("qry.max_#{prop}")
  pvalues = eval("qry.#{prop}")
  return if pmin == 0 && pmax == 0 && pvalues.empty?
  src << "int"
  src << "Ok#{prop.capitalize}#{ndx}(Molecule_to_Match& target) {"
  src << "  const auto #{prop} = target.#{fname}();"
  if pmin > 0
    src << "  if (#{prop} < #{pmin}) {"
    src << "    return 0;"
    src << "  }"
  end
  if pmax > 0
    src << "  if (#{prop} > #{pmax}) {"
    src << "    return 0;"
    src << "  }"
  end
  if ! pvalues.empty?
    src << "  static const resizable_array<int> okvalues{#{pvalues.join(',')}};"
    src << "  if (! okvalues.contains(#{prop})) {"
    src << "    return 0;"
    src << "  }"
  end
  src << "  return 1;"
  src << "}"
  src
end

def add_atoms_in_spinach(qry, ndx, src)
  return if qry.min_atoms_in_spinach == 0 && qry.max_atoms_in_spinach == 0 && qry.atoms_in_spinach.empty?
  src << "int\nOkAtomsInSpinach#{ndx}(Molecule_to_Match& target) {";
    src << "  Molecule& m = *target.molecule();"
    src << "  const int matoms = m.natoms();"
    src << "  std::unique_ptr<int[]> spinach(new_int(matoms));"
    src << "  const int atoms_in_spinach = m.identify_spinach(spinach.get());"
    if qry.min_atoms_in_spinach > 0
      src << "  if (atoms_in_spinach < #{qry.min_atoms_in_spinach}) {"
      src << "    return 0;";
      src << "  }"
    end
    if qry.max_atoms_in_spinach > 0
      src << "  if (atoms_in_spinach > #{qry.max_atoms_in_spinach})"
      src << "     return 0;";
      src << "  }"
    end
    if ! qry.atoms_in_spinach.empty?
      src << "  static const resizable_array<int> ok_values{#{qry.atoms_in_spinach.join(',') }};"
      src << "  if (! ok_values.contains(atoms_in_spinach) {"
      src << "    return 0;"
      src << "  }"
    end
  src << "  return 1;"
  src << "}"
end

def add_inter_ring_atoms(qry, ndx, src)
  return if qry.min_inter_ring_atoms == 0 && qry.max_inter_ring_atoms == 0 && qry.inter_ring_atoms.empty?
  src << "int\nOkAtomsInSpinach#{ndx}(Molecule_to_Match& target) {";
    src << "  Molecule& m = *target.molecule();"
    src << "  const int matoms = m.natoms();"
    src << "  std::unique_ptr<int[]> spinach(new_int(matoms));"
    src << "  const int atoms_in_spinach = m.identify_spinach(spinach.get());"
    src << "  int inter_ring_atoms = 0;"
    src << "  for (int i = 0; i < matoms; ++i) {"
    src << "    if (spinach[i] == 0 && target[i].is_non_ring_atom()) {"
    src << "      ++inter_ring_atoms;"
    src << "    }"
    src << "  }"
    if qry.min_inter_ring_atoms > 0
      src << "  if (inter_ring_atoms < #{qry.min_inter_ring_atoms}) {"
      src << "    return 0;";
      src << "  }"
    end
    if qry.max_inter_ring_atoms > 0
      src << "  if (inter_ring_atoms > #{qry.max_inter_ring_atoms})"
      src << "     return 0;";
      src << "  }"
    end
    if qry.inter_ring_atoms.empty?
      src << "  return 1;"
    end

    if qry.inter_ring_atoms.size == 1
      src << "  return inter_ring_atoms == #{qry.inter_ring_atoms[0]};"
    end

    src << "  static const resizable_array<int> ok_values{#{qry.inter_ring_atoms.join(',') }};"
    src << "  return ok_values.contains(inter_ring_atoms);"
  src << "}"
end

def add_any_net_formal_charge(qry, ndx, src)
  return unless qry.net_formal_charge
  # Note that this is not implemented.
  src << "int\nOkAnyNetFormalCharge#{ndx}(Molecule_to_Match& target) {"
  src << "  // not implemented, todo"
  src << "  return 1;"
  src << "  // return target.net_formal_charge() != 0;"
  src << "}"
  src

end

def add_aromatic_atoms(qry, ndx, src)
  return unless qry.min_aromatic_atoms > 0 && qry.max_aromatic_atoms > 0 && qry.aromatic_atoms.empty?
  src << "int\nOkAromaticAtoms#{ndx}(Molecule_to_Match& target){"
  src << "  const int count = target.molecule()->aromatic_atom_count();"
  if qry.min_aromatic_atoms > 0
    src << "  if (count < #{qry.min_aromatic_atoms}) {"
    src << "    return 0;"
    src << "  }"
  end
  if qry.max_aromatic_atoms > 0
    src << "  if (count > #{qry.max_aromatic_atoms}) {"
    src << "    return 0;"
    src << "  }"
  end
  src << "  return 1;"
  src << "}"
  src
end

def add_elements_needed(qry, ndx, src)
  return if qry.elements_needed.empty?
  qry.elements_needed.each_with_index do |e, i|
    src << "int\nOkElementsNeeded#{ndx}_#{i}(Molecule_to_Match& target) {"
    if e.atomic_number > 0
      src << "  const Element* e = get_element_from_atomic_number(#{e.atomic_number});"
    else
      src << "  const Element* e = get_element_from_symbol_no_case_conversion(#{e.atomic_symbol});"
    end
    src << "  const int count = target.molecule()->natoms(e);";
    if e.min_hits_needed > 0
      src << "  if (count < #{e.min_hits_needed}) {"
      src << "    return 0;"
      src << "  }"
    end
    if e.max_hits_needed
      src << "  if (count > #{e.max_hits_needed}) {"
      src << "    return 0;"
      src << "  }"
    end
    unless e.hits_needed.empty?
      src << "  static const resizable_array<int> okvalues{#{e.hits_needed.join(',')}};"
      src << "  if (! okvalues.contains(count)) {"
      src << "    return 0;"
      src << "  }"
    end
    src << "  return 1;"
    src << "}"
  end
  src
end

def add_ring_specifiers(qry, ndx, src)
  return if qry.ring_specifier.empty?
  qry.ring_specifier.each_with_index do |r, i|
    $stderr << "What is aromatic #{r.aromatic} has #{r.has_aromatic?}\n"
    $stderr << " in #{r}\n"
    src << "int\nOkRingSpecifier#{ndx}_#{i}(Molecule_to_Match& target) {";
    src << "  Molecule& m = *target.molecule();"
    if ! r.base.ncon.empty?
      src << "  static const resizable_array<int> ok_ncon{#{r.base.ncon.join(',')}};"
    end
    if ! r.ring_size.empty?
      src << "  static const resizable_array<int> okringsize{#{r.ring_size.join(',')}};"
    end
    if ! r.base.hits_needed.empty?
      src << "  static const resizable_array<int> oknhits{#{r.base.hits_needed.join(',')}};"
    end
    if ! r.fused.empty?
      src << "  static const resizable_array<int> okfused{#{r.fused.join(',')}};"
    end
    if ! r.fused_aromatic_neighbours.empty?
      src << "  static const resizable_array<int> okfused_arom{#{r.fused_aromatic_neighbours.join(',')}};"
    end
    if ! r.fused_non_aromatic_neighbours.empty?
      src << "  static const resizable_array<int> okfused_non_arom{#{r.fused_non_aromatic_neighbours.join(',')}};"
    end
    if ! r.base.attached_heteroatom_count.empty?
      src << "  static const resizable_array<int> ok_attached_heteroatom_count{#{r.base.attached_heteroatom_count.join(',')}};"
    end
    if ! r.base.within_ring_unsaturation.empty?
      src << "  static const resizable_array<int> ok_within_ring_unsaturation{#{r.base.within_ring_unsaturation.join(',')}};"
    end
    if ! r.base.largest_number_of_bonds_shared_with_another_ring.empty?
      src << "  static const resizable_array<int> ok_largest_number_of_bonds_shared_with_another_ring{#{r.base.largest_number_of_bonds_shared_with_another_ring.join(',')}};"
    end
    if ! r.base.atoms_with_pi_electrons.empty?
      src << "  static const resizable_array<int> ok_atoms_with_pi_electrons{#{r.base.atoms_with_pi_electrons.join(',')}};"
    end
    if ! r.base.strongly_fused_ring_neighbours.empty?
      src << "  static const resizable_array<int> ok_strongly_fused_ring_neighbours{#{r.base.strongly_fused_ring_neighbours.join(',')}};"
    end
    if r.base.all_hits_in_same_fragment == true
      src << "  extending_resizable_array<int> hits_in_fragment;"
    end
    src << "  int nhits = 0;"
    src << "  for (const Ring* r : target.molecule()->sssr_rings()) {"
    if r.min_ring_size > 0
      src << "    if (r-size() < #{r.min_ring_size}) {"
      src << "      continue;"
      src << "    }"
    end
    if r.max_ring_size > 0
      src << "    if (r-size() > #{r.max_ring_size}) {"
      src << "      continue;"
      src << "    }"
    end
    if ! r.ring_size.empty?
      src << "    if (! okringsize.contains(r->size())) {"
      src << "      continue;"
      src << "    }"
    end
    if ! r.has_aromatic?
    elsif r.aromatic == true
      src << "    if (! r->is_aromatic()) {"
      src << "      continue;"
      src << "    }"
    elsif r.aromatic == false
      src << "    if (r->is_aromatic()) {"
      src << "      continue;"
      src << "    }"
    end
    if r.min_fused > 0 || r.max_fused > 0 || ! r.fused.empty?
      src << "    const int fused = r->fused_ring_neighbours();"
    end
    if r.min_fused > 0
      src << "    if (fused < #{r.min_fused}) {"
      src << "      continue;"
      src << "    }"
    end
    if r.max_fused > 0
      src << "    if (fused > #{r.max_fused}) {"
      src << "      continue;"
      src << "    }"
    end
    if ! r.fused.empty?
      src << "    if (! okfused.contains(fused)) {"
      src << "      continue;"
      src << "    }"
    end
    if r.min_fused_aromatic_neighbours > 0 || r.max_fused_aromatic_neighbours > 0 || ! r.fused_aromatic_neighbours.empty?
      src << "    const int fused_aromatic_neighbours = ss_ring::fused_aromatic_neighbours(*r);"
    end
    if r.min_fused_aromatic_neighbours > 0
      src << "    if (fused_aromatic_neighbours < #{r.min_fused_aromatic_neighbours}) {"
      src << "      continue;"
      src << "    }"
    end
    if r.max_fused_aromatic_neighbours > 0
      src << "    if (fused_aromatic_neighbours > #{r.max_fused_aromatic_neighbours}) {"
      src << "      continue;"
      src << "    }"
    end
    if ! r.fused_aromatic_neighbours.empty?
      src << "    if (! okfused_arom.contains(fused_aromatic_neighbours)) {"
      src << "      continue;"
      src << "    }"
    end
    if r.min_fused_non_aromatic_neighbours > 0 || r.max_fused_non_aromatic_neighbours > 0 || ! r.fused_non_aromatic_neighbours.empty?
      src << "    const int fused_non_aromatic_neighbours = ss_ring::fused_non_aromatic_neighbours(*r);"
    end
    if r.min_fused_non_aromatic_neighbours > 0
      src << "    if (fused_non_aromatic_neighbours < #{r.min_fused_non_aromatic_neighbours}) {"
      src << "      continue;"
      src << "    }"
    end
    if r.max_fused_non_aromatic_neighbours > 0
      src << "    if (fused_non_aromatic_neighbours > #{r.max_fused_non_aromatic_neighbours}) {"
      src << "      continue;"
      src << "    }"
    end
    if ! r.fused_non_aromatic_neighbours.empty?
      src << "    if (! okfused_non_arom.contains(fused_non_aromatic_neighbours)) {"
      src << "      continue;"
      src << "    }"
    end

    need_attached_heteroatom_count = false
    if r.base.min_attached_heteroatom_count > 0 || r.base.max_attached_heteroatom_count > 0 || ! r.base.attached_heteroatom_count.empty?
      src << "    int attached_heteroatom_count = 0;"
      need_attached_heteroatom_count = true
    end
    need_ncon = false
    if r.base.min_ncon > 0 || r.base.max_ncon > 0 || ! r.base.ncon.empty?
      src << "    int ncon = 0;"
      need_ncon = true
    end
    need_heteroatom_count = false
    if r.base.min_heteroatom_count > 0 || r.base.max_heteroatom_count > 0 || ! r.base.heteroatom_count.empty?
      src << "    int heteroatom_count = 0;"
      need_heteroatom_count = true
    end
    if need_ncon || need_heteroatom_count || need_attached_heteroatom_count
      src << "    const Set_of_Atoms& s = *r;"
      src << "    for (auto a : s) {"
      if need_ncon
        src << "      if (int x = m.ncon(a); x > 2) { ncon += x - 2;}"
      end
      if need_heteroatom_count
        src << "      if (m.atomic_number(a) != 6) { ++attached_heteroatom_count;}"
      end
      if need_attached_heteroatom_count
        src << "      const Atom& atom = m.atom(a);"
        src << "      for (const Bond* b : atom) {"
        src << "        atom_number_t j = b->other(a);"
        src << "        if (r->contains(j)) continue;"
        src << "        if (m.atomic_number(j) != 6) {++attached_heteroatom_count;}"
        src << "      }"
      end
      src << "    }"
      if r.base.min_attached_heteroatom_count > 0
        src << "    if (attached_heteroatom_count < #{r.base.min_attached_heteroatom_count}) {"
        src << "      continue;"
        src << "    }"
      end
      if r.base.max_attached_heteroatom_count > 0
        src << "    if (attached_heteroatom_count > #{r.base.max_attached_heteroatom_count}) {"
        src << "      continue;"
        src << "    }"
      end
      if ! r.base.attached_heteroatom_count.empty?
        src << "    if (! ok_attached_heteroatom_count.contains(attached_heteroatom_count)) {"
        src << "      continue;"
        src << "    }"
      end
      if r.base.min_heteroatom_count > 0
        src << "    if (heteroatom_count < #{r.base.min_heteroatom_count}) {"
        src << "      continue;"
        src << "    }"
      end
      if r.base.max_heteroatom_count > 0
        src << "    if (heteroatom_count > #{r.base.max_heteroatom_count}) {"
        src << "      continue;"
        src << "    }"
      end
      if ! r.base.heteroatom_count.empty?
        src << "    if (! ok_heteroatom_count.contains(heteroatom_count)) }"
        src << "      continue;"
        src << "    }"
      end
      if r.base.min_ncon > 0
        src << "    if (ncon < #{r.base.min_ncon}) {"
        src << "      continue;"
        src << "    }"
      end
      if r.base.max_ncon > 0
        src << "    if (ncon > #{r.base.max_ncon}) {"
        src << "      continue;"
        src << "    }"
      end
      if ! r.base.ncon.empty?
        src << "    if (! ok_ncon.contains(ncon)) {"
        src << "      continue;"
        src << "    }"
      end
    end
    if r.base.min_within_ring_unsaturation > 0 || r.base.max_within_ring_unsaturation > 0 || ! r.base.within_ring_unsaturation.empty?
      src << "    const int within_ring_unsaturation = ss_ring::compute_within_ring_unsaturation(*r, target);"
    end
    if r.base.min_within_ring_unsaturation > 0
      src << "    if (within_ring_unsaturation < #{r.base.min_within_ring_unsaturation}) {"
      src << "      continue;"
      src << "    }"
    end
    if r.base.max_within_ring_unsaturation > 0
      src << "    if (within_ring_unsaturation > #{r.base.max_within_ring_unsaturation}) {"
      src << "      continue;"
      src << "    }"
    end
    if ! r.base.within_ring_unsaturation.empty?
      src << "    if (! ok_within_ring_unsaturation.contains(within_ring_unsaturation)) {"
      src << "      continue;"
      src << "    }"
    end
    if r.base.min_largest_number_of_bonds_shared_with_another_ring > 0 || r.base.max_largest_number_of_bonds_shared_with_another_ring > 0 || ! r.base.largest_number_of_bonds_shared_with_another_ring.empty?
      src << "    const int largest_number_of_bonds_shared_with_another_ring = r->largest_number_of_bonds_shared_with_another_ring();"
    end
    if r.base.min_largest_number_of_bonds_shared_with_another_ring > 0
      src << "    if (largest_number_of_bonds_shared_with_another_ring < #{r.base.min_largest_number_of_bonds_shared_with_another_ring}) {"
      src << "      continue;"
      src << "    }"
    end
    if r.base.max_largest_number_of_bonds_shared_with_another_ring > 0
      src << "    if (largest_number_of_bonds_shared_with_another_ring > #{r.base.max_largest_number_of_bonds_shared_with_another_ring}) {"
      src << "      continue;"
      src << "    }"
    end
    if ! r.base.largest_number_of_bonds_shared_with_another_ring.empty?
      src << "    if (! ok_largest_number_of_bonds_shared_with_another_ring.contains(largest_number_of_bonds_shared_with_another_ring)) {"
      src << "      continue;"
      src << "    }"
    end
    if r.base.min_atoms_with_pi_electrons > 0 || r.base.max_atoms_with_pi_electrons > 0 || ! r.base.atoms_with_pi_electrons.empty?
      src << "    const int atoms_with_pi_electrons = ss_ring::compute_atoms_with_pi_electrons(r, target);"
    end
    if r.base.min_atoms_with_pi_electrons > 0
      src << "    if (atoms_with_pi_electrons < #{r.base.min_atoms_with_pi_electrons}) {"
      src << "      continue;"
      src << "    }"
    end
    if r.base.max_atoms_with_pi_electrons > 0
      src << "    if (atoms_with_pi_electrons > #{r.base.max_atoms_with_pi_electrons}) {"
      src << "      continue;"
      src << "    }"
    end
    if ! r.base.atoms_with_pi_electrons.empty?
      src << "    if (! ok_atoms_with_pi_electrons.contains(atoms_with_pi_electrons)) {"
      src << "      continue;"
      src << "    }"
    end
    if r.base.min_strongly_fused_ring_neighbours > 0 || r.base.max_strongly_fused_ring_neighbours > 0 || ! r.base.strongly_fused_ring_neighbours.empty?
      src << "    const int strongly_fused_ring_neighbours = r->strongly_fused_ring_neighbours();"
    end
    if r.base.min_strongly_fused_ring_neighbours > 0
      src << "    if (strongly_fused_ring_neighbours < #{r.base.min_strongly_fused_ring_neighbours}) {"
      src << "      continue;"
      src << "    }"
    end
    if r.base.max_strongly_fused_ring_neighbours > 0
      src << "    if (strongly_fused_ring_neighbours > #{r.base.max_strongly_fused_ring_neighbours}) {"
      src << "      continue;"
      src << "    }"
    end
    if ! r.base.strongly_fused_ring_neighbours.empty?
      src << "    if (! ok_strongly_fused_ring_neighbours.contains(strongly_fused_ring_neighbours)) {"
      src << "      continue;"
      src << "    }"
    end

    if ! r.base.environment.empty?
      $stderr << "Ring environments not supported #{r.base.environment}, ignored\n"
    end

    if r.base.all_hits_in_same_fragment == true
      src << "    ++hits_in_fragment[r->item(0)];"
    end
    src << "    ++nhits;"
    src << "  }"

    src << "  if (nhits == 0) {"
    if r.base.match_as_match == true
      src << "    return 0;"
    else
      src << "    return 1;"
    end
    src << "  }"

    if r.base.all_hits_in_same_fragment
      # todo implement this
    end

    if r.base.min_hits_needed > 0
      src << "  if (nhits < #{r.base.min_hits_needed}) {"
      src << "    return 0;"
      src << "  }"
    end
    if r.base.max_hits_needed > 0
      src << "    if (nhits > #{r.base.max_hits_needed}) {"
      src << "      return 0;"
      src << "    }"
    end
    if ! r.base.hits_needed.empty?
      src << "  if (! oknhits.contains(nhits)) {"
      src << "    return 0;"
      src << "  }"
    end
    src << "  return 1;"
    src << "}"
  end
  src
end

def add_ring_system_specifiers(qry, ndx, src)
  return if qry.ring_system_specifier.empty?
  $stderr << "Ring system specifiers not implemente\n"
end

# Add a check for a possibly multi-valued value.
# The expression to compute the value is `var`.
# The allowed values are in `values`, which might contain 0... entries.
# If it is empty, we do nothing.
# If it contains 1 value, we do a scalar check against that value.
# If it contains 2 values, we do a scalar check against both values.
# If it contains more values, we instantiate a static resizable_array
# and use the contains method. In that case, we use `vname` to name
# the variable.
def add_check_values(vname, var, values, src)
  return if values.empty?
  if values.size == 1
    src << "  if (#{var} != #{values[0]}) return 0;"
    return
  end
  if values.size == 2
    src << "  if (#{var} != #{values[0]} && #{var} != #{values[1]}) return 0;"
  end

  ok_values = "ok_#{vname}"
  src << "  static const resizable_array<int> #{ok_values}{#{values.join(',')}};"
  src << "  if (! #{ok_values}.contains(#{var}) return 0;"
end

# `qry` is a SubstructureAtomSpecifier, a matching function for it.
def add_specifier_matcher(qry, uid, src)
  src << "int\nMatchAtomSpec#{uid}(Target_Atom& atom) {"
  add_check_values('atomic_number', 'atom.atomic_number()', qry.atomic_number, src)

  src << "  if (atom.ncon() < #{qry.min_ncon}) return 0;" if qry.has_min_ncon?
  src << "  if (atom.ncon() > #{qry.max_ncon}) return 0;" if qry.has_max_ncon?
  add_check_values('ncon', 'atom.ncon()', qry.ncon, src)

  src << "  if (atom.nbonds() < #{qry.min_nbonds}) return 0;" if qry.has_min_nbonds?
  src << "  if (atom.nbonds() > #{qry.max_nbonds}) return 0;" if qry.has_max_nbonds?
  add_check_values('nbonds', 'atom.nbonds()', qry.nbonds, src)

  src << "  if (atom.formal_charge() < #{qry.min_formal_charge}) return 0;" if qry.has_min_formal_charge? 
  src << "  if (atom.formal_charge() > #{qry.max_formal_charge}) return 0;" if qry.has_max_formal_charge?
  add_check_values('formal_charge', 'atom.formal_charge()', qry.formal_charge, src)

  src << "  if (atom.nrings() < #{qry.min_nrings}) return 0;" if qry.has_min_nrings? 
  src << "  if (atom.nrings() > #{qry.max_nrings}) return 0;" if qry.has_max_nrings?
  add_check_values('nrings', 'atom.nrings()', qry.nrings, src)

  src << "  if (atom.ring_bond_count() < #{qry.min_ring_bond_count}) return 0;" if qry.has_min_ring_bond_count? 
  src << "  if (atom.ring_bond_count() > #{qry.max_ring_bond_count}) return 0;" if qry.has_max_ring_bond_count?
  add_check_values('ring_bond_count', 'atom.ring_bond_count()', qry.ring_bond_count, src)

  src << "  if (atom.ring_size() < #{qry.min_ring_size}) return 0;" if qry.has_min_ring_size? 
  src << "  if (atom.ring_size() > #{qry.max_ring_size}) return 0;" if qry.has_max_ring_size?
  add_check_values('ring_size', 'atom.ring_size()', qry.ring_size, src)

  src << "  if (atom.hcount() < #{qry.min_hcount}) return 0;" if qry.has_min_hcount? 
  src << "  if (atom.hcount() > #{qry.max_hcount}) return 0;" if qry.has_max_hcount?
  add_check_values('hcount', 'atom.hcount()', qry.hcount, src)

  if ! qry.has_aromatic?
  elsif qry.aromatic
    src << "  if (! atom.aromatic()) return 0;"
  else
    src << "  if (atom.aromatic()) return 0;"
  end

  if ! qry.has_chirality?
  elsif qry.chirality
    src << "  if (! atom.chiral_centre()) return 0;"
  else
    src << "  if (atom.chiral_centre()) return 0;"
  end

  src << "  if (atom.aromatic_ring_size() < #{qry.min_aromatic_ring_size}) return 0;" if qry.has_min_aromatic_ring_size? 
  src << "  if (atom.aromatic_ring_size() > #{qry.max_aromatic_ring_size}) return 0;" if qry.has_max_aromatic_ring_size?
  add_check_values('aromatic_ring_size', 'atom.aromatic_ring_size()', qry.aromatic_ring_size, src)

  src << "  if (atom.aliphatic_ring_size() < #{qry.min_aliphatic_ring_size}) return 0;" if qry.has_min_aliphatic_ring_size? 
  src << "  if (atom.aliphatic_ring_size() > #{qry.max_aliphatic_ring_size}) return 0;" if qry.has_max_aliphatic_ring_size?
  add_check_values('aliphatic_ring_size', 'atom.aliphatic_ring_size()', qry.aliphatic_ring_size, src)

  src << "  if (atom.attached_heteroatom_count() < #{qry.min_attached_heteroatom_count}) return 0;" if qry.has_min_attached_heteroatom_count? 
  src << "  if (atom.attached_heteroatom_count() > #{qry.max_attached_heteroatom_count}) return 0;" if qry.has_max_attached_heteroatom_count?
  add_check_values('attached_heteroatom_count', 'atom.attached_heteroatom_count()', qry.attached_heteroatom_count, src)

  src << "  if (atom.lone_pair_count() < #{qry.min_lone_pair_count}) return 0;" if qry.has_min_lone_pair_count? 
  src << "  if (atom.lone_pair_count() > #{qry.max_lone_pair_count}) return 0;" if qry.has_max_lone_pair_count?
  add_check_values('lone_pair_count', 'atom.lone_pair_count()', qry.lone_pair_count, src)

  src << "  if (atom.unsaturation() < #{qry.min_unsaturation}) return 0;" if qry.has_min_unsaturation? 
  src << "  if (atom.unsaturation.() > #{qry.max_unsaturation}) return 0;" if qry.has_max_unsaturation?
  add_check_values('unsaturation', 'atom.unsaturation()', qry.unsaturation, src)

  src << "  if (atom.daylight_x() < #{qry.min_daylight_x}) return 0;" if qry.has_min_daylight_x? 
  src << "  if (atom.daylight_x() > #{qry.max_daylight_x}) return 0;" if qry.has_max_daylight_x?
  add_check_values('daylight_x', 'atom.daylight_x()', qry.daylight_x, src)

  src << "  if (atom.isotope() < #{qry.min_isotope}) return 0;" if qry.has_min_isotope? 
  src << "  if (atom.isotope() > #{qry.max_isotope}) return 0;" if qry.has_max_isotope?
  add_check_values('isotope', 'atom.isotope()', qry.isotope, src)

  src << "  if (atom.aryl() < #{qry.min_aryl}) return 0;" if qry.has_min_aryl? 
  src << "  if (atom.aryl() > #{qry.max_aryl}) return 0;" if qry.has_max_aryl?
  add_check_values('aryl', 'atom.aryl()', qry.aryl, src)

  src << "  if (atom.vinyl() < #{qry.min_vinyl}) return 0;" if qry.has_min_vinyl? 
  src << "  if (atom.vinyl() > #{qry.max_vinyl}) return 0;" if qry.has_max_vinyl?
  add_check_values('vinyl', 'atom.vinyl()', qry.vinyl, src)

  src << "  if (atom.fused_system_size() < #{qry.min_fused_system_size}) return 0;" if qry.has_min_fused_system_size? 
  src << "  if (atom.fused_system_size() > #{qry.max_fused_system_size}) return 0;" if qry.has_max_fused_system_size?
  add_check_values('fused_system_size', 'atom.fused_system_size()', qry.fused_system_size, src)

  # all_rings_kekule????

  src << "  if (atom.heteroatoms_in_ring() < #{qry.min_heteroatoms_in_ring}) return 0;" if qry.has_min_heteroatoms_in_ring? 
  src << "  if (atom.heteroatoms_in_ring() > #{qry.max_heteroatoms_in_ring}) return 0;" if qry.has_max_heteroatoms_in_ring?
  add_check_values('heteroatoms_in_ring', 'atom.heteroatoms_in_ring()', qry.heteroatoms_in_ring, src)

  # match spinach only????

  src << "  if (atom.scaffold_bonds_attached_to_ring() < #{qry.min_scaffold_bonds_attached_to_ring}) return 0;" if qry.has_min_scaffold_bonds_attached_to_ring? 
  src << "  if (atom.scaffold_bonds_attached_to_ring() > #{qry.max_scaffold_bonds_attached_to_ring}) return 0;" if qry.has_max_scaffold_bonds_attached_to_ring?
  add_check_values('scaffold_bonds_attached_to_ring', 'atom.scaffold_bonds_attached_to_ring()', qry.scaffold_bonds_attached_to_ring, src)

  src << "  if (atom.symmetry_degree() < #{qry.min_symmetry_degree}) return 0;" if qry.has_min_symmetry_degree? 
  src << "  if (atom.symmetry_degree() > #{qry.max_symmetry_degree}) return 0;" if qry.has_max_symmetry_degree?
    add_check_values('symmetry_degree', 'atom.symmetry_degree()', qry.symmetry_degree, src)

  # symmetry group?
  # atom_type, user_atom_type????

  src << "  if (atom.valence.() < #{qry.min_valence}) return 0;" if qry.has_min_valence? 
  src << "  if (atom.valence.() > #{qry.max_valence}) return 0;" if qry.has_max_valence?
  add_check_values('valence', 'atom.valence()', qry.valence, src)

  src << "  return 1;"
  src << "}"
end

# Qry is a SubstructureAtom. Build a function that evaluates it.
def add_query_atom(qry, ndx1, ndx2, src)
  # if there is just one Specifier, things are easy
  if qry.atom_properties.size == 1
    src << "int\nMatchAtom#{ndx1}_#{ndx2}(Molecule_to_Match& target) {"
    src << "}"
    return
  end

  # If there are more than one Atom specifiers present, we must build the logical expression.
  $stderr << "SubstructureAtom has #{qry.atom_properties.size} specifiers\n"
  initialise_logexp_func = ""
  if qry.atom_properties.size > 1
  initialise_logexp_func = "LigExpAtomProp#{ndx1}_#{ndx2}"
    src << "IW_Logical_Expression\n#{initialise_logexp_func}() {"
    src << "  IW_Logical_Expression result;"
      qry.atom_properties.each_with_index do |p, i|
        next if i == 0
        $stderr << "Examining #{p} #{p.logical_operator.class}\n"
        if p.logical_operator == :SS_OR
          src << "  result.add_operator(IW_LOGEXP_OR);"
        elsif p.logical_operator == :SS_AND
          src << "  result.add_operator(IW_LOGEXP_AND);"
        elsif p.logical_operator == :SS_XOR
          src << "  result.add_operator(IW_LOGEXP_XOR);"
        elsif p.logical_operator == :SS_LP_AND
          src << "  result.add_operator(IW_LOGEXP_LOW_PRIORITY_AND);"
        else
          $stderr << "What is this operator #{p.logical_operator}\n"
        end
      end
    src << "  return result;"
    src << "}"
  end
  # First write a function for each of the specifiers.
   specifier_fns = []
   qry.atom_properties.each_with_index do |q, i|
     uid = "#{ndx1}_#{ndx2}_#{i}"
     specifier_fns << "MatchAtomSpec#{uid}"
     add_specifier_matcher(q, uid, src)
   end
   # Now write the Substructure_Atom matcher
   # Handle the case of just one atom_properties separately
   if qry.atom_properties.size == 1
     src << "  const int res = #{specifier_fns[0]}(target);"
     if qry.match_as_match
       src << "  return res;"
     else
       src << "  return !res;"
     end
     src << "}"
     return
   end
   # argh, need an array of functions in order to do this as a loop.
   # maybe should not use a c++ loop...

   src << "int\nMatchAtom#{ndx1}_#{ndx2}(Molecule_to_Match& target, int zatom) {"
   src << "  static IW_Logical_Expression logexp = #{initialise_logexp_func}();"
   src << "  int outcome = 0;"
   specifier_fns.each_with_index do |fn, i|
     src << "  if (logexp.result_needed(#{i})) {"
     src << "    int res = #{fn}(target[zatom]);"
     src << "    logexp.set_result(#{i}, res);"
     src << "    if (logexp.evaluate(outcome)) {"
     if ! qry.has_match_as_match? || qry.match_as_match
       src << "      return outcome;"
     else
       src << "      return !outcome;"
     end
     src << "    }"
     src << "  }"
   end
   src << "  // control should not come to here."
   src << "  return 0;"
   src << "}"
end

def add_query_atoms(qry, ndx, src)
  qry.query_atom.each_with_index do |q, i|
    add_query_atom(q, ndx, i, src)
  end
end

def generate_code(qry, ndx)
  result = Array.new

  add_property(qry, ndx, 'natoms', 'natoms', result)
  add_property(qry, ndx, 'nrings', 'nrings', result)
  add_property(qry, ndx, 'heteroatoms_in_molecule', 'heteroatom_count', result)
  add_property(qry, ndx, 'fused_rings', 'fused_rings', result)
  add_property(qry, ndx, 'strongly_fused_rings', 'strongly_fused_ring_count', result)
  add_property(qry, ndx, 'isolated_rings', 'isolated_ring_count', result)
  add_property(qry, ndx, 'isolated_ring_objects', 'ring_object_count', result)
  add_property(qry, ndx, 'aromatic_rings', 'aromatic_ring_count', result)
  add_property(qry, ndx, 'non_aromatic_rings', 'non_aromatic_ring_count', result)
  add_property(qry, ndx, 'number_isotopic_atoms', 'number_isotopic_atoms', result)
  add_property(qry, ndx, 'number_fragments', 'number_fragments', result)
  add_property(qry, ndx, 'net_formal_charge', 'net_formal_charge', result)
  add_any_net_formal_charge(qry, ndx, result)
  add_atoms_in_spinach(qry, ndx, result)
  add_inter_ring_atoms(qry, ndx, result)
  add_elements_needed(qry, ndx, result)
  add_ring_specifiers(qry, ndx, result)
  add_ring_system_specifiers(qry, ndx, result)
  add_query_atoms(qry, ndx, result)

  result
end

def write_header(output)
  header = %Q(
// Generated by #{$0}
#include <iostream>

#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/path.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

namespace qry2cpp {
)
  output << header
end

def write_closing(output)
  output << "}  // namespace qry2cpp\n"
end

def main
  cmdline = IWCmdline.new("-v-output=s")
  if ARGV.empty?
    $stderr << "Must specify proto file as argument\n"
  exit(1)
end

  binary_data = File.read(ARGV[0])

  qry = SubstructureSearch::SubstructureQuery.decode(binary_data)
  #$stderr << "Read #{qry}\n"

  per_query = []
  qry.query.each_with_index do |q, ndx|
    per_query << generate_code(q, ndx)
  end
  write_header($stdout)
  per_query.each do |q|
    $stdout << q.join("\n") << "\n"
  end
  write_closing($stdout)
end

main
