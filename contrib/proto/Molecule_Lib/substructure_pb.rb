# Generated by the protocol buffer compiler.  DO NOT EDIT!
# source: Molecule_Lib/substructure.proto

require 'google/protobuf'

require 'Molecule_Lib/geometric_constraints_pb'

Google::Protobuf::DescriptorPool.generated_pool.build do
  add_file("Molecule_Lib/substructure.proto", :syntax => :proto2) do
    add_message "SubstructureSearch.AtomNumberHydrogenLonePair" do
      oneof :AtomNumberOr do
        optional :atom_number, :uint32, 1
        optional :h_or_lp, :enum, 2, "SubstructureSearch.AtomNumberHydrogenLonePair.HorLP"
      end
    end
    add_enum "SubstructureSearch.AtomNumberHydrogenLonePair.HorLP" do
      value :HYDROGEN, 1
      value :LONEPAIR, 2
    end
    add_message "SubstructureSearch.SubstructureChiralCenter" do
      optional :center, :uint32, 1
      optional :top_front, :message, 2, "SubstructureSearch.AtomNumberHydrogenLonePair"
      optional :top_back, :message, 3, "SubstructureSearch.AtomNumberHydrogenLonePair"
      optional :left_down, :message, 4, "SubstructureSearch.AtomNumberHydrogenLonePair"
      optional :right_down, :message, 5, "SubstructureSearch.AtomNumberHydrogenLonePair"
    end
    add_message "SubstructureSearch.SubstructureBond" do
      repeated :bond_type, :enum, 1, "SubstructureSearch.BondType"
      optional :other_end, :uint32, 2
    end
    add_message "SubstructureSearch.SubstructureEnvironmentBond" do
      repeated :bond_type, :enum, 1, "SubstructureSearch.BondType"
      repeated :other_end, :uint32, 2
    end
    add_message "SubstructureSearch.ElementsNeeded" do
      repeated :hits_needed, :uint32, 3
      optional :min_hits_needed, :uint32, 4
      optional :max_hits_needed, :uint32, 5
      oneof :ElementSpecifier do
        optional :atomic_number, :uint32, 1
        optional :atomic_symbol, :string, 2
      end
    end
    add_message "SubstructureSearch.NoMatchedAtomsBetween" do
      optional :a1, :uint32, 1
      optional :a2, :uint32, 2
      optional :qualifier, :string, 3
    end
    add_message "SubstructureSearch.LinkAtoms" do
      optional :a1, :uint32, 1
      optional :a2, :uint32, 2
      repeated :distance, :uint32, 3
      optional :min_distance, :uint32, 4
      optional :max_distance, :uint32, 5
    end
    add_message "SubstructureSearch.EnvironmentAttachment" do
      repeated :attachment_point, :uint32, 1
      repeated :btype, :enum, 2, "SubstructureSearch.BondType"
      optional :substructure_bond, :string, 3
    end
    add_message "SubstructureSearch.SubstructureRingEnvironment" do
      optional :substructure_atom, :message, 1, "SubstructureSearch.SubstructureAtom"
      optional :min_hits_needed, :uint32, 2
      optional :max_hits_needed, :uint32, 3
    end
    add_message "SubstructureSearch.SubstructureRingBase" do
      optional :match_as_match, :bool, 1
      repeated :hits_needed, :uint32, 2
      optional :min_hits_needed, :uint32, 3
      optional :max_hits_needed, :uint32, 4
      repeated :attached_heteroatom_count, :uint32, 5
      optional :min_attached_heteroatom_count, :uint32, 6
      optional :max_attached_heteroatom_count, :uint32, 7
      repeated :heteroatom_count, :uint32, 8
      optional :min_heteroatom_count, :uint32, 9
      optional :max_heteroatom_count, :uint32, 10
      repeated :ncon, :uint32, 11
      optional :min_ncon, :uint32, 12
      optional :max_ncon, :uint32, 13
      optional :all_hits_in_same_fragment, :bool, 14
      repeated :within_ring_unsaturation, :uint32, 16
      optional :min_within_ring_unsaturation, :uint32, 17
      optional :max_within_ring_unsaturation, :uint32, 18
      repeated :largest_number_of_bonds_shared_with_another_ring, :uint32, 19
      optional :min_largest_number_of_bonds_shared_with_another_ring, :uint32, 20
      optional :max_largest_number_of_bonds_shared_with_another_ring, :uint32, 21
      repeated :atoms_with_pi_electrons, :uint32, 26
      optional :min_atoms_with_pi_electrons, :uint32, 27
      optional :max_atoms_with_pi_electrons, :uint32, 28
      repeated :strongly_fused_ring_neighbours, :uint32, 29
      optional :min_strongly_fused_ring_neighbours, :uint32, 30
      optional :max_strongly_fused_ring_neighbours, :uint32, 31
      optional :environment, :string, 22
      optional :environment_can_match_in_ring_atoms, :bool, 23
      optional :set_global_id, :uint32, 24
    end
    add_message "SubstructureSearch.SubstructureRingSpecification" do
      optional :base, :message, 1, "SubstructureSearch.SubstructureRingBase"
      repeated :ring_size, :uint32, 2
      optional :min_ring_size, :uint32, 3
      optional :max_ring_size, :uint32, 4
      optional :aromatic, :bool, 5
      repeated :fused, :uint32, 6
      optional :min_fused, :uint32, 7
      optional :max_fused, :uint32, 8
      repeated :fused_aromatic_neighbours, :uint32, 9
      optional :min_fused_aromatic_neighbours, :uint32, 10
      optional :max_fused_aromatic_neighbours, :uint32, 11
      repeated :fused_non_aromatic_neighbours, :uint32, 12
      optional :min_fused_non_aromatic_neighbours, :uint32, 13
      optional :max_fused_non_aromatic_neighbours, :uint32, 14
    end
    add_message "SubstructureSearch.RingSizeRequirement" do
      optional :ring_size, :uint32, 1
      repeated :count, :uint32, 2
      optional :min_count, :uint32, 3
      optional :max_count, :uint32, 4
    end
    add_message "SubstructureSearch.SubstructureRingSystemSpecification" do
      optional :base, :message, 1, "SubstructureSearch.SubstructureRingBase"
      repeated :rings_in_system, :uint32, 2
      optional :min_rings_in_system, :uint32, 3
      optional :max_rings_in_system, :uint32, 4
      repeated :ring_sizes, :uint32, 5
      optional :min_ring_sizes, :uint32, 6
      optional :max_ring_sizes, :uint32, 7
      repeated :ring_size_count, :message, 8, "SubstructureSearch.RingSizeRequirement"
      repeated :aromatic_ring_count, :uint32, 11
      optional :min_aromatic_ring_count, :uint32, 12
      optional :max_aromatic_ring_count, :uint32, 13
      repeated :non_aromatic_ring_count, :uint32, 14
      optional :min_non_aromatic_ring_count, :uint32, 15
      optional :max_non_aromatic_ring_count, :uint32, 16
      repeated :degree_of_fusion, :uint32, 17
      optional :min_degree_of_fusion, :uint32, 18
      optional :max_degree_of_fusion, :uint32, 19
      repeated :atoms_in_system, :uint32, 20
      optional :min_atoms_in_system, :uint32, 21
      optional :max_atoms_in_system, :uint32, 22
      repeated :number_spinach_groups, :uint32, 23
      optional :min_number_spinach_groups, :uint32, 24
      optional :max_number_spinach_groups, :uint32, 25
      repeated :number_non_spinach_groups, :uint32, 26
      optional :min_number_non_spinach_groups, :uint32, 27
      optional :max_number_non_spinach_groups, :uint32, 28
      repeated :atoms_in_spinach_group, :uint32, 29
      optional :min_atoms_in_spinach_group, :uint32, 30
      optional :max_atoms_in_spinach_group, :uint32, 31
      repeated :length_of_spinach_group, :uint32, 32
      optional :min_length_of_spinach_group, :uint32, 33
      optional :max_length_of_spinach_group, :uint32, 34
      repeated :distance_to_another_ring, :uint32, 35
      optional :min_distance_to_another_ring, :uint32, 36
      optional :max_distance_to_another_ring, :uint32, 37
      repeated :strongly_fused_ring_count, :uint32, 38
      optional :min_strongly_fused_ring_count, :uint32, 39
      optional :max_strongly_fused_ring_count, :uint32, 40
    end
    add_message "SubstructureSearch.SubstructureAtomSpecifier" do
      repeated :atomic_symbol, :string, 1
      repeated :atomic_number, :uint32, 2
      repeated :ncon, :uint32, 3
      optional :min_ncon, :uint32, 4
      optional :max_ncon, :uint32, 5
      repeated :ncon2, :uint32, 6
      optional :min_ncon2, :uint32, 7
      optional :max_ncon2, :uint32, 8
      repeated :nbonds, :uint32, 9
      optional :min_nbonds, :uint32, 10
      optional :max_nbonds, :uint32, 11
      repeated :formal_charge, :int32, 12
      optional :min_formal_charge, :int32, 13
      optional :max_formal_charge, :int32, 14
      repeated :nrings, :uint32, 15
      optional :min_nrings, :uint32, 16
      optional :max_nrings, :uint32, 17
      repeated :ring_bond_count, :uint32, 18
      optional :min_ring_bond_count, :uint32, 19
      optional :max_ring_bond_count, :uint32, 20
      repeated :ring_size, :uint32, 21
      optional :min_ring_size, :uint32, 22
      optional :max_ring_size, :uint32, 23
      repeated :hcount, :uint32, 24
      optional :min_hcount, :uint32, 25
      optional :max_hcount, :uint32, 26
      optional :aromatic, :bool, 27
      optional :chirality, :bool, 28
      repeated :aromatic_ring_size, :uint32, 30
      optional :min_aromatic_ring_size, :uint32, 31
      optional :max_aromatic_ring_size, :uint32, 32
      repeated :aliphatic_ring_size, :uint32, 33
      optional :min_aliphatic_ring_size, :uint32, 34
      optional :max_aliphatic_ring_size, :uint32, 35
      repeated :attached_heteroatom_count, :uint32, 36
      optional :min_attached_heteroatom_count, :uint32, 37
      optional :max_attached_heteroatom_count, :uint32, 38
      repeated :lone_pair_count, :uint32, 39
      optional :min_lone_pair_count, :uint32, 40
      optional :max_lone_pair_count, :uint32, 41
      repeated :unsaturation, :uint32, 42
      optional :min_unsaturation, :uint32, 43
      optional :max_unsaturation, :uint32, 44
      repeated :daylight_x, :uint32, 45
      optional :min_daylight_x, :uint32, 46
      optional :max_daylight_x, :uint32, 47
      repeated :isotope, :uint32, 48
      optional :min_isotope, :uint32, 49
      optional :max_isotope, :uint32, 50
      repeated :aryl, :uint32, 51
      optional :min_aryl, :uint32, 52
      optional :max_aryl, :uint32, 53
      repeated :vinyl, :uint32, 54
      optional :min_vinyl, :uint32, 55
      optional :max_vinyl, :uint32, 56
      repeated :fused_system_size, :uint32, 57
      optional :min_fused_system_size, :uint32, 58
      optional :max_fused_system_size, :uint32, 59
      optional :all_rings_kekule, :bool, 60
      repeated :heteroatoms_in_ring, :uint32, 61
      optional :min_heteroatoms_in_ring, :uint32, 62
      optional :max_heteroatoms_in_ring, :uint32, 63
      optional :match_spinach_only, :int32, 64
      repeated :scaffold_bonds_attached_to_ring, :uint32, 65
      optional :min_scaffold_bonds_attached_to_ring, :uint32, 66
      optional :max_scaffold_bonds_attached_to_ring, :uint32, 67
      optional :preference_value, :int32, 68
      repeated :symmetry_degree, :uint32, 69
      optional :min_symmetry_degree, :uint32, 70
      optional :max_symmetry_degree, :uint32, 71
      optional :symmetry_group, :int32, 72
      optional :logical_operator, :enum, 76, "SubstructureSearch.Operator"
      optional :user_atom_type, :uint32, 77
      optional :atom_type, :uint32, 78
      repeated :valence, :uint32, 79
      optional :min_valence, :uint32, 80
      optional :max_valence, :uint32, 81
    end
    add_message "SubstructureSearch.SubstructureAtomEnvironment" do
      optional :id, :uint32, 1
      repeated :substructure_atom, :message, 2, "SubstructureSearch.SubstructureAtom"
      optional :op, :string, 3
    end
    add_message "SubstructureSearch.SubstructureAtom" do
      optional :id, :int32, 1
      optional :match_as_match, :bool, 2
      optional :text_identifier, :string, 3
      optional :atom_map_number, :uint32, 4
      optional :initial_atom_number, :uint32, 5
      repeated :atom_properties, :message, 7, "SubstructureSearch.SubstructureAtomSpecifier"
      optional :ring_id, :int32, 9
      optional :fused_system_id, :uint32, 10
      optional :fragment_id, :int32, 11
      optional :numeric_value, :double, 12
      optional :include_in_embedding, :bool, 13
      repeated :environment, :message, 17, "SubstructureSearch.SubstructureAtomEnvironment"
      repeated :query_bond, :message, 21, "SubstructureSearch.SubstructureBond"
      optional :bond_smarts, :string, 22
      repeated :single_bond, :uint32, 25
      repeated :double_bond, :uint32, 26
      repeated :triple_bond, :uint32, 27
      repeated :aromatic_bond, :uint32, 28
      repeated :bond, :uint32, 29
      repeated :preference, :message, 23, "SubstructureSearch.SubstructureAtomSpecifier"
      optional :sum_all_preference_hits, :bool, 24
      repeated :unmatched_atoms_attached, :uint32, 30
      optional :min_unmatched_atoms_attached, :uint32, 31
      optional :max_unmatched_atoms_attached, :uint32, 32
      optional :atom_type_group, :uint32, 33
      optional :global_match_id, :uint32, 34
      repeated :children, :message, 35, "SubstructureSearch.SubstructureAtom"
      oneof :SmilesOrSmarts do
        optional :smarts, :string, 14
        optional :atom_smarts, :string, 15
        optional :smiles, :string, 16
      end
    end
    add_message "SubstructureSearch.SubstructureEnvironment" do
      optional :id, :uint32, 1
      repeated :smarts, :string, 3
      repeated :smiles, :string, 4
      repeated :query_atom, :message, 5, "SubstructureSearch.SubstructureAtom"
      optional :attachment, :message, 6, "SubstructureSearch.EnvironmentAttachment"
      repeated :bond, :string, 7
      optional :or_id, :uint32, 8
      optional :and_id, :uint32, 9
      repeated :hits_needed, :uint32, 10
      optional :min_hits_needed, :uint32, 11
      optional :max_hits_needed, :uint32, 12
      optional :no_other_substituents_allowed, :bool, 13
      optional :env_matches_can_share_attachment_points, :bool, 15
      optional :max_matches_to_find, :uint32, 16
      optional :hydrogen_ok, :bool, 17
      optional :max_env_matches_per_anchor, :uint32, 18
    end
    add_message "SubstructureSearch.MatchedAtomMatch" do
      repeated :atom, :int32, 1
      repeated :smarts, :string, 2
      optional :logexp, :enum, 3, "SubstructureSearch.Operator"
    end
    add_message "SubstructureSearch.SeparatedAtoms" do
      optional :a1, :uint32, 1
      optional :a2, :uint32, 2
      repeated :bonds_between, :uint32, 3
      optional :min_bonds_between, :uint32, 4
      optional :max_bonds_between, :uint32, 5
    end
    add_message "SubstructureSearch.SingleSubstructureQuery" do
      optional :id, :int32, 1
      optional :label, :string, 2
      optional :comment, :string, 3
      optional :one_embedding_per_start_atom, :bool, 4
      optional :normalise_rc_per_hits_needed, :uint32, 5
      optional :subtract_from_rc, :uint32, 6
      optional :max_matches_to_find, :uint32, 8
      optional :save_matched_atoms, :bool, 9
      optional :ncon_ignore_singly_connected, :bool, 10
      optional :perceive_symmetric_equivalents, :bool, 11
      optional :implicit_ring_condition, :uint32, 12
      optional :all_hits_in_same_fragment, :bool, 13
      optional :only_match_largest_fragment, :bool, 14
      optional :embeddings_do_not_overlap, :bool, 15
      optional :sort_by_preference_value, :bool, 16
      repeated :numeric_value, :double, 19
      repeated :no_matched_atoms_between, :message, 20, "SubstructureSearch.NoMatchedAtomsBetween"
      optional :no_matched_atoms_between_exhaustive, :bool, 21
      repeated :link_atoms, :message, 22, "SubstructureSearch.LinkAtoms"
      optional :fail_if_embeddings_too_close, :bool, 23
      optional :distance_between_hits_ncheck, :uint32, 24
      optional :sort_matches, :string, 25
      repeated :attached_heteroatom_count, :uint32, 26
      optional :min_attached_heteroatom_count, :uint32, 27
      optional :max_attached_heteroatom_count, :uint32, 28
      repeated :hits_needed, :uint32, 29
      optional :min_hits_needed, :uint32, 30
      optional :max_hits_needed, :uint32, 31
      repeated :ring_atoms_matched, :uint32, 32
      optional :min_ring_atoms_matched, :uint32, 33
      optional :max_ring_atoms_matched, :uint32, 34
      repeated :heteroatoms_matched, :uint32, 35
      optional :min_heteroatoms_matched, :uint32, 36
      optional :max_heteroatoms_matched, :uint32, 37
      repeated :heteroatoms_in_molecule, :uint32, 38
      optional :min_heteroatoms_in_molecule, :uint32, 39
      optional :max_heteroatoms_in_molecule, :uint32, 40
      repeated :natoms, :uint32, 41
      optional :min_natoms, :uint32, 42
      optional :max_natoms, :uint32, 43
      repeated :nrings, :uint32, 44
      optional :min_nrings, :uint32, 45
      optional :max_nrings, :uint32, 46
      repeated :ncon, :uint32, 47
      optional :min_ncon, :uint32, 48
      optional :max_ncon, :uint32, 49
      repeated :fused_rings, :uint32, 50
      optional :min_fused_rings, :uint32, 51
      optional :max_fused_rings, :uint32, 52
      repeated :strongly_fused_rings, :uint32, 53
      optional :min_strongly_fused_rings, :uint32, 54
      optional :max_strongly_fused_rings, :uint32, 55
      repeated :isolated_rings, :uint32, 56
      optional :min_isolated_rings, :uint32, 57
      optional :max_isolated_rings, :uint32, 58
      repeated :isolated_ring_objects, :uint32, 59
      optional :min_isolated_ring_objects, :uint32, 60
      optional :max_isolated_ring_objects, :uint32, 61
      repeated :aromatic_rings, :uint32, 62
      optional :min_aromatic_rings, :uint32, 63
      optional :max_aromatic_rings, :uint32, 64
      repeated :non_aromatic_rings, :uint32, 65
      optional :min_non_aromatic_rings, :uint32, 66
      optional :max_non_aromatic_rings, :uint32, 67
      repeated :distance_between_hits, :uint32, 68
      optional :min_distance_between_hits, :uint32, 69
      optional :max_distance_between_hits, :uint32, 70
      repeated :number_isotopic_atoms, :uint32, 71
      optional :min_number_isotopic_atoms, :uint32, 72
      optional :max_number_isotopic_atoms, :uint32, 73
      repeated :number_fragments, :uint32, 74
      optional :min_number_fragments, :uint32, 75
      optional :max_number_fragments, :uint32, 76
      repeated :distance_between_root_atoms, :uint32, 77
      optional :min_distance_between_root_atoms, :uint32, 78
      optional :max_distance_between_root_atoms, :uint32, 79
      repeated :atoms_in_spinach, :uint32, 80
      optional :min_atoms_in_spinach, :uint32, 81
      optional :max_atoms_in_spinach, :uint32, 82
      repeated :inter_ring_atoms, :uint32, 83
      optional :min_inter_ring_atoms, :uint32, 84
      optional :max_inter_ring_atoms, :uint32, 85
      repeated :unmatched_atoms, :uint32, 86
      optional :min_unmatched_atoms, :uint32, 87
      optional :max_unmatched_atoms, :uint32, 88
      repeated :net_formal_charge, :int32, 89
      optional :min_net_formal_charge, :int32, 90
      optional :max_net_formal_charge, :int32, 91
      optional :any_net_formal_charge, :bool, 92
      optional :min_fraction_atoms_matched, :float, 93
      optional :max_fraction_atoms_matched, :float, 94
      repeated :environment, :message, 95, "SubstructureSearch.SubstructureEnvironment"
      repeated :environment_no_match, :message, 96, "SubstructureSearch.SubstructureEnvironment"
      optional :environment_must_match_unmatched_atoms, :bool, 97
      optional :env_matches_can_share_attachment_points, :bool, 98
      repeated :matched_atom_must_be, :message, 122, "SubstructureSearch.MatchedAtomMatch"
      repeated :ring_specifier, :message, 99, "SubstructureSearch.SubstructureRingSpecification"
      repeated :ring_specification_logexp, :enum, 100, "SubstructureSearch.Operator"
      repeated :ring_system_specifier, :message, 101, "SubstructureSearch.SubstructureRingSystemSpecification"
      repeated :ring_system_specifier_logexp, :enum, 102, "SubstructureSearch.Operator"
      repeated :element_hits_needed, :message, 103, "SubstructureSearch.ElementsNeeded"
      repeated :elements_needed, :message, 104, "SubstructureSearch.ElementsNeeded"
      repeated :aromatic_atoms, :uint32, 105
      optional :min_aromatic_atoms, :uint32, 106
      optional :max_aromatic_atoms, :uint32, 107
      optional :unique_embeddings_only, :bool, 110
      repeated :heteroatoms, :uint32, 112
      optional :respect_initial_atom_numbering, :bool, 113
      optional :compress_embeddings, :bool, 114
      optional :environments_can_share_attachment_points, :bool, 115
      repeated :query_atom, :message, 116, "SubstructureSearch.SubstructureAtom"
      repeated :chiral_centre, :message, 117, "SubstructureSearch.SubstructureChiralCenter"
      optional :atom_type, :string, 119
      repeated :geometric_constraints, :message, 120, "GeometricConstraints.SetOfConstraints"
      repeated :separated_atoms, :message, 121, "SubstructureSearch.SeparatedAtoms"
      oneof :smiles_or_smarts do
        optional :smiles, :string, 17
        optional :smarts, :string, 18
      end
    end
    add_message "SubstructureSearch.SubstructureQuery" do
      optional :comment, :string, 1
      optional :name, :string, 2
      repeated :query, :message, 3, "SubstructureSearch.SingleSubstructureQuery"
      repeated :logexp, :enum, 4, "SubstructureSearch.Operator"
      optional :match_each_component, :int32, 5
    end
    add_message "SubstructureSearch.MinMaxSpecifierInt" do
      repeated :value, :int32, 1
      optional :min, :int32, 2
      optional :max, :int32, 3
    end
    add_message "SubstructureSearch.MinMaxSpecifierUInt" do
      repeated :value, :uint32, 1
      optional :min, :uint32, 2
      optional :max, :uint32, 3
    end
    add_message "SubstructureSearch.QueryMatchResults" do
      optional :smiles, :string, 1
      optional :name, :string, 2
      repeated :matches, :message, 3, "SubstructureSearch.QueryMatchResults.Matches"
    end
    add_message "SubstructureSearch.QueryMatchResults.Matches" do
      optional :name, :string, 1
      optional :nhits, :uint32, 2
    end
    add_enum "SubstructureSearch.Aromaticity" do
      value :SS_ALIPHATIC, 1
      value :SS_AROMATIC, 2
    end
    add_enum "SubstructureSearch.BondType" do
      value :SS_SINGLE_BOND, 3
      value :SS_DOUBLE_BOND, 4
      value :SS_TRIPLE_BOND, 5
      value :SS_AROMATIC_BOND, 6
      value :SS_BOND, 7
    end
    add_enum "SubstructureSearch.Operator" do
      value :SS_OR, 8
      value :SS_AND, 9
      value :SS_XOR, 10
      value :SS_LP_AND, 11
    end
  end
end

module SubstructureSearch
  AtomNumberHydrogenLonePair = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.AtomNumberHydrogenLonePair").msgclass
  AtomNumberHydrogenLonePair::HorLP = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.AtomNumberHydrogenLonePair.HorLP").enummodule
  SubstructureChiralCenter = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.SubstructureChiralCenter").msgclass
  SubstructureBond = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.SubstructureBond").msgclass
  SubstructureEnvironmentBond = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.SubstructureEnvironmentBond").msgclass
  ElementsNeeded = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.ElementsNeeded").msgclass
  NoMatchedAtomsBetween = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.NoMatchedAtomsBetween").msgclass
  LinkAtoms = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.LinkAtoms").msgclass
  EnvironmentAttachment = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.EnvironmentAttachment").msgclass
  SubstructureRingEnvironment = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.SubstructureRingEnvironment").msgclass
  SubstructureRingBase = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.SubstructureRingBase").msgclass
  SubstructureRingSpecification = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.SubstructureRingSpecification").msgclass
  RingSizeRequirement = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.RingSizeRequirement").msgclass
  SubstructureRingSystemSpecification = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.SubstructureRingSystemSpecification").msgclass
  SubstructureAtomSpecifier = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.SubstructureAtomSpecifier").msgclass
  SubstructureAtomEnvironment = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.SubstructureAtomEnvironment").msgclass
  SubstructureAtom = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.SubstructureAtom").msgclass
  SubstructureEnvironment = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.SubstructureEnvironment").msgclass
  MatchedAtomMatch = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.MatchedAtomMatch").msgclass
  SeparatedAtoms = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.SeparatedAtoms").msgclass
  SingleSubstructureQuery = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.SingleSubstructureQuery").msgclass
  SubstructureQuery = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.SubstructureQuery").msgclass
  MinMaxSpecifierInt = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.MinMaxSpecifierInt").msgclass
  MinMaxSpecifierUInt = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.MinMaxSpecifierUInt").msgclass
  QueryMatchResults = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.QueryMatchResults").msgclass
  QueryMatchResults::Matches = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.QueryMatchResults.Matches").msgclass
  Aromaticity = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.Aromaticity").enummodule
  BondType = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.BondType").enummodule
  Operator = ::Google::Protobuf::DescriptorPool.generated_pool.lookup("SubstructureSearch.Operator").enummodule
end