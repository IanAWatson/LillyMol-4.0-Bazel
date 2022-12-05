# -*- coding: utf-8 -*-
# Generated by the protocol buffer compiler.  DO NOT EDIT!
# source: Molecule_Lib/substructure.proto
"""Generated protocol buffer code."""
from google.protobuf.internal import enum_type_wrapper
from google.protobuf import descriptor as _descriptor
from google.protobuf import descriptor_pool as _descriptor_pool
from google.protobuf import message as _message
from google.protobuf import reflection as _reflection
from google.protobuf import symbol_database as _symbol_database
# @@protoc_insertion_point(imports)

_sym_db = _symbol_database.Default()


from Molecule_Lib import geometric_constraints_pb2 as Molecule__Lib_dot_geometric__constraints__pb2


DESCRIPTOR = _descriptor_pool.Default().AddSerializedFile(b'\n\x1fMolecule_Lib/substructure.proto\x12\x12SubstructureSearch\x1a(Molecule_Lib/geometric_constraints.proto\"\xb1\x01\n\x1a\x41tomNumberHydrogenLonePair\x12\x15\n\x0b\x61tom_number\x18\x01 \x01(\rH\x00\x12G\n\x07h_or_lp\x18\x02 \x01(\x0e\x32\x34.SubstructureSearch.AtomNumberHydrogenLonePair.HorLPH\x00\"#\n\x05HorLP\x12\x0c\n\x08HYDROGEN\x10\x01\x12\x0c\n\x08LONEPAIR\x10\x02\x42\x0e\n\x0c\x41tomNumberOr\"\xb6\x02\n\x18SubstructureChiralCenter\x12\x0e\n\x06\x63\x65nter\x18\x01 \x01(\r\x12\x41\n\ttop_front\x18\x02 \x01(\x0b\x32..SubstructureSearch.AtomNumberHydrogenLonePair\x12@\n\x08top_back\x18\x03 \x01(\x0b\x32..SubstructureSearch.AtomNumberHydrogenLonePair\x12\x41\n\tleft_down\x18\x04 \x01(\x0b\x32..SubstructureSearch.AtomNumberHydrogenLonePair\x12\x42\n\nright_down\x18\x05 \x01(\x0b\x32..SubstructureSearch.AtomNumberHydrogenLonePair\"V\n\x10SubstructureBond\x12/\n\tbond_type\x18\x01 \x03(\x0e\x32\x1c.SubstructureSearch.BondType\x12\x11\n\tother_end\x18\x02 \x01(\r\"a\n\x1bSubstructureEnvironmentBond\x12/\n\tbond_type\x18\x01 \x03(\x0e\x32\x1c.SubstructureSearch.BondType\x12\x11\n\tother_end\x18\x02 \x03(\r\"\x9d\x01\n\x0e\x45lementsNeeded\x12\x17\n\ratomic_number\x18\x01 \x01(\rH\x00\x12\x17\n\ratomic_symbol\x18\x02 \x01(\tH\x00\x12\x13\n\x0bhits_needed\x18\x03 \x03(\r\x12\x17\n\x0fmin_hits_needed\x18\x04 \x01(\r\x12\x17\n\x0fmax_hits_needed\x18\x05 \x01(\rB\x12\n\x10\x45lementSpecifier\"B\n\x15NoMatchedAtomsBetween\x12\n\n\x02\x61\x31\x18\x01 \x01(\r\x12\n\n\x02\x61\x32\x18\x02 \x01(\r\x12\x11\n\tqualifier\x18\x03 \x01(\t\"a\n\tLinkAtoms\x12\n\n\x02\x61\x31\x18\x01 \x01(\r\x12\n\n\x02\x61\x32\x18\x02 \x01(\r\x12\x10\n\x08\x64istance\x18\x03 \x03(\r\x12\x14\n\x0cmin_distance\x18\x04 \x01(\r\x12\x14\n\x0cmax_distance\x18\x05 \x01(\r\"y\n\x15\x45nvironmentAttachment\x12\x18\n\x10\x61ttachment_point\x18\x01 \x03(\r\x12+\n\x05\x62type\x18\x02 \x03(\x0e\x32\x1c.SubstructureSearch.BondType\x12\x19\n\x11substructure_bond\x18\x03 \x01(\t\"\x90\x01\n\x1bSubstructureRingEnvironment\x12?\n\x11substructure_atom\x18\x01 \x01(\x0b\x32$.SubstructureSearch.SubstructureAtom\x12\x17\n\x0fmin_hits_needed\x18\x02 \x01(\r\x12\x17\n\x0fmax_hits_needed\x18\x03 \x01(\r\"\xf9\x07\n\x14SubstructureRingBase\x12\x16\n\x0ematch_as_match\x18\x01 \x01(\x08\x12\x13\n\x0bhits_needed\x18\x02 \x03(\r\x12\x17\n\x0fmin_hits_needed\x18\x03 \x01(\r\x12\x17\n\x0fmax_hits_needed\x18\x04 \x01(\r\x12!\n\x19\x61ttached_heteroatom_count\x18\x05 \x03(\r\x12%\n\x1dmin_attached_heteroatom_count\x18\x06 \x01(\r\x12%\n\x1dmax_attached_heteroatom_count\x18\x07 \x01(\r\x12\x18\n\x10heteroatom_count\x18\x08 \x03(\r\x12\x1c\n\x14min_heteroatom_count\x18\t \x01(\r\x12\x1c\n\x14max_heteroatom_count\x18\n \x01(\r\x12\x0c\n\x04ncon\x18\x0b \x03(\r\x12\x10\n\x08min_ncon\x18\x0c \x01(\r\x12\x10\n\x08max_ncon\x18\r \x01(\r\x12!\n\x19\x61ll_hits_in_same_fragment\x18\x0e \x01(\x08\x12 \n\x18within_ring_unsaturation\x18\x10 \x03(\r\x12$\n\x1cmin_within_ring_unsaturation\x18\x11 \x01(\r\x12$\n\x1cmax_within_ring_unsaturation\x18\x12 \x01(\r\x12\x38\n0largest_number_of_bonds_shared_with_another_ring\x18\x13 \x03(\r\x12<\n4min_largest_number_of_bonds_shared_with_another_ring\x18\x14 \x01(\r\x12<\n4max_largest_number_of_bonds_shared_with_another_ring\x18\x15 \x01(\r\x12\x1f\n\x17\x61toms_with_pi_electrons\x18\x1a \x03(\r\x12#\n\x1bmin_atoms_with_pi_electrons\x18\x1b \x01(\r\x12#\n\x1bmax_atoms_with_pi_electrons\x18\x1c \x01(\r\x12&\n\x1estrongly_fused_ring_neighbours\x18\x1d \x03(\r\x12*\n\"min_strongly_fused_ring_neighbours\x18\x1e \x01(\r\x12*\n\"max_strongly_fused_ring_neighbours\x18\x1f \x01(\r\x12\x13\n\x0b\x65nvironment\x18\x16 \x01(\t\x12+\n#environment_can_match_in_ring_atoms\x18\x17 \x01(\x08\x12\x15\n\rset_global_id\x18\x18 \x01(\r\"\xcd\x03\n\x1dSubstructureRingSpecification\x12\x36\n\x04\x62\x61se\x18\x01 \x01(\x0b\x32(.SubstructureSearch.SubstructureRingBase\x12\x11\n\tring_size\x18\x02 \x03(\r\x12\x15\n\rmin_ring_size\x18\x03 \x01(\r\x12\x15\n\rmax_ring_size\x18\x04 \x01(\r\x12\x10\n\x08\x61romatic\x18\x05 \x01(\x08\x12\r\n\x05\x66used\x18\x06 \x03(\r\x12\x11\n\tmin_fused\x18\x07 \x01(\r\x12\x11\n\tmax_fused\x18\x08 \x01(\r\x12!\n\x19\x66used_aromatic_neighbours\x18\t \x03(\r\x12%\n\x1dmin_fused_aromatic_neighbours\x18\n \x01(\r\x12%\n\x1dmax_fused_aromatic_neighbours\x18\x0b \x01(\r\x12%\n\x1d\x66used_non_aromatic_neighbours\x18\x0c \x03(\r\x12)\n!min_fused_non_aromatic_neighbours\x18\r \x01(\r\x12)\n!max_fused_non_aromatic_neighbours\x18\x0e \x01(\r\"]\n\x13RingSizeRequirement\x12\x11\n\tring_size\x18\x01 \x01(\r\x12\r\n\x05\x63ount\x18\x02 \x03(\r\x12\x11\n\tmin_count\x18\x03 \x01(\r\x12\x11\n\tmax_count\x18\x04 \x01(\r\"\xb1\n\n#SubstructureRingSystemSpecification\x12\x36\n\x04\x62\x61se\x18\x01 \x01(\x0b\x32(.SubstructureSearch.SubstructureRingBase\x12\x17\n\x0frings_in_system\x18\x02 \x03(\r\x12\x1b\n\x13min_rings_in_system\x18\x03 \x01(\r\x12\x1b\n\x13max_rings_in_system\x18\x04 \x01(\r\x12\x12\n\nring_sizes\x18\x05 \x03(\r\x12\x16\n\x0emin_ring_sizes\x18\x06 \x01(\r\x12\x16\n\x0emax_ring_sizes\x18\x07 \x01(\r\x12@\n\x0fring_size_count\x18\x08 \x03(\x0b\x32\'.SubstructureSearch.RingSizeRequirement\x12\x1b\n\x13\x61romatic_ring_count\x18\x0b \x03(\r\x12\x1f\n\x17min_aromatic_ring_count\x18\x0c \x01(\r\x12\x1f\n\x17max_aromatic_ring_count\x18\r \x01(\r\x12\x1f\n\x17non_aromatic_ring_count\x18\x0e \x03(\r\x12#\n\x1bmin_non_aromatic_ring_count\x18\x0f \x01(\r\x12#\n\x1bmax_non_aromatic_ring_count\x18\x10 \x01(\r\x12\x18\n\x10\x64\x65gree_of_fusion\x18\x11 \x03(\r\x12\x1c\n\x14min_degree_of_fusion\x18\x12 \x01(\r\x12\x1c\n\x14max_degree_of_fusion\x18\x13 \x01(\r\x12\x17\n\x0f\x61toms_in_system\x18\x14 \x03(\r\x12\x1b\n\x13min_atoms_in_system\x18\x15 \x01(\r\x12\x1b\n\x13max_atoms_in_system\x18\x16 \x01(\r\x12\x1d\n\x15number_spinach_groups\x18\x17 \x03(\r\x12!\n\x19min_number_spinach_groups\x18\x18 \x01(\r\x12!\n\x19max_number_spinach_groups\x18\x19 \x01(\r\x12!\n\x19number_non_spinach_groups\x18\x1a \x03(\r\x12%\n\x1dmin_number_non_spinach_groups\x18\x1b \x01(\r\x12%\n\x1dmax_number_non_spinach_groups\x18\x1c \x01(\r\x12\x1e\n\x16\x61toms_in_spinach_group\x18\x1d \x03(\r\x12\"\n\x1amin_atoms_in_spinach_group\x18\x1e \x01(\r\x12\"\n\x1amax_atoms_in_spinach_group\x18\x1f \x01(\r\x12\x1f\n\x17length_of_spinach_group\x18  \x03(\r\x12#\n\x1bmin_length_of_spinach_group\x18! \x01(\r\x12#\n\x1bmax_length_of_spinach_group\x18\" \x01(\r\x12 \n\x18\x64istance_to_another_ring\x18# \x03(\r\x12$\n\x1cmin_distance_to_another_ring\x18$ \x01(\r\x12$\n\x1cmax_distance_to_another_ring\x18% \x01(\r\x12!\n\x19strongly_fused_ring_count\x18& \x03(\r\x12%\n\x1dmin_strongly_fused_ring_count\x18\' \x01(\r\x12%\n\x1dmax_strongly_fused_ring_count\x18( \x01(\r\"\xa1\x0f\n\x19SubstructureAtomSpecifier\x12\x15\n\ratomic_symbol\x18\x01 \x03(\t\x12\x15\n\ratomic_number\x18\x02 \x03(\r\x12\x0c\n\x04ncon\x18\x03 \x03(\r\x12\x10\n\x08min_ncon\x18\x04 \x01(\r\x12\x10\n\x08max_ncon\x18\x05 \x01(\r\x12\r\n\x05ncon2\x18\x06 \x03(\r\x12\x11\n\tmin_ncon2\x18\x07 \x01(\r\x12\x11\n\tmax_ncon2\x18\x08 \x01(\r\x12\x0e\n\x06nbonds\x18\t \x03(\r\x12\x12\n\nmin_nbonds\x18\n \x01(\r\x12\x12\n\nmax_nbonds\x18\x0b \x01(\r\x12\x15\n\rformal_charge\x18\x0c \x03(\x05\x12\x19\n\x11min_formal_charge\x18\r \x01(\x05\x12\x19\n\x11max_formal_charge\x18\x0e \x01(\x05\x12\x0e\n\x06nrings\x18\x0f \x03(\r\x12\x12\n\nmin_nrings\x18\x10 \x01(\r\x12\x12\n\nmax_nrings\x18\x11 \x01(\r\x12\x17\n\x0fring_bond_count\x18\x12 \x03(\r\x12\x1b\n\x13min_ring_bond_count\x18\x13 \x01(\r\x12\x1b\n\x13max_ring_bond_count\x18\x14 \x01(\r\x12\x11\n\tring_size\x18\x15 \x03(\r\x12\x15\n\rmin_ring_size\x18\x16 \x01(\r\x12\x15\n\rmax_ring_size\x18\x17 \x01(\r\x12\x0e\n\x06hcount\x18\x18 \x03(\r\x12\x12\n\nmin_hcount\x18\x19 \x01(\r\x12\x12\n\nmax_hcount\x18\x1a \x01(\r\x12\x10\n\x08\x61romatic\x18\x1b \x01(\x08\x12\x11\n\tchirality\x18\x1c \x01(\x08\x12\x1a\n\x12\x61romatic_ring_size\x18\x1e \x03(\r\x12\x1e\n\x16min_aromatic_ring_size\x18\x1f \x01(\r\x12\x1e\n\x16max_aromatic_ring_size\x18  \x01(\r\x12\x1b\n\x13\x61liphatic_ring_size\x18! \x03(\r\x12\x1f\n\x17min_aliphatic_ring_size\x18\" \x01(\r\x12\x1f\n\x17max_aliphatic_ring_size\x18# \x01(\r\x12!\n\x19\x61ttached_heteroatom_count\x18$ \x03(\r\x12%\n\x1dmin_attached_heteroatom_count\x18% \x01(\r\x12%\n\x1dmax_attached_heteroatom_count\x18& \x01(\r\x12\x17\n\x0flone_pair_count\x18\' \x03(\r\x12\x1b\n\x13min_lone_pair_count\x18( \x01(\r\x12\x1b\n\x13max_lone_pair_count\x18) \x01(\r\x12\x14\n\x0cunsaturation\x18* \x03(\r\x12\x18\n\x10min_unsaturation\x18+ \x01(\r\x12\x18\n\x10max_unsaturation\x18, \x01(\r\x12\x12\n\ndaylight_x\x18- \x03(\r\x12\x16\n\x0emin_daylight_x\x18. \x01(\r\x12\x16\n\x0emax_daylight_x\x18/ \x01(\r\x12\x0f\n\x07isotope\x18\x30 \x03(\r\x12\x13\n\x0bmin_isotope\x18\x31 \x01(\r\x12\x13\n\x0bmax_isotope\x18\x32 \x01(\r\x12\x0c\n\x04\x61ryl\x18\x33 \x03(\r\x12\x10\n\x08min_aryl\x18\x34 \x01(\r\x12\x10\n\x08max_aryl\x18\x35 \x01(\r\x12\r\n\x05vinyl\x18\x36 \x03(\r\x12\x11\n\tmin_vinyl\x18\x37 \x01(\r\x12\x11\n\tmax_vinyl\x18\x38 \x01(\r\x12\x19\n\x11\x66used_system_size\x18\x39 \x03(\r\x12\x1d\n\x15min_fused_system_size\x18: \x01(\r\x12\x1d\n\x15max_fused_system_size\x18; \x01(\r\x12\x18\n\x10\x61ll_rings_kekule\x18< \x01(\x08\x12\x1b\n\x13heteroatoms_in_ring\x18= \x03(\r\x12\x1f\n\x17min_heteroatoms_in_ring\x18> \x01(\r\x12\x1f\n\x17max_heteroatoms_in_ring\x18? \x01(\r\x12\x1a\n\x12match_spinach_only\x18@ \x01(\x05\x12\'\n\x1fscaffold_bonds_attached_to_ring\x18\x41 \x03(\r\x12+\n#min_scaffold_bonds_attached_to_ring\x18\x42 \x01(\r\x12+\n#max_scaffold_bonds_attached_to_ring\x18\x43 \x01(\r\x12\x18\n\x10preference_value\x18\x44 \x01(\x05\x12\x17\n\x0fsymmetry_degree\x18\x45 \x03(\r\x12\x1b\n\x13min_symmetry_degree\x18\x46 \x01(\r\x12\x1b\n\x13max_symmetry_degree\x18G \x01(\r\x12\x16\n\x0esymmetry_group\x18H \x01(\x05\x12\x36\n\x10logical_operator\x18L \x01(\x0e\x32\x1c.SubstructureSearch.Operator\x12\x16\n\x0euser_atom_type\x18M \x01(\r\x12\x11\n\tatom_type\x18N \x01(\r\x12\x0f\n\x07valence\x18O \x03(\r\x12\x13\n\x0bmin_valence\x18P \x01(\r\x12\x13\n\x0bmax_valence\x18Q \x01(\r\"v\n\x1bSubstructureAtomEnvironment\x12\n\n\x02id\x18\x01 \x01(\r\x12?\n\x11substructure_atom\x18\x02 \x03(\x0b\x32$.SubstructureSearch.SubstructureAtom\x12\n\n\x02op\x18\x03 \x01(\t\"\x8b\x07\n\x10SubstructureAtom\x12\n\n\x02id\x18\x01 \x01(\x05\x12\x16\n\x0ematch_as_match\x18\x02 \x01(\x08\x12\x17\n\x0ftext_identifier\x18\x03 \x01(\t\x12\x17\n\x0f\x61tom_map_number\x18\x04 \x01(\r\x12\x1b\n\x13initial_atom_number\x18\x05 \x01(\r\x12\x46\n\x0f\x61tom_properties\x18\x07 \x03(\x0b\x32-.SubstructureSearch.SubstructureAtomSpecifier\x12\x0f\n\x07ring_id\x18\t \x01(\x05\x12\x17\n\x0f\x66used_system_id\x18\n \x01(\r\x12\x13\n\x0b\x66ragment_id\x18\x0b \x01(\x05\x12\x15\n\rnumeric_value\x18\x0c \x01(\x01\x12\x1c\n\x14include_in_embedding\x18\r \x01(\x08\x12\x10\n\x06smarts\x18\x0e \x01(\tH\x00\x12\x15\n\x0b\x61tom_smarts\x18\x0f \x01(\tH\x00\x12\x10\n\x06smiles\x18\x10 \x01(\tH\x00\x12\x44\n\x0b\x65nvironment\x18\x11 \x03(\x0b\x32/.SubstructureSearch.SubstructureAtomEnvironment\x12\x38\n\nquery_bond\x18\x15 \x03(\x0b\x32$.SubstructureSearch.SubstructureBond\x12\x13\n\x0b\x62ond_smarts\x18\x16 \x01(\t\x12\x13\n\x0bsingle_bond\x18\x19 \x03(\r\x12\x13\n\x0b\x64ouble_bond\x18\x1a \x03(\r\x12\x13\n\x0btriple_bond\x18\x1b \x03(\r\x12\x15\n\raromatic_bond\x18\x1c \x03(\r\x12\x0c\n\x04\x62ond\x18\x1d \x03(\r\x12\x41\n\npreference\x18\x17 \x03(\x0b\x32-.SubstructureSearch.SubstructureAtomSpecifier\x12\x1f\n\x17sum_all_preference_hits\x18\x18 \x01(\x08\x12 \n\x18unmatched_atoms_attached\x18\x1e \x03(\r\x12$\n\x1cmin_unmatched_atoms_attached\x18\x1f \x01(\r\x12$\n\x1cmax_unmatched_atoms_attached\x18  \x01(\r\x12\x17\n\x0f\x61tom_type_group\x18! \x01(\r\x12\x17\n\x0fglobal_match_id\x18\" \x01(\rB\x10\n\x0eSmilesOrSmarts\"\xe0\x03\n\x17SubstructureEnvironment\x12\n\n\x02id\x18\x01 \x01(\r\x12\x0e\n\x06smarts\x18\x03 \x03(\t\x12\x0e\n\x06smiles\x18\x04 \x03(\t\x12\x38\n\nquery_atom\x18\x05 \x03(\x0b\x32$.SubstructureSearch.SubstructureAtom\x12=\n\nattachment\x18\x06 \x01(\x0b\x32).SubstructureSearch.EnvironmentAttachment\x12\x0c\n\x04\x62ond\x18\x07 \x03(\t\x12\r\n\x05or_id\x18\x08 \x01(\r\x12\x0e\n\x06\x61nd_id\x18\t \x01(\r\x12\x13\n\x0bhits_needed\x18\n \x03(\r\x12\x17\n\x0fmin_hits_needed\x18\x0b \x01(\r\x12\x17\n\x0fmax_hits_needed\x18\x0c \x01(\r\x12%\n\x1dno_other_substituents_allowed\x18\r \x01(\x08\x12/\n\'env_matches_can_share_attachment_points\x18\x0f \x01(\x08\x12\x1b\n\x13max_matches_to_find\x18\x10 \x01(\r\x12\x13\n\x0bhydrogen_ok\x18\x11 \x01(\x08\x12\"\n\x1amax_env_matches_per_anchor\x18\x12 \x01(\r\"^\n\x10MatchedAtomMatch\x12\x0c\n\x04\x61tom\x18\x01 \x03(\x05\x12\x0e\n\x06smarts\x18\x02 \x03(\t\x12,\n\x06logexp\x18\x03 \x01(\x0e\x32\x1c.SubstructureSearch.Operator\"u\n\x0eSeparatedAtoms\x12\n\n\x02\x61\x31\x18\x01 \x01(\r\x12\n\n\x02\x61\x32\x18\x02 \x01(\r\x12\x15\n\rbonds_between\x18\x03 \x03(\r\x12\x19\n\x11min_bonds_between\x18\x04 \x01(\r\x12\x19\n\x11max_bonds_between\x18\x05 \x01(\r\"\x88 \n\x17SingleSubstructureQuery\x12\n\n\x02id\x18\x01 \x01(\x05\x12\r\n\x05label\x18\x02 \x01(\t\x12\x0f\n\x07\x63omment\x18\x03 \x01(\t\x12$\n\x1cone_embedding_per_start_atom\x18\x04 \x01(\x08\x12$\n\x1cnormalise_rc_per_hits_needed\x18\x05 \x01(\r\x12\x18\n\x10subtract_from_rc\x18\x06 \x01(\r\x12\x1b\n\x13max_matches_to_find\x18\x08 \x01(\r\x12\x1a\n\x12save_matched_atoms\x18\t \x01(\x08\x12$\n\x1cncon_ignore_singly_connected\x18\n \x01(\x08\x12&\n\x1eperceive_symmetric_equivalents\x18\x0b \x01(\x08\x12\x1f\n\x17implicit_ring_condition\x18\x0c \x01(\r\x12!\n\x19\x61ll_hits_in_same_fragment\x18\r \x01(\x08\x12#\n\x1bonly_match_largest_fragment\x18\x0e \x01(\x08\x12!\n\x19\x65mbeddings_do_not_overlap\x18\x0f \x01(\x08\x12 \n\x18sort_by_preference_value\x18\x10 \x01(\x08\x12\x10\n\x06smiles\x18\x11 \x01(\tH\x00\x12\x10\n\x06smarts\x18\x12 \x01(\tH\x00\x12\x15\n\rnumeric_value\x18\x13 \x03(\x01\x12K\n\x18no_matched_atoms_between\x18\x14 \x03(\x0b\x32).SubstructureSearch.NoMatchedAtomsBetween\x12+\n#no_matched_atoms_between_exhaustive\x18\x15 \x01(\x08\x12\x31\n\nlink_atoms\x18\x16 \x03(\x0b\x32\x1d.SubstructureSearch.LinkAtoms\x12$\n\x1c\x66\x61il_if_embeddings_too_close\x18\x17 \x01(\x08\x12$\n\x1c\x64istance_between_hits_ncheck\x18\x18 \x01(\r\x12\x14\n\x0csort_matches\x18\x19 \x01(\t\x12!\n\x19\x61ttached_heteroatom_count\x18\x1a \x03(\r\x12%\n\x1dmin_attached_heteroatom_count\x18\x1b \x01(\r\x12%\n\x1dmax_attached_heteroatom_count\x18\x1c \x01(\r\x12\x13\n\x0bhits_needed\x18\x1d \x03(\r\x12\x17\n\x0fmin_hits_needed\x18\x1e \x01(\r\x12\x17\n\x0fmax_hits_needed\x18\x1f \x01(\r\x12\x1a\n\x12ring_atoms_matched\x18  \x03(\r\x12\x1e\n\x16min_ring_atoms_matched\x18! \x01(\r\x12\x1e\n\x16max_ring_atoms_matched\x18\" \x01(\r\x12\x1b\n\x13heteroatoms_matched\x18# \x03(\r\x12\x1f\n\x17min_heteroatoms_matched\x18$ \x01(\r\x12\x1f\n\x17max_heteroatoms_matched\x18% \x01(\r\x12\x1f\n\x17heteroatoms_in_molecule\x18& \x03(\r\x12#\n\x1bmin_heteroatoms_in_molecule\x18\' \x01(\r\x12#\n\x1bmax_heteroatoms_in_molecule\x18( \x01(\r\x12\x0e\n\x06natoms\x18) \x03(\r\x12\x12\n\nmin_natoms\x18* \x01(\r\x12\x12\n\nmax_natoms\x18+ \x01(\r\x12\x0e\n\x06nrings\x18, \x03(\r\x12\x12\n\nmin_nrings\x18- \x01(\r\x12\x12\n\nmax_nrings\x18. \x01(\r\x12\x0c\n\x04ncon\x18/ \x03(\r\x12\x10\n\x08min_ncon\x18\x30 \x01(\r\x12\x10\n\x08max_ncon\x18\x31 \x01(\r\x12\x13\n\x0b\x66used_rings\x18\x32 \x03(\r\x12\x17\n\x0fmin_fused_rings\x18\x33 \x01(\r\x12\x17\n\x0fmax_fused_rings\x18\x34 \x01(\r\x12\x1c\n\x14strongly_fused_rings\x18\x35 \x03(\r\x12 \n\x18min_strongly_fused_rings\x18\x36 \x01(\r\x12 \n\x18max_strongly_fused_rings\x18\x37 \x01(\r\x12\x16\n\x0eisolated_rings\x18\x38 \x03(\r\x12\x1a\n\x12min_isolated_rings\x18\x39 \x01(\r\x12\x1a\n\x12max_isolated_rings\x18: \x01(\r\x12\x1d\n\x15isolated_ring_objects\x18; \x03(\r\x12!\n\x19min_isolated_ring_objects\x18< \x01(\r\x12!\n\x19max_isolated_ring_objects\x18= \x01(\r\x12\x16\n\x0e\x61romatic_rings\x18> \x03(\r\x12\x1a\n\x12min_aromatic_rings\x18? \x01(\r\x12\x1a\n\x12max_aromatic_rings\x18@ \x01(\r\x12\x1a\n\x12non_aromatic_rings\x18\x41 \x03(\r\x12\x1e\n\x16min_non_aromatic_rings\x18\x42 \x01(\r\x12\x1e\n\x16max_non_aromatic_rings\x18\x43 \x01(\r\x12\x1d\n\x15\x64istance_between_hits\x18\x44 \x03(\r\x12!\n\x19min_distance_between_hits\x18\x45 \x01(\r\x12!\n\x19max_distance_between_hits\x18\x46 \x01(\r\x12\x1d\n\x15number_isotopic_atoms\x18G \x03(\r\x12!\n\x19min_number_isotopic_atoms\x18H \x01(\r\x12!\n\x19max_number_isotopic_atoms\x18I \x01(\r\x12\x18\n\x10number_fragments\x18J \x03(\r\x12\x1c\n\x14min_number_fragments\x18K \x01(\r\x12\x1c\n\x14max_number_fragments\x18L \x01(\r\x12#\n\x1b\x64istance_between_root_atoms\x18M \x03(\r\x12\'\n\x1fmin_distance_between_root_atoms\x18N \x01(\r\x12\'\n\x1fmax_distance_between_root_atoms\x18O \x01(\r\x12\x18\n\x10\x61toms_in_spinach\x18P \x03(\r\x12\x1c\n\x14min_atoms_in_spinach\x18Q \x01(\r\x12\x1c\n\x14max_atoms_in_spinach\x18R \x01(\r\x12\x18\n\x10inter_ring_atoms\x18S \x03(\r\x12\x1c\n\x14min_inter_ring_atoms\x18T \x01(\r\x12\x1c\n\x14max_inter_ring_atoms\x18U \x01(\r\x12\x17\n\x0funmatched_atoms\x18V \x03(\r\x12\x1b\n\x13min_unmatched_atoms\x18W \x01(\r\x12\x1b\n\x13max_unmatched_atoms\x18X \x01(\r\x12\x19\n\x11net_formal_charge\x18Y \x03(\x05\x12\x1d\n\x15min_net_formal_charge\x18Z \x01(\x05\x12\x1d\n\x15max_net_formal_charge\x18[ \x01(\x05\x12\x1d\n\x15\x61ny_net_formal_charge\x18\\ \x01(\x08\x12\"\n\x1amin_fraction_atoms_matched\x18] \x01(\x02\x12\"\n\x1amax_fraction_atoms_matched\x18^ \x01(\x02\x12@\n\x0b\x65nvironment\x18_ \x03(\x0b\x32+.SubstructureSearch.SubstructureEnvironment\x12I\n\x14\x65nvironment_no_match\x18` \x03(\x0b\x32+.SubstructureSearch.SubstructureEnvironment\x12.\n&environment_must_match_unmatched_atoms\x18\x61 \x01(\x08\x12/\n\'env_matches_can_share_attachment_points\x18\x62 \x01(\x08\x12\x42\n\x14matched_atom_must_be\x18z \x03(\x0b\x32$.SubstructureSearch.MatchedAtomMatch\x12I\n\x0ering_specifier\x18\x63 \x03(\x0b\x32\x31.SubstructureSearch.SubstructureRingSpecification\x12?\n\x19ring_specification_logexp\x18\x64 \x03(\x0e\x32\x1c.SubstructureSearch.Operator\x12V\n\x15ring_system_specifier\x18\x65 \x03(\x0b\x32\x37.SubstructureSearch.SubstructureRingSystemSpecification\x12\x42\n\x1cring_system_specifier_logexp\x18\x66 \x03(\x0e\x32\x1c.SubstructureSearch.Operator\x12?\n\x13\x65lement_hits_needed\x18g \x03(\x0b\x32\".SubstructureSearch.ElementsNeeded\x12;\n\x0f\x65lements_needed\x18h \x03(\x0b\x32\".SubstructureSearch.ElementsNeeded\x12\x16\n\x0e\x61romatic_atoms\x18i \x03(\r\x12\x1a\n\x12min_aromatic_atoms\x18j \x01(\r\x12\x1a\n\x12max_aromatic_atoms\x18k \x01(\r\x12\x1e\n\x16unique_embeddings_only\x18n \x01(\x08\x12\x13\n\x0bheteroatoms\x18p \x03(\r\x12&\n\x1erespect_initial_atom_numbering\x18q \x01(\x08\x12\x1b\n\x13\x63ompress_embeddings\x18r \x01(\x08\x12\x30\n(environments_can_share_attachment_points\x18s \x01(\x08\x12\x38\n\nquery_atom\x18t \x03(\x0b\x32$.SubstructureSearch.SubstructureAtom\x12\x43\n\rchiral_centre\x18u \x03(\x0b\x32,.SubstructureSearch.SubstructureChiralCenter\x12\x11\n\tatom_type\x18w \x01(\t\x12\x45\n\x15geometric_constraints\x18x \x03(\x0b\x32&.GeometricConstraints.SetOfConstraints\x12;\n\x0fseparated_atoms\x18y \x03(\x0b\x32\".SubstructureSearch.SeparatedAtomsB\x12\n\x10smiles_or_smarts\"\xba\x01\n\x11SubstructureQuery\x12\x0f\n\x07\x63omment\x18\x01 \x01(\t\x12\x0c\n\x04name\x18\x02 \x01(\t\x12:\n\x05query\x18\x03 \x03(\x0b\x32+.SubstructureSearch.SingleSubstructureQuery\x12,\n\x06logexp\x18\x04 \x03(\x0e\x32\x1c.SubstructureSearch.Operator\x12\x1c\n\x14match_each_component\x18\x05 \x01(\x05\"=\n\x12MinMaxSpecifierInt\x12\r\n\x05value\x18\x01 \x03(\x05\x12\x0b\n\x03min\x18\x02 \x01(\x05\x12\x0b\n\x03max\x18\x03 \x01(\x05\">\n\x13MinMaxSpecifierUInt\x12\r\n\x05value\x18\x01 \x03(\r\x12\x0b\n\x03min\x18\x02 \x01(\r\x12\x0b\n\x03max\x18\x03 \x01(\r\"\x99\x01\n\x11QueryMatchResults\x12\x0e\n\x06smiles\x18\x01 \x01(\t\x12\x0c\n\x04name\x18\x02 \x01(\t\x12>\n\x07matches\x18\x03 \x03(\x0b\x32-.SubstructureSearch.QueryMatchResults.Matches\x1a&\n\x07Matches\x12\x0c\n\x04name\x18\x01 \x01(\t\x12\r\n\x05nhits\x18\x02 \x01(\r*0\n\x0b\x41romaticity\x12\x10\n\x0cSS_ALIPHATIC\x10\x01\x12\x0f\n\x0bSS_AROMATIC\x10\x02*i\n\x08\x42ondType\x12\x12\n\x0eSS_SINGLE_BOND\x10\x03\x12\x12\n\x0eSS_DOUBLE_BOND\x10\x04\x12\x12\n\x0eSS_TRIPLE_BOND\x10\x05\x12\x14\n\x10SS_AROMATIC_BOND\x10\x06\x12\x0b\n\x07SS_BOND\x10\x07*<\n\x08Operator\x12\t\n\x05SS_OR\x10\x08\x12\n\n\x06SS_AND\x10\t\x12\n\n\x06SS_XOR\x10\n\x12\r\n\tSS_LP_AND\x10\x0b')

_AROMATICITY = DESCRIPTOR.enum_types_by_name['Aromaticity']
Aromaticity = enum_type_wrapper.EnumTypeWrapper(_AROMATICITY)
_BONDTYPE = DESCRIPTOR.enum_types_by_name['BondType']
BondType = enum_type_wrapper.EnumTypeWrapper(_BONDTYPE)
_OPERATOR = DESCRIPTOR.enum_types_by_name['Operator']
Operator = enum_type_wrapper.EnumTypeWrapper(_OPERATOR)
SS_ALIPHATIC = 1
SS_AROMATIC = 2
SS_SINGLE_BOND = 3
SS_DOUBLE_BOND = 4
SS_TRIPLE_BOND = 5
SS_AROMATIC_BOND = 6
SS_BOND = 7
SS_OR = 8
SS_AND = 9
SS_XOR = 10
SS_LP_AND = 11


_ATOMNUMBERHYDROGENLONEPAIR = DESCRIPTOR.message_types_by_name['AtomNumberHydrogenLonePair']
_SUBSTRUCTURECHIRALCENTER = DESCRIPTOR.message_types_by_name['SubstructureChiralCenter']
_SUBSTRUCTUREBOND = DESCRIPTOR.message_types_by_name['SubstructureBond']
_SUBSTRUCTUREENVIRONMENTBOND = DESCRIPTOR.message_types_by_name['SubstructureEnvironmentBond']
_ELEMENTSNEEDED = DESCRIPTOR.message_types_by_name['ElementsNeeded']
_NOMATCHEDATOMSBETWEEN = DESCRIPTOR.message_types_by_name['NoMatchedAtomsBetween']
_LINKATOMS = DESCRIPTOR.message_types_by_name['LinkAtoms']
_ENVIRONMENTATTACHMENT = DESCRIPTOR.message_types_by_name['EnvironmentAttachment']
_SUBSTRUCTURERINGENVIRONMENT = DESCRIPTOR.message_types_by_name['SubstructureRingEnvironment']
_SUBSTRUCTURERINGBASE = DESCRIPTOR.message_types_by_name['SubstructureRingBase']
_SUBSTRUCTURERINGSPECIFICATION = DESCRIPTOR.message_types_by_name['SubstructureRingSpecification']
_RINGSIZEREQUIREMENT = DESCRIPTOR.message_types_by_name['RingSizeRequirement']
_SUBSTRUCTURERINGSYSTEMSPECIFICATION = DESCRIPTOR.message_types_by_name['SubstructureRingSystemSpecification']
_SUBSTRUCTUREATOMSPECIFIER = DESCRIPTOR.message_types_by_name['SubstructureAtomSpecifier']
_SUBSTRUCTUREATOMENVIRONMENT = DESCRIPTOR.message_types_by_name['SubstructureAtomEnvironment']
_SUBSTRUCTUREATOM = DESCRIPTOR.message_types_by_name['SubstructureAtom']
_SUBSTRUCTUREENVIRONMENT = DESCRIPTOR.message_types_by_name['SubstructureEnvironment']
_MATCHEDATOMMATCH = DESCRIPTOR.message_types_by_name['MatchedAtomMatch']
_SEPARATEDATOMS = DESCRIPTOR.message_types_by_name['SeparatedAtoms']
_SINGLESUBSTRUCTUREQUERY = DESCRIPTOR.message_types_by_name['SingleSubstructureQuery']
_SUBSTRUCTUREQUERY = DESCRIPTOR.message_types_by_name['SubstructureQuery']
_MINMAXSPECIFIERINT = DESCRIPTOR.message_types_by_name['MinMaxSpecifierInt']
_MINMAXSPECIFIERUINT = DESCRIPTOR.message_types_by_name['MinMaxSpecifierUInt']
_QUERYMATCHRESULTS = DESCRIPTOR.message_types_by_name['QueryMatchResults']
_QUERYMATCHRESULTS_MATCHES = _QUERYMATCHRESULTS.nested_types_by_name['Matches']
_ATOMNUMBERHYDROGENLONEPAIR_HORLP = _ATOMNUMBERHYDROGENLONEPAIR.enum_types_by_name['HorLP']
AtomNumberHydrogenLonePair = _reflection.GeneratedProtocolMessageType('AtomNumberHydrogenLonePair', (_message.Message,), {
  'DESCRIPTOR' : _ATOMNUMBERHYDROGENLONEPAIR,
  '__module__' : 'Molecule_Lib.substructure_pb2'
  # @@protoc_insertion_point(class_scope:SubstructureSearch.AtomNumberHydrogenLonePair)
  })
_sym_db.RegisterMessage(AtomNumberHydrogenLonePair)

SubstructureChiralCenter = _reflection.GeneratedProtocolMessageType('SubstructureChiralCenter', (_message.Message,), {
  'DESCRIPTOR' : _SUBSTRUCTURECHIRALCENTER,
  '__module__' : 'Molecule_Lib.substructure_pb2'
  # @@protoc_insertion_point(class_scope:SubstructureSearch.SubstructureChiralCenter)
  })
_sym_db.RegisterMessage(SubstructureChiralCenter)

SubstructureBond = _reflection.GeneratedProtocolMessageType('SubstructureBond', (_message.Message,), {
  'DESCRIPTOR' : _SUBSTRUCTUREBOND,
  '__module__' : 'Molecule_Lib.substructure_pb2'
  # @@protoc_insertion_point(class_scope:SubstructureSearch.SubstructureBond)
  })
_sym_db.RegisterMessage(SubstructureBond)

SubstructureEnvironmentBond = _reflection.GeneratedProtocolMessageType('SubstructureEnvironmentBond', (_message.Message,), {
  'DESCRIPTOR' : _SUBSTRUCTUREENVIRONMENTBOND,
  '__module__' : 'Molecule_Lib.substructure_pb2'
  # @@protoc_insertion_point(class_scope:SubstructureSearch.SubstructureEnvironmentBond)
  })
_sym_db.RegisterMessage(SubstructureEnvironmentBond)

ElementsNeeded = _reflection.GeneratedProtocolMessageType('ElementsNeeded', (_message.Message,), {
  'DESCRIPTOR' : _ELEMENTSNEEDED,
  '__module__' : 'Molecule_Lib.substructure_pb2'
  # @@protoc_insertion_point(class_scope:SubstructureSearch.ElementsNeeded)
  })
_sym_db.RegisterMessage(ElementsNeeded)

NoMatchedAtomsBetween = _reflection.GeneratedProtocolMessageType('NoMatchedAtomsBetween', (_message.Message,), {
  'DESCRIPTOR' : _NOMATCHEDATOMSBETWEEN,
  '__module__' : 'Molecule_Lib.substructure_pb2'
  # @@protoc_insertion_point(class_scope:SubstructureSearch.NoMatchedAtomsBetween)
  })
_sym_db.RegisterMessage(NoMatchedAtomsBetween)

LinkAtoms = _reflection.GeneratedProtocolMessageType('LinkAtoms', (_message.Message,), {
  'DESCRIPTOR' : _LINKATOMS,
  '__module__' : 'Molecule_Lib.substructure_pb2'
  # @@protoc_insertion_point(class_scope:SubstructureSearch.LinkAtoms)
  })
_sym_db.RegisterMessage(LinkAtoms)

EnvironmentAttachment = _reflection.GeneratedProtocolMessageType('EnvironmentAttachment', (_message.Message,), {
  'DESCRIPTOR' : _ENVIRONMENTATTACHMENT,
  '__module__' : 'Molecule_Lib.substructure_pb2'
  # @@protoc_insertion_point(class_scope:SubstructureSearch.EnvironmentAttachment)
  })
_sym_db.RegisterMessage(EnvironmentAttachment)

SubstructureRingEnvironment = _reflection.GeneratedProtocolMessageType('SubstructureRingEnvironment', (_message.Message,), {
  'DESCRIPTOR' : _SUBSTRUCTURERINGENVIRONMENT,
  '__module__' : 'Molecule_Lib.substructure_pb2'
  # @@protoc_insertion_point(class_scope:SubstructureSearch.SubstructureRingEnvironment)
  })
_sym_db.RegisterMessage(SubstructureRingEnvironment)

SubstructureRingBase = _reflection.GeneratedProtocolMessageType('SubstructureRingBase', (_message.Message,), {
  'DESCRIPTOR' : _SUBSTRUCTURERINGBASE,
  '__module__' : 'Molecule_Lib.substructure_pb2'
  # @@protoc_insertion_point(class_scope:SubstructureSearch.SubstructureRingBase)
  })
_sym_db.RegisterMessage(SubstructureRingBase)

SubstructureRingSpecification = _reflection.GeneratedProtocolMessageType('SubstructureRingSpecification', (_message.Message,), {
  'DESCRIPTOR' : _SUBSTRUCTURERINGSPECIFICATION,
  '__module__' : 'Molecule_Lib.substructure_pb2'
  # @@protoc_insertion_point(class_scope:SubstructureSearch.SubstructureRingSpecification)
  })
_sym_db.RegisterMessage(SubstructureRingSpecification)

RingSizeRequirement = _reflection.GeneratedProtocolMessageType('RingSizeRequirement', (_message.Message,), {
  'DESCRIPTOR' : _RINGSIZEREQUIREMENT,
  '__module__' : 'Molecule_Lib.substructure_pb2'
  # @@protoc_insertion_point(class_scope:SubstructureSearch.RingSizeRequirement)
  })
_sym_db.RegisterMessage(RingSizeRequirement)

SubstructureRingSystemSpecification = _reflection.GeneratedProtocolMessageType('SubstructureRingSystemSpecification', (_message.Message,), {
  'DESCRIPTOR' : _SUBSTRUCTURERINGSYSTEMSPECIFICATION,
  '__module__' : 'Molecule_Lib.substructure_pb2'
  # @@protoc_insertion_point(class_scope:SubstructureSearch.SubstructureRingSystemSpecification)
  })
_sym_db.RegisterMessage(SubstructureRingSystemSpecification)

SubstructureAtomSpecifier = _reflection.GeneratedProtocolMessageType('SubstructureAtomSpecifier', (_message.Message,), {
  'DESCRIPTOR' : _SUBSTRUCTUREATOMSPECIFIER,
  '__module__' : 'Molecule_Lib.substructure_pb2'
  # @@protoc_insertion_point(class_scope:SubstructureSearch.SubstructureAtomSpecifier)
  })
_sym_db.RegisterMessage(SubstructureAtomSpecifier)

SubstructureAtomEnvironment = _reflection.GeneratedProtocolMessageType('SubstructureAtomEnvironment', (_message.Message,), {
  'DESCRIPTOR' : _SUBSTRUCTUREATOMENVIRONMENT,
  '__module__' : 'Molecule_Lib.substructure_pb2'
  # @@protoc_insertion_point(class_scope:SubstructureSearch.SubstructureAtomEnvironment)
  })
_sym_db.RegisterMessage(SubstructureAtomEnvironment)

SubstructureAtom = _reflection.GeneratedProtocolMessageType('SubstructureAtom', (_message.Message,), {
  'DESCRIPTOR' : _SUBSTRUCTUREATOM,
  '__module__' : 'Molecule_Lib.substructure_pb2'
  # @@protoc_insertion_point(class_scope:SubstructureSearch.SubstructureAtom)
  })
_sym_db.RegisterMessage(SubstructureAtom)

SubstructureEnvironment = _reflection.GeneratedProtocolMessageType('SubstructureEnvironment', (_message.Message,), {
  'DESCRIPTOR' : _SUBSTRUCTUREENVIRONMENT,
  '__module__' : 'Molecule_Lib.substructure_pb2'
  # @@protoc_insertion_point(class_scope:SubstructureSearch.SubstructureEnvironment)
  })
_sym_db.RegisterMessage(SubstructureEnvironment)

MatchedAtomMatch = _reflection.GeneratedProtocolMessageType('MatchedAtomMatch', (_message.Message,), {
  'DESCRIPTOR' : _MATCHEDATOMMATCH,
  '__module__' : 'Molecule_Lib.substructure_pb2'
  # @@protoc_insertion_point(class_scope:SubstructureSearch.MatchedAtomMatch)
  })
_sym_db.RegisterMessage(MatchedAtomMatch)

SeparatedAtoms = _reflection.GeneratedProtocolMessageType('SeparatedAtoms', (_message.Message,), {
  'DESCRIPTOR' : _SEPARATEDATOMS,
  '__module__' : 'Molecule_Lib.substructure_pb2'
  # @@protoc_insertion_point(class_scope:SubstructureSearch.SeparatedAtoms)
  })
_sym_db.RegisterMessage(SeparatedAtoms)

SingleSubstructureQuery = _reflection.GeneratedProtocolMessageType('SingleSubstructureQuery', (_message.Message,), {
  'DESCRIPTOR' : _SINGLESUBSTRUCTUREQUERY,
  '__module__' : 'Molecule_Lib.substructure_pb2'
  # @@protoc_insertion_point(class_scope:SubstructureSearch.SingleSubstructureQuery)
  })
_sym_db.RegisterMessage(SingleSubstructureQuery)

SubstructureQuery = _reflection.GeneratedProtocolMessageType('SubstructureQuery', (_message.Message,), {
  'DESCRIPTOR' : _SUBSTRUCTUREQUERY,
  '__module__' : 'Molecule_Lib.substructure_pb2'
  # @@protoc_insertion_point(class_scope:SubstructureSearch.SubstructureQuery)
  })
_sym_db.RegisterMessage(SubstructureQuery)

MinMaxSpecifierInt = _reflection.GeneratedProtocolMessageType('MinMaxSpecifierInt', (_message.Message,), {
  'DESCRIPTOR' : _MINMAXSPECIFIERINT,
  '__module__' : 'Molecule_Lib.substructure_pb2'
  # @@protoc_insertion_point(class_scope:SubstructureSearch.MinMaxSpecifierInt)
  })
_sym_db.RegisterMessage(MinMaxSpecifierInt)

MinMaxSpecifierUInt = _reflection.GeneratedProtocolMessageType('MinMaxSpecifierUInt', (_message.Message,), {
  'DESCRIPTOR' : _MINMAXSPECIFIERUINT,
  '__module__' : 'Molecule_Lib.substructure_pb2'
  # @@protoc_insertion_point(class_scope:SubstructureSearch.MinMaxSpecifierUInt)
  })
_sym_db.RegisterMessage(MinMaxSpecifierUInt)

QueryMatchResults = _reflection.GeneratedProtocolMessageType('QueryMatchResults', (_message.Message,), {

  'Matches' : _reflection.GeneratedProtocolMessageType('Matches', (_message.Message,), {
    'DESCRIPTOR' : _QUERYMATCHRESULTS_MATCHES,
    '__module__' : 'Molecule_Lib.substructure_pb2'
    # @@protoc_insertion_point(class_scope:SubstructureSearch.QueryMatchResults.Matches)
    })
  ,
  'DESCRIPTOR' : _QUERYMATCHRESULTS,
  '__module__' : 'Molecule_Lib.substructure_pb2'
  # @@protoc_insertion_point(class_scope:SubstructureSearch.QueryMatchResults)
  })
_sym_db.RegisterMessage(QueryMatchResults)
_sym_db.RegisterMessage(QueryMatchResults.Matches)

if _descriptor._USE_C_DESCRIPTORS == False:

  DESCRIPTOR._options = None
  _AROMATICITY._serialized_start=12548
  _AROMATICITY._serialized_end=12596
  _BONDTYPE._serialized_start=12598
  _BONDTYPE._serialized_end=12703
  _OPERATOR._serialized_start=12705
  _OPERATOR._serialized_end=12765
  _ATOMNUMBERHYDROGENLONEPAIR._serialized_start=98
  _ATOMNUMBERHYDROGENLONEPAIR._serialized_end=275
  _ATOMNUMBERHYDROGENLONEPAIR_HORLP._serialized_start=224
  _ATOMNUMBERHYDROGENLONEPAIR_HORLP._serialized_end=259
  _SUBSTRUCTURECHIRALCENTER._serialized_start=278
  _SUBSTRUCTURECHIRALCENTER._serialized_end=588
  _SUBSTRUCTUREBOND._serialized_start=590
  _SUBSTRUCTUREBOND._serialized_end=676
  _SUBSTRUCTUREENVIRONMENTBOND._serialized_start=678
  _SUBSTRUCTUREENVIRONMENTBOND._serialized_end=775
  _ELEMENTSNEEDED._serialized_start=778
  _ELEMENTSNEEDED._serialized_end=935
  _NOMATCHEDATOMSBETWEEN._serialized_start=937
  _NOMATCHEDATOMSBETWEEN._serialized_end=1003
  _LINKATOMS._serialized_start=1005
  _LINKATOMS._serialized_end=1102
  _ENVIRONMENTATTACHMENT._serialized_start=1104
  _ENVIRONMENTATTACHMENT._serialized_end=1225
  _SUBSTRUCTURERINGENVIRONMENT._serialized_start=1228
  _SUBSTRUCTURERINGENVIRONMENT._serialized_end=1372
  _SUBSTRUCTURERINGBASE._serialized_start=1375
  _SUBSTRUCTURERINGBASE._serialized_end=2392
  _SUBSTRUCTURERINGSPECIFICATION._serialized_start=2395
  _SUBSTRUCTURERINGSPECIFICATION._serialized_end=2856
  _RINGSIZEREQUIREMENT._serialized_start=2858
  _RINGSIZEREQUIREMENT._serialized_end=2951
  _SUBSTRUCTURERINGSYSTEMSPECIFICATION._serialized_start=2954
  _SUBSTRUCTURERINGSYSTEMSPECIFICATION._serialized_end=4283
  _SUBSTRUCTUREATOMSPECIFIER._serialized_start=4286
  _SUBSTRUCTUREATOMSPECIFIER._serialized_end=6239
  _SUBSTRUCTUREATOMENVIRONMENT._serialized_start=6241
  _SUBSTRUCTUREATOMENVIRONMENT._serialized_end=6359
  _SUBSTRUCTUREATOM._serialized_start=6362
  _SUBSTRUCTUREATOM._serialized_end=7269
  _SUBSTRUCTUREENVIRONMENT._serialized_start=7272
  _SUBSTRUCTUREENVIRONMENT._serialized_end=7752
  _MATCHEDATOMMATCH._serialized_start=7754
  _MATCHEDATOMMATCH._serialized_end=7848
  _SEPARATEDATOMS._serialized_start=7850
  _SEPARATEDATOMS._serialized_end=7967
  _SINGLESUBSTRUCTUREQUERY._serialized_start=7970
  _SINGLESUBSTRUCTUREQUERY._serialized_end=12074
  _SUBSTRUCTUREQUERY._serialized_start=12077
  _SUBSTRUCTUREQUERY._serialized_end=12263
  _MINMAXSPECIFIERINT._serialized_start=12265
  _MINMAXSPECIFIERINT._serialized_end=12326
  _MINMAXSPECIFIERUINT._serialized_start=12328
  _MINMAXSPECIFIERUINT._serialized_end=12390
  _QUERYMATCHRESULTS._serialized_start=12393
  _QUERYMATCHRESULTS._serialized_end=12546
  _QUERYMATCHRESULTS_MATCHES._serialized_start=12508
  _QUERYMATCHRESULTS_MATCHES._serialized_end=12546
# @@protoc_insertion_point(module_scope)
