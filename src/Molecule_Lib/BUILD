package(default_visibility = ["//visibility:public"])

#genrule(
#  name = 'gfortran_iw',
#  srcs = [
#    'u3b.f',
#  ],
#  outs = ['u3b.o'],
#  cmd = 'gfortran -O3 -c -o $(OUTS) u3b.f'
#)

proto_library(
    name = "mol2graph_proto",
    srcs = [
        "mol2graph.proto",
    ],
)

cc_proto_library(
    name = "mol2graph_cc_proto",
    deps = [
        ":mol2graph_proto",
    ],
)

proto_library(
    name = "geometric_constraints_proto",
    srcs = [
        "geometric_constraints.proto",
    ],
)

cc_proto_library(
    name = "geometric_constraints_cc_proto",
    deps = [
        ":geometric_constraints_proto",
    ],
)


proto_library(
    name = "substructure_proto",
    srcs = [
        "substructure.proto"
    ],
    deps = [
        "geometric_constraints_proto",
    ],
)

cc_proto_library(
    name = "substructure_cc_proto",
    deps = [
        ":substructure_proto",
    ],
)

proto_library(
    name = "reaction_proto",
    srcs = [
        "reaction.proto",
    ],
    deps = [
       ":substructure_proto",
       ":toggle_kekule_form",
    ],
)

cc_proto_library(
    name = "reaction_cc_proto",
    deps = [
        "reaction_proto",
    ],
)

proto_library(
    name = "toggle_kekule_form",
    srcs = [
        "toggle_kekule_form.proto",
    ],
    deps = [
       ":substructure_proto",
    ],
)

cc_proto_library(
    name = "toggle_kekule_form_cc_proto",
    deps = [
        ":toggle_kekule_form",
    ],
)

cc_library(
    name = "coordinate_box",
    srcs = [
        "coordinate_box.cc",
    ],
    hdrs = [
        "coordinate_box.h",
    ],
    deps = [
        ":space_vector",
        "//Foundational/iwstring:iwstring",
    ],
)

cc_library(
    name = "iwmolecule",
    visibility = [
        "//visibility:public",
    ],
    srcs = [
        "allowed_elements.cc",
        "aromatic.cc",
        "atom.cc",
        "atom_alias.cc",
        "atom_typing.cc",
        "bond.cc",
        "bond_list.cc",
        "careful_frag.cc",
        "charge_assigner.cc",
        "charge_calculation.cc",
        "charmm.cc",
        "chiral_centre.cc",
        "cif.cc",
        "cis_trans_bond.cc",
        "coordinates.cc",
        "coordinate_box.cc",
        "count_atoms_in_smiles.cc",
        "csubstructure.cc",
        "dihedral.cc",
        "donor_acceptor.cc",
        "element.cc",
        "element_hits_needed.cc",
        "ematch.cc",
        "etrans.cc",
        "frag.cc",
        "geometric_constraints.cc",
        "ISIS_Atom_List.cc",
        "is_actually_chiral.cc",
        "iwrcb.cc",
        "iwqry_wstats.cc",
        "iwrnm.cc",
        "iwsubstructure.cc",
        "mdl_atom_record.cc",
        "mdl_file_data.cc",
        "mdl_molecule.cc",
        "misc2.cc",
        "mol2graph.cc",
        "mol2graph_proto.cc",
        "molecule.cc",
        "molecule_to_query.cc",
        "moleculeb.cc",
        "moleculed.cc",
        "moleculeh.cc",
        "moleculer.cc",
        "molecule_smarts.cc",
        "mpr.cc",
        "numass.cc",
        "ostream_and_type.cc",
        "output.cc",
        "parse_smarts_tmp.cc",
        "path.cc",
        "path_around_ring.cc",
        "path_scoring.cc",
        "pearlman.cc",
        "ring_bond_iterator.cc",
        "rmele.cc",
        "rotbond_common.cc",
        "rwsubstructure.cc",
        "rwsubstructure_proto.cc",
        "set_of_atoms.cc",
        "smiles.cc",
        "smiles_support.cc",
        "ss_atom_env.cc",
        "ss_bonds.cc",
        "ss_ring.cc",
        "ss_ring_base.cc",
        "ss_ring_sys.cc",
        "substructure_chiral.cc",
        "substructure_results.cc",
        "symm_class_can_rank.cc",
        "standardise.cc",
        "substructure_a.cc",
        "substructure_nmab.cc",
        "substructure_spec.cc",
        "substructure_env.cc",
        "symmetry.cc",
        "target.cc",
        "temp_detach_atoms.cc",
        "toggle_kekule_form.cc",
        "tokenise_atomic_smarts.cc",
        "unique.cc",
        "u3b.c",
        "_substructure.cc",
        "_qry_wstats.cc",
    ],
    hdrs = [
        "allowed_elements.h",
        "aromatic.h",
        "atom.h",
        "atom_alias.h",
        "atom_typing.h",
        "bond.h",
        "bond_list.h",
        "charge_assigner.h",
        "charge_calculation.h",
        "cis_trans_bond.h",
        "coordinates.h",
        "coordinate_box.h",
        "chiral_centre.h",
        "collection_template.h",
        "dihedral.h",
        "donor_acceptor.h",
        "element.h",
        "ematch.h",
        "etrans.h",
        "istream_and_type.h",
        "geometric_constraints.h",
        "is_actually_chiral.h",
        "iwmtypes.h",
        "iwrcb.h",
        "iwrnm.h",
        "standardise.h",
        "marvin.h",
        "mdl.h",
        "mdl_atom_record.h",
        "mdl_file_data.h",
        "mdl_molecule.h",
        "misc2.h",
        "mol2graph.h",
        "mol2graph_proto.h",
        "molecule.h",
        "molecule_to_query.h",
        "moleculej.h",
        "molecule_arom.h",
        "molecule_mdl.h",
        "molecule_smi.h",
        "mpr.h",
        "numass.h",
        "ostream_and_type.h",
        "output.h",
        "path.h",
        "path_around_ring.h",
        "path_scoring.h",
        "parse_smarts_tmp.h",
        "pearlman.h",
        "qry_wstats.h",
        "rotbond_common.h",
        "rmele.h",
        "rwsubstructure.h",
        "set_of_atoms.h",
        "smiles.h",
        "space_vector.h",
        "substructure.h",
        "symmetry.h",
        "target.h",
        "temp_detach_atoms.h",
        "tmpsssr.h",
        "toggle_kekule_form.h",
        "tokenise_atomic_smarts.h",
        "u3b.h",
    ],
    deps = [
        ":mol2graph_cc_proto",
        ":space_vector",
        ":substructure_cc_proto",
        ":toggle_kekule_form_cc_proto",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwconfig",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:iwre2",
        "//Foundational/iwbits:iwbits",
        "//Foundational/iwqsort:iwqsort",
        "//Foundational/iwstring:iwstring",
        "//Foundational/iw_tdt:iw_tdt",
        "//Foundational/xmlParser:xmlparser",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_library(
    name = "iwsubstructure",
    srcs = [
    ],
    hdrs = [
    ],
    deps = [
    ],
)

cc_library(
    name = "moleculeio",
    visibility = [
        "//visibility:public",
    ],
    srcs = [
        "charmm.cc",
        "cif.cc",
        "csv.cc",
        "inchi_dummy.cc",
        "marvin.cc",
        "mdl.cc",
        "mdl_v30.cc",
        "mmod.cc",
        "mrk.cc",
        "msi.cc",
        "pdb.cc",
        "rwmolecule.cc",
        "smi.cc",
        "tripos.cc",
    ],
    hdrs = [
        "readmdl.h",
        "rwmolecule.h",
        "tripos.h",
    ],
    deps = [
        ":iwmolecule",
    ],
    linkopts = [
    ],
)

cc_library(
    name = "iwreaction",
    srcs = [
        "iwreaction.cc",
        "iwreaction2.cc",
        "reaction_iterator.cc",
        "reaction_from_smirks.cc",
        "reaction_match_conditions.cc",
        "reaction_proto.cc",
        "rxnfile.cc",
        "rxnfile2.cc",
        "rxnfile3.cc",
        "rxn_3d_replace.cc",
        "rxn_bond_angle.cc",
        "rxn_bond_length.cc",
        "rxn_dihedral_angle.cc",
        "rxn_matched_sors.cc",
    ],
    hdrs = [
        "iwreaction.h",
        "reaction_match_conditions.h",
        "rxn_file.h",
    ],
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":reaction_cc_proto",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Foundational/iwmisc:iwmisc",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_library(
    name = "atom_pair_fingerprint",
    srcs = [
        "atom_pair_fingerprint.cc",
    ],
    hdrs = [
        "atom_pair_fingerprint.h",
    ],
    deps = [
        ":iwmolecule",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
)

cc_library(
    name = "ec_fingerprint",
    srcs = [
        "ec_fingerprint.cc",
    ],
    hdrs = [
        "ec_fingerprint.h",
    ],
    deps = [
        ":iwmolecule",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
)

cc_library(
    name = "geometric_constraints",
    srcs = [
        "geometric_constraints.cc",
    ],
    hdrs = [
        "geometric_constraints.h",
    ],
    deps = [
        ":iwmolecule",
        ":geometric_constraints_cc_proto",
    ],
)

cc_library(
    name = "linear_fingerprint",
    srcs = [
        "linear_fingerprint.cc",
    ],
    hdrs = [
        "linear_fingerprint.h",
    ],
    deps = [
        ":iwmolecule",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
)

cc_library(
    name = "space_vector",
    srcs = [
    ],
    hdrs = [
        "space_vector.h",
    ],
    deps = [
        "//Foundational/iwmisc:iwmisc",
    ],
)

cc_test(
    name = "aromatic_test",
    srcs = ["aromatic_test.cc"],
    timeout = "short",
    deps = [
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "atom_test",
    srcs = ["atom_test.cc"],
    timeout = "short",
    deps = [
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "linear_fingerprint_test",
    srcs = ["linear_fingerprint_test.cc"],
    timeout = "short",
    deps = [
        ":linear_fingerprint",
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "ec_fingerprint_test",
    srcs = ["ec_fingerprint_test.cc"],
    timeout = "short",
    deps = [
        ":ec_fingerprint",
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "frag_test",
    srcs = ["frag_test.cc"],
    timeout = "short",
    deps = [
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "geometric_constraints_test",
    srcs = ["geometric_constraints_test.cc"],
    timeout = "short",
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":geometric_constraints",
        ":geometric_constraints_cc_proto",
        "@googletest//:gtest_main",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "atom_pair_fingerprint_test",
    srcs = ["atom_pair_fingerprint_test.cc"],
    timeout = "short",
    deps = [
        ":atom_pair_fingerprint",
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "coordinate_box_test",
    srcs = [
        "coordinate_box_test.cc",
    ],
    deps = [
        ":coordinate_box",
        "@googletest//:gtest_main",
    ],
    timeout="short",
    linkopts = [
        "-lz -lm",
    ],
)


cc_test(
    name = "molecule_each_index_test",
    srcs = [
        "molecule_each_index_test.cc",
    ],
    timeout="short",
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest_main",
    ],
    linkopts = [
        "-lz",
    ]
)

cc_test(
    name = "ring_bond_count_test",
    srcs = [
        "ring_bond_count_test.cc",
    ],
    timeout="short",
    deps = [
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
    ],
    linkopts = [
        "-lz",
    ]
)

cc_test(
    name = "smi_test",
    srcs = [
        "smi_test.cc",
    ],
    timeout="short",
    deps = [
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "smiles_test",
    srcs = [
        "smiles_test.cc",
    ],
    timeout="short",
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest_main",
    ],
    linkopts = [
        "-lz",
    ]
)

cc_test(
    name = "space_vector_test",
    srcs = [
        "space_vector_test.cc",
    ],
    timeout="short",
    deps = [
        ":space_vector",
        "@googletest//:gtest_main",
    ],
    linkopts = [
        "-lz",
    ]
)

cc_test(
    name = "substructure_sepma_test",
    srcs = ["substructure_sepma_test.cc"],
    timeout="short",
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest_main",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "substructure_test",
    srcs = ["substructure_test.cc"],
    timeout="short",
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest_main",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "substructure_spec_test",
    srcs = ["substructure_spec_test.cc"],
    timeout="short",
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "substructure_single_query_test",
    srcs = ["substructure_single_query_test.cc"],
    timeout="short",
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "substructure_env_test",
    srcs = ["substructure_env_test.cc"],
    timeout="short",
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "substructure_ring_test",
    srcs = ["substructure_ring_test.cc"],
    timeout="short",
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "substructure_ring_system_test",
    srcs = ["substructure_ring_system_test.cc"],
    timeout="short",
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "substructure_smarts_test",
    srcs = ["substructure_smarts_test.cc"],
    timeout="short",
    deps = [
        ":iwmolecule",
        ":moleculeio",
        ":substructure_cc_proto",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "standardise_test",
    srcs = ["standardise_test.cc"],
    timeout="short",
    deps = [
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "reaction_test",
    srcs = ["iwreaction_test.cc"],
    timeout="short",
    deps = [
        ":iwmolecule",
        ":iwreaction",
        ":moleculeio",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "rxnfile_test",
    srcs = ["rxnfile_test.cc"],
    timeout="short",
    deps = [
        ":iwmolecule",
        ":iwreaction",
        ":moleculeio",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "rxnfile2_test",
    srcs = ["rxnfile2_test.cc"],
    timeout="short",
    deps = [
        ":iwmolecule",
        ":iwreaction",
        ":moleculeio",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "substructure_nmab_test",
    srcs = ["substructure_nmab_test.cc"],
    timeout="short",
    deps = [
        ":iwmolecule",
        ":iwreaction",
        ":moleculeio",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "symmetry_test",
    srcs = ["symmetry_test.cc"],
    timeout="short",
    deps = [
        ":iwmolecule",
        ":moleculeio",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "misc2_test",
    srcs = ["misc2_test.cc"],
    timeout = "short",
    deps = [
        ":iwmolecule",
        ":iwreaction",
        ":moleculeio",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "unique_test",
    srcs = ["unique_test.cc"],
    timeout = "short",
    deps = [
        ":iwmolecule",
        ":iwreaction",
        ":moleculeio",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_library(
    name = "iwmfingerprint",
    srcs = [
        "iwmfingerprint.cc",
    ],
    hdrs = [
        "iwmfingerprint.h",
    ],
    deps = [
        ":iwmolecule",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
)

cc_library(
    name = "iw_vdw",
    srcs = [
        "van_der_waals.cc",
    ],
    hdrs = [
        "iw_vdw.h",
    ],
    deps = [
        ":iwmolecule",
    ],
)
