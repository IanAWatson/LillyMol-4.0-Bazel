
proto_library(
    name = "embedded_fragment_proto",
    srcs = [
        "embedded_fragment.proto",
    ],
)

cc_proto_library(
    name = "embedded_fragment_cc_proto",
    deps = [
        ":embedded_fragment_proto",
    ],
)


proto_library(
    name = "molecule_database_options_proto",
    srcs = [
        "molecule_database_options.proto",
    ],
    deps = [
        "//Molecule_Lib:mol2graph_proto",
    ]
)

cc_proto_library(
    name = "molecule_database_options_cc_proto",
    deps = [
        ":molecule_database_options_proto",
    ],
)

proto_library(
    name = "substructure_results_proto",
    srcs = [
        "substructure_results.proto",
    ],
    deps = [
    ],
)

cc_proto_library(
    name = "substructure_results_cc_proto",
    deps = [
        ":substructure_results_proto",
    ],
)

cc_binary(
    name = "activity_consistency",
    srcs = [
        "activity_consistency.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
        "//Foundational/accumulator:accumulator",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "atom_pair_fingerprint",
    srcs = [
        "atom_pair_fingerprint_main.cc",
    ],
    deps = [
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:atom_pair_fingerprint",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "buildsmidb_leveldb",
    srcs = [
        "buildsmidb_leveldb.cc",
    ],
    deps = [
        ":molecule_database_options",
        ":molecule_database_options_cc_proto",
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@leveldb//:leveldb",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "common_names",
    srcs = [
        "common_names.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline_v2:cmdline_v2",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "dbf",
    srcs = [
        "dbf.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline_v2:cmdline_v2",
        "//Foundational/histogram:iwhistogram",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "dicer",
    srcs = [
        "dicer.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline_v2:cmdline_v2",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "ec_fingerprint",
    srcs = [
        "ec_fingerprint_main.cc",
    ],
    deps = [
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:ec_fingerprint",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "echoqry",
    srcs = [
        "echoqry.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
    linkopts = [
        "-lz"
    ],
)

cc_binary(
    name = "fileconv",
    srcs = [
        "fileconv.cc",
    ],
    deps = [
        ":fileconv_lib",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline_v2:cmdline_v2",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "fragment_replacement",
    srcs = [
        "fragment_replacement.cc",
    ],
    deps = [
        ":atom_separations",
        ":embedded_fragment_cc_proto",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@com_google_protobuf//:protobuf",
        "@leveldb//:leveldb",
    ],
)

cc_binary(
    name = "generate_atom_separations",
    srcs = [
        "generate_atom_separations.cc",
    ],
    deps = [
        ":atom_separations",
        ":embedded_fragment_cc_proto",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@leveldb//:leveldb",
    ],
)

cc_binary(
    name = "generate_atom_separations_lmdb",
    srcs = [
        "generate_atom_separations_lmdb.cc",
    ],
    deps = [
        ":atom_separations",
        ":embedded_fragment_cc_proto",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
    linkopts = [
        "-llmdb",
    ],
)

cc_binary(
    name = "geometric_descriptors",
    srcs = [
        "geometric_descriptors.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
    linkopts = [
        "-lz",
    ]
)

cc_binary(
    name = "hydrophobic_sections",
    srcs = [
        "hydrophobic_sections.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray:iwaray",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
    linkopts = [
        "-lz"
    ],
)


cc_binary(
    name = "in_database_leveldb",
    srcs = [
        "in_database_leveldb.cc",
    ],
    deps = [
        ":molecule_database_options",
        ":molecule_database_options_cc_proto",
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@leveldb//:leveldb",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "iwdescr",
    srcs = [
        "iwdescr.cc",
        "qry_wcharge.cc",
        "qry_wcharge.h",
    ],
    deps = [
        ":nvrtspsa",
        ":partial_symmetry",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline_v2:cmdline_v2",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "iwecfp",
    srcs = [
        "iwecfp.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "iwfp",
    srcs = [
        "iwfp.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline_v2:cmdline_v2",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "linear_fingerprint",
    srcs = [
        "linear_fingerprint_main.cc",
    ],
    deps = [
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:linear_fingerprint",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "maccskeys_fn5",
    srcs = [
        "maccskeys_fn5_main.cc",
    ],
    deps = [
        ":maccskeys_fn5_lib",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:iwreaction",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "make_these_molecules",
    srcs = [
        "make_these_molecules.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:iwreaction",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "mc_summarise",
    srcs = [
        "mc_summarise.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_binary(
    name = "mol2qry",
    srcs = [
        "mol2qry.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@com_google_protobuf//:protobuf",
    ],
    linkopts = [
        "-lz",
    ],
)


cc_binary(
    name = "molecular_abstractions",
    srcs = [
        "mav2.cc",
        "molecular_abstraction_functions.cc",
        "molecular_abstraction_functions.h",
        "molecular_abstraction_specifications.cc",
        "molecular_abstraction_specifications.h",
    ],
    deps = [
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:iwreaction",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "molecular_scaffold",
    srcs = [
        "molecular_scaffold.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:iwreaction",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "molecular_transformations",
    srcs = [
        "molecular_transformations.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:iwreaction",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "molecular_variants",
    srcs = [
        "molecular_variants.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:iwreaction",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "molecule_subset",
    srcs = [
        "molecule_subset.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:iwreaction",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "msort",
    srcs = [
        "msort.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt:iw_tdt",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "pharmacaphore_2d",
    srcs = [
        "pharmacaphore_2d.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Foundational/cmdline:iwcmdline",
        "@com_google_protobuf//:protobuf",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "preferred_smiles",
    srcs = [
        "preferred_smiles.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt:iw_tdt",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "random_molecular_permutations",
    srcs = [
        "random_molecular_permutations.cc",
        "set_of_target_molecules.cc",
        "set_of_target_molecules.h",
        "swap_kekule_forms.cc",
    ],
    deps = [
        ":linear_path_fingerprint",
        ":maccskeys_fn5_lib",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Foundational/iw_tdt:iw_tdt",
        "//Utilities/GFP_Tools:gfp_standard",
        "//Utilities/GFP_Tools:sparse_collection",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "random_smiles",
    srcs = [
        "random_smiles.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt:iw_tdt",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "remove_and_label",
    srcs = [
        "remove_and_label.cc"
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "retrosynthesis",
    srcs = [
        "retrosynthesis.cc"
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "rgroup",
    srcs = [
        "rgroup.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt:iw_tdt",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "ring_substitution",
    srcs = [
        "ring_substitution.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Foundational/iw_tdt:iw_tdt",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "ring_extraction",
    srcs = [
        "ring_extraction.cc"
    ],
    deps = [
        ":reaction_duplicate",
        ":ring_ext_rep",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
        "@com_google_protobuf//:protobuf",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "ring_replacement",
    srcs = [
        "ring_replacement.cc"
    ],
    deps = [
        ":reaction_duplicate",
        ":ring_ext_rep",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "ring_trimming",
    srcs = [
        "ring_trimming.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt:iw_tdt",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "rotatable_bonds",
    srcs = [
        "rotatable_bonds.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Foundational/cmdline:iwcmdline",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "rotatable_bond_fingerprint",
    srcs = [
        "rotatable_bond_fingerprint.cc",
    ],
    deps = [
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Foundational/cmdline:iwcmdline",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "rule_of_five",
    srcs = [
        "rule_of_five.cc",
    ],
    deps = [
        ":nvrtspsa",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Foundational/cmdline:iwcmdline",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "rxn_fingerprint",
    srcs = [
        "rxn_fingerprint.cc",
        "extended_connectivity_fp.cc",
        "extended_connectivity_fp.h",
        "linear_path_fingerprint.cc",
        "linear_path_fingerprint.h",
        "ct_cache.h",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Foundational/iw_tdt:iw_tdt",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "rxn_reverse",
    srcs = [
        "rxn_reverse.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Foundational/iw_tdt:iw_tdt",
        "//Molecule_Lib:iwmfingerprint",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "rxn_signature",
    srcs = [
        "rxn_signature.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "rxn_substructure_search",
    srcs = [
        "rxn_substructure_search.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "rxn_standardize",
    srcs = [
        "rxn_standardize.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
    ],
    linkopts = [
        "-lz",
    ],
)

#cc_library(
#    name = "rxn_to_openrxn_lib",
#    srcs = [
#        "rxn_to_openrxn.cc",
#        "rxn_to_openrxn.h",
#    ],
#    deps = [
#        "//Foundational/cmdline:iwcmdline",
#        "//Foundational/iwaray:iwaray",
#        "//Foundational/iwmisc:iwmisc",
#        "//Molecule_Lib:iwmolecule",
#        "//Molecule_Lib:iwreaction",
#        "//Molecule_Lib:moleculeio",
#        "//Molecule_Lib:substructure_cc_proto",
#        "@ord_schema//:dataset_cc_proto",
#        "@ord_schema//:reaction_cc_proto",
#    ],
#)

#cc_binary(
#    name = "rxn_to_openrxn",
#    srcs = [
#        "rxn_to_openrxn_main.cc"
#    ],
#    deps = [
#        ":rxn_to_openrxn_lib",
#        "@ord_schema//:dataset_cc_proto",
#        "@ord_schema//:reaction_cc_proto",
#    ],
#    linkopts = [
#        "-lz",
#    ]
#)

cc_binary(
    name = "smiles_mutation",
    srcs = [
        "smiles_mutation.cc",
        "random_reactions.cc",
        "random_reactions.h",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "sp3_filter",
    srcs = [
        "sp3_filter.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt:iw_tdt",
    ],
    linkopts = [
        "-lz",
    ],
)

#cc_library(
#    name = "surface_area_molvol",
#    srcs = [
#        "surface_area_molvol.cc",
#    ],
#    hdrs = [
#        "surface_area_molvol.h"
#    ],
#    deps = [
#        "//Molecule_Lib:iwmolecule",
#        "//Molecule_Lib:moleculeio",
#        "//Molecule_Lib:iw_vdw",
#    ],
#)

cc_binary(
    name = "tautomer_generation",
    srcs = [
        "tautomer_generation.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "temperature",
    srcs = [
        "temperature.cc",
    ],
    deps = [
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:linear_fingerprint",
    ],
    linkopts = [
        "-lz",
    ]
)

cc_binary(
    name = "tnass",
    srcs = [
        "tnass.cc",
        "nass.cc",
        "nass.h",
    ],
    deps = [
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:linear_fingerprint",
    ],
    linkopts = [
        "-lz",
    ]
)

cc_binary(
    name = "topological_torsion",
    srcs = [
        "topological_torsion_main.cc",
    ],
    deps = [
        ":topological_torsion_lib",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "tp_first_pass",
    srcs = [
        "tp_first_pass.cc",
    ],
    deps = [
        ":tp_first_pass_lib",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Foundational/cmdline:iwcmdline",
    ],
    linkopts = [
        "-lz"
    ],
)

cc_binary(
    name = "trxn",
    srcs = [
        "trxn.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:iwreaction",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "tsmiles",
    srcs = [
        "tsmiles.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "tsubstructure",
    srcs = [
        "tsubstructure.cc",
        "tsubstructure_fp.cc",
        "tsubstructure_fp.h",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline_v2:cmdline_v2",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Foundational/iw_tdt:iw_tdt",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "unique_molecules",
    srcs = [
        "unique_molecules.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:iwreaction",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_library(
    name = "atom_separations",
    srcs = [
        "atom_separations.cc",
    ],
    hdrs = [
        "atom_separations.h",
    ],
    deps = [
        "//Foundational/iwstring:iwstring",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)

cc_library(
    name = "demerit",
    srcs = [
        "demerit.cc",
    ],
    hdrs = [
        "demerit.h",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "extended_pi",
    srcs = [
        "extended_pi.cc",
    ],
    hdrs = [
        "extended_pi.h",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "fileconv_lib",
    hdrs = [
        "fileconv_opts.h",
        "fix_structures.h",
        "known_fragment_data.h",
        "remove_duplicate_fragments.h",
    ],
    srcs = [
        "fileconv_opts.cc",
        "fix_structures.cc",
        "known_fragment_data.cc",
        "remove_duplicate_fragments.cc",
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "iwdemerit_lib",
    srcs = [
        "iwdemerit_lib.cc",
    ],
    hdrs = [
        "iwdemerit_lib.h",
    ],
    deps = [
        ":demerit",
        "//Molecule_Lib:iwmolecule",
    ],
)


cc_library(
    name = "linear_path_fingerprint",
    srcs = [
        "linear_path_fingerprint.cc",
    ],
    hdrs = [
        "linear_path_fingerprint.h",
        "ct_cache.h",
    ],
    deps = [
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "maccskeys_fn5_lib",
    srcs = [
        "maccskeys_fn5.cc",
    ],
    hdrs = [
        "maccskeys_fn5.h"
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
)

cc_library(
    name = "molecule_database_options",
    srcs = [
        "molecule_database_options.cc",
    ],
    hdrs = [
        "molecule_database_options.h",
    ],
    deps = [
        ":molecule_database_options_cc_proto",
        "//Foundational/iwstring:iwstring",
        "//Molecule_Lib:iwmolecule",
        "@leveldb//:leveldb",
    ],
)

cc_library(
    name = "nvrtspsa",
    srcs = [
        "nvrtspsa.cc",
    ],
    hdrs = [
        "nvrtspsa.h"
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "partial_symmetry",
    srcs = [
        "partial_symmetry.cc",
    ],
    hdrs = [
        "partial_symmetry.h"
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "reaction_duplicate",
    srcs = [
        "reaction_duplicate.cc",
    ],
    hdrs = [
        "reaction_duplicate.h"
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "ring_ext_rep",
    srcs = [
        "ring_ext_rep.cc",
    ],
    hdrs = [
        "ring_ext_rep.h"
    ],
    deps = [
        "//Molecule_Lib:iwmolecule",
    ],
)

cc_library(
    name = "topological_torsion_lib",
    srcs = [
        "topological_torsion.cc",
    ],
    hdrs = [
        "topological_torsion.h",
    ],
    deps = [
        "//Foundational/iwstring:iwstring",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Molecule_Lib:iwmolecule"
    ],
    linkopts = [
        "-lz"
    ],
)

cc_library(
    name = "tp_first_pass_lib",
    srcs = [
        "tp_first_pass_lib.cc",
    ],
    hdrs = [
        "tp_first_pass_lib.h",
    ],
    deps = [
        "//Foundational/iwstring:iwstring",
        "//Foundational/iwmisc:iwmisc",
        "//Molecule_Lib:iwmolecule"
    ],
    linkopts = [
        "-lz"
    ],
)

cc_test(
    name = "atom_separation_test",
    srcs = [
        "atom_separation_test.cc",
    ],
    deps = [
        ":atom_separations",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ],
)

cc_test(
    name = "extended_pi_test",
    srcs = [
        "extended_pi_test.cc",
    ],
    deps = [
        ":extended_pi",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
        "//Molecule_Lib:moleculeio",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test( 
    name = "fileconv_test",
    srcs = [
        "fileconv_test.cc",
    ],
    deps = [
        ":fileconv_lib",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
        "@com_google_protobuf//:protobuf",
        "//Molecule_Lib:moleculeio",
    ],
    timeout="short",
)

cc_test(
    name = "iwdemerit_lib_test",
    srcs = [
        "iwdemerit_lib_test.cc",
    ],
    deps = [
        ":demerit",
        ":iwdemerit_lib",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
        "//Molecule_Lib:moleculeio",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "partial_symmetry_test",
    srcs = [
        "partial_symmetry_test.cc",
    ],
    deps = [
        ":partial_symmetry",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
        "//Molecule_Lib:moleculeio",
    ],
    timeout="short",
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "topological_torsion_test",
    srcs = [
        "topological_torsion_test.cc",
    ],
    deps = [
        ":topological_torsion_lib",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
        "//Molecule_Lib:moleculeio",
    ],
    timeout="short",
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "tp_first_pass_test",
    srcs = [
        "tp_first_pass_test.cc",
    ],
    deps = [
        ":tp_first_pass_lib",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
        "//Molecule_Lib:moleculeio",
    ],
    timeout="short",
    linkopts = [
        "-lz",
    ],
)

#cc_test(
#    name = "rxn_to_openrxn_test",
#    srcs = [
#        "rxn_to_openrxn_test.cc",
#    ],
#    deps = [
#        "rxn_to_openrxn_lib",
#        "//Molecule_Lib:iwmolecule",
#        "//Molecule_Lib:iwreaction",
#        "//Molecule_Lib:moleculeio",
#        "@googletest//:gtest_main",
#        "@googletest//:gtest",
#    ],
#    linkopts = [
#    ],
#)
