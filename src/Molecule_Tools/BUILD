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
    name = "fileconv",
    srcs = [
        "fileconv.cc",
        "fix_structures.cc",
        "fix_structures.h",
        "do_remove_duplicate_fragments.cc",
        "do_remove_duplicate_fragments.h",
        "known_fragment_data.cc",
        "known_fragment_data.h",
    ],
    deps = [
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
    name = "iwdescr",
    srcs = [
        "iwdescr.cc",
        "qry_wcharge.cc",
        "qry_wcharge.h",
    ],
    deps = [
        ":nvrtspsa",
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
    name = "mol2qry",
    srcs = [
        "mol2qry.cc",
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
    name = "tp_first_pass",
    srcs = [
        "tp_first_pass.cc",
    ],
    deps = [
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
