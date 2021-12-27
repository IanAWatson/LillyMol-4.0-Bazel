
proto_library(
    name = "dataset_proto",
    srcs = [
        "dataset.proto",
    ],
)

cc_proto_library(
    name = "dataset_cc_proto",
    deps = [
        ":dataset_proto",
    ],
)

cc_binary(
    name = "get_3d_smiles",
    srcs = [
        "get_3d_smiles.cc",
    ],
    deps = [
        "dataset_cc_proto",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwtfdata_record",
         "//Molecule_Lib:iwmolecule",
         "//Molecule_Lib:moleculeio",
    ],
)
cc_binary(
    name = "smu",
    srcs = [
        "smu_main.cc",
        "smu.cc",
        "smu.h",
        "smu_results.cc",
        "smu_results.h",
    ],
    deps = [
         "//Foundational/cmdline:iwcmdline",
         "//Foundational/iwaray:iwaray",
         "//Foundational/iwmisc:iwmisc",
         "//Molecule_Lib:iwmolecule",
         "//Molecule_Lib:moleculeio",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "smu_duplicates",
    srcs = [
        "smu_duplicates.cc",
    ],
    deps = [
         "//Foundational/cmdline:iwcmdline",
         "//Foundational/iwaray:iwaray",
         "//Foundational/iwmisc:iwmisc",
         "//Molecule_Lib:iwmolecule",
         "//Molecule_Lib:moleculeio",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "smu_ring_membership_change",
    srcs = [
        "smu_ring_membership_change.cc",
    ],
    deps = [
         "//Foundational/cmdline:iwcmdline",
         "//Foundational/iwaray:iwaray",
         "//Foundational/iwmisc:iwmisc",
         "//Molecule_Lib:iwmolecule",
         "//Molecule_Lib:moleculeio",
    ],
    linkopts = [
        "-lz",
    ],
)



cc_binary(
    name = "smu_rotatable_bond",
    srcs = [
        "smu_rotatable_bond.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
         "//Foundational/cmdline:iwcmdline",
         "//Foundational/iwaray:iwaray",
         "//Foundational/iwmisc:iwmisc",
         "//Molecule_Lib:iwmolecule",
         "//Molecule_Lib:moleculeio",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "smu_superimpose_3d_matches",
    srcs = [
        "smu_superimpose_3d_matches.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
         "//Foundational/cmdline:iwcmdline",
         "//Foundational/iwaray:iwaray",
         "//Foundational/iwmisc:iwmisc",
         "//Molecule_Lib:iwmolecule",
         "//Molecule_Lib:moleculeio",
    ],
    linkopts = [
        "-lz",
    ],
)
