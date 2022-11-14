proto_library(
    name = "bond_length_distribution_proto",
    srcs = [
        "bond_length_distribution.proto",
    ],
    deps = [
    ],
)

cc_proto_library(
    name = "bond_length_distribution_cc_proto",
    deps = [
        ":bond_length_distribution_proto",
    ],
)

cc_binary(
    name = "topology_from_geometry",
    srcs = [
        "topology_from_geometry_main.cc",
    ],
    deps = [
        ":topology_from_geometry_lib",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:proto_support",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ]
)

cc_library(
    name = "empirical_bond_length_distribution",
    srcs = [
        "empirical_bond_length_distribution.cc",
        "empirical_bond_length_distribution.h",
    ],
    deps = [
        ":bond_length_distribution_cc_proto",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:iwmisc",
    ],
)

cc_library(
    name = "topology_from_geometry_lib",
    srcs = [
         "topology_from_geometry.cc",
    ],
    hdrs = [
         "topology_from_geometry.h",
    ],
    deps = [
        ":bond_length_distribution_cc_proto",
        ":empirical_bond_length_distribution",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:proto_support",
        "//Molecule_Lib:iwmolecule"
    ],
    visibility = [
        "//visibility:public",
    ],
)

cc_test(
    name = "empirical_bond_length_distribution_test",

    srcs = [
        "empirical_bond_length_distribution_test.cc",
    ],
    deps = [
        ":bond_length_distribution_cc_proto",
        ":empirical_bond_length_distribution",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
        "@com_google_protobuf//:protobuf",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
)
