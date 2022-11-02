proto_library(
    name = "dicer_fragment_proto",
    srcs = [
        "dicer_fragment.proto",
    ],
)
    
cc_proto_library(
    name = "dicer_fragment_cc_proto",
    deps = [
        ":dicer_fragment_proto",
    ],
)

proto_library(
    name = "replacement_ring_proto",
    srcs = [
        "replacement_ring.proto",
    ],
)
    
cc_proto_library(
    name = "replacement_ring_cc_proto",
    deps = [
        ":replacement_ring_proto",
    ],
)

cc_binary(
    name = "ring_extraction",
    srcs = [
        "ring_extraction.cc"
    ],
    deps = [
        ":replacement_ring_cc_proto",
        "//Molecule_Tools:reaction_duplicate",
        "//Molecule_Tools:ring_ext_rep",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "ring_replacement",
    srcs = [
        "ring_replacement.cc"
    ],
    deps = [
        ":replacement_ring_cc_proto",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwmisc:iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "//Molecule_Lib:substructure_cc_proto",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_binary(
    name = "get_substituent",
    srcs = [
        "get_substituent.cc",
    ],
    deps = [
        ":dicer_fragment_cc_proto",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
        "@com_google_protobuf//:protobuf",
    ],
)
