cc_library(
    name = "iwbits",
    visibility = ["//visibility:public"],
    srcs = [
        "bits_in_common.c",
        "du_bin2ascii.cc",
        "fastprint.cc",
        "fixed_size_counted_fingerprint.cc",
        "hexbits.cc",
        "iwbits.cc",
        "swap_byte_order.cc",
        "fixed_bit_vector.cc",
    ],
    hdrs = [
        "iwbits.h",
        "dy_fingerprint.h",
        "iwbits_support.h",
        "fixed_bit_vector.h",
    ],
    deps = [
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwconfig",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwstring:iwstring",
    ],
)

cc_binary(
    name = "fixed_bit_vector_benchmark",
    srcs = [
        "fixed_bit_vector_benchmark.cc",
    ],
    deps = [
        "@benchmark//:benchmark",
        ":iwbits",
    ],
)

# cc_binary(
#     name = "fppop",
#     srcs = [
#         "fppop.cc",
#     ],
#     deps = [
#         ":iwbits",
#         "//Foundational/cmdline:iwcmdline",
#         "//Foundational/data_source:iwstring_data_source",
#         "//Foundational/iwmisc:iwconfig",
#         "//Foundational/iwmisc:iwmisc",
#         "//Foundational/iwstring:iwstring",
#         "//Molecule_Lib:iwmolecule",
#     ]
# )

cc_test(
    name = "fixed_bit_vector_test",
    srcs = [
        "fixed_bit_vector_test.cc",
    ],
    deps = [
        ":iwbits",
	"@com_google_protobuf//:protobuf",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ],
    linkopts = [
        "-lz",
    ],
)
