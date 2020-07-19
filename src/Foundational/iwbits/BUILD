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
    ],
    hdrs = [
        "iwbits.h",
        "dy_fingerprint.h",
        "iwbits_support.h",
    ],
    deps = [
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwconfig",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwstring:iwstring",
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
