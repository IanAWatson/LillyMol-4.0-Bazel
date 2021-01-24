cc_library(
    name = "iwmisc",
    visibility = ["//visibility:public"],
    srcs = [
        "append_first_token_of_name.cc",
        "dash_d.cc",
        "dash_f.cc",
        "dash_s.cc",
        "dash_x.cc",
        "endian.cc",
        "fraction_as_string.cc",
        "int_comparator.cc",
        "iwdigits.cc",
        "logical_expression.cc",
        "msi_object.cc",
        "new_int.cc",
        "new_float.cc",
        "report_progress.cc",
        "write_space_suppressed_string.cc",
        "_numeric_data_from_file_int.cc",
        "_numeric_data_from_file_float.cc",
        "_set_or_unset_float.cc",
    ],
    hdrs = [
        "iwarchive.h",
        "iwdigits.h",
        "iwminmax.h",
        "iw_time.h",
        "logical_expression.h",
        "minmaxspc.h",
        "misc.h",
        "msi_object.h",
        "numeric_data_from_file.h",
        "primes.h",
        "report_progress.h",
        "set_or_unset.h",
        "timsort.hpp",
    ],
    deps = [
        ":internal",
        ":iwconfig",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
    ],
)

cc_library(
    name = "internal",
    hdrs = [
        "new_array_.h",
    ],
)

cc_library(
    name = "iwconfig",
    visibility = ["//visibility:public"],
    hdrs = [
        "iwconfig.h",
    ],
)

cc_library(
    name = "sparse_fp_creator",
    visibility = ["//visibility:public"],
    srcs = [
        "sparse_fp_creator.cc",
        "sparse_fp_creator_support.cc",
        "md5c.c",
    ],
    hdrs = [
        "md5.h",
        "md5_global.h",
        "sparse_fp_creator.h",
    ],
    deps = [
        ":iwmisc",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwbits:iwbits",
        "//Foundational/iwqsort:iwqsort",
        "//Foundational/iwstring:iwstring",
    ],
)

cc_library(
    name = "iw_tabular_data",
    visibility = ["//visibility:public"],
    srcs = [
        "iw_tabular_data.cc",
    ],
    hdrs = [
        "iw_tabular_data.h",
    ],
    deps = [
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwbits:iwbits",
    ],
)

cc_library(
    name = "iwre2",
    visibility = ["//visibility:public"],
    srcs = [
        "iwre2.cc",
    ],
    hdrs = [
        "iwre2.h",
    ],
    deps = [
      "@re2//:re2",
      "//Foundational/iwstring:iwstring",
    ],
)

cc_test  (
    name = "iwre2_test",
    srcs = [
        "iwre2_test.cc",
    ],
    deps = [
        ":iwre2",
        "@googletest//:gtest_main",
    ],
    linkopts = [
        "-lz",
    ],
)
