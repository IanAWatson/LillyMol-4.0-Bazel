cc_library(
    name = "activity_from_file",
    visibility = ["//visibility:public"],
    srcs = [
    ],
    hdrs = [
         "activity_data_from_file.h",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwstring:iwstring",
    ],
)

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
        "fetch_via_regexp.cc",
        "fraction_as_string.cc",
        "int_comparator.cc",
        "iwdigits.cc",
        "iwdirname.cc",
        "logical_expression.cc",
        "msi_object.cc",
        "new_double.cc",
        "new_float.cc",
        "new_int.cc",
        "new_unsigned_int.cc",
        "report_progress.cc",
        "restore_state.cc",
        "write_space_suppressed_string.cc",
        "_numeric_data_from_file_int.cc",
        "_numeric_data_from_file_float.cc",
        "_set_or_unset_float.cc",
    ],
    hdrs = [
        "combinations.h",
        "fetch_via_regexp.h",
        "iwarchive.h",
        "iwdigits.h",
        "iwminmax.h",
        "iw_time.h",
        "logical_expression.h",
        "matcher.h",
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

cc_library(
    name = "memoized_floats",
    visibility = ["//visibility:public"],
    srcs = [
        "memoized_floats.cc",
    ],
    hdrs = [
        "memoized_floats.h",
    ],
    deps = [
        ":iwmisc",
    ],
    linkopts = [
    ],
)

cc_library(
    name = "proto_support",
    visibility = ["//visibility:public"],
    srcs = [
        "proto_support.cc",
    ],
    hdrs = [
        "proto_support.h",
    ],
    deps = [
	"@com_google_protobuf//:protobuf",
        "@com_google_absl//absl/strings:strings",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwstring:iwstring",
    ],
)

cc_binary(
    name = "test_iwdigits",
    srcs = [
        "test_iwdigits.cc",
    ],
    deps = [
        ":iwmisc",
        "//Foundational/cmdline_v2:cmdline_v2",
    ],
    linkopts = [
    ],
)

cc_test(
    name = "combinations_test",
    srcs = [
        "combinations_test.cc",
    ],
    deps = [
        ":iwmisc",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "fraction_as_string_test",
    srcs = [
        "fraction_as_string_test.cc",
    ],
    deps = [
        ":iwmisc",
        "@googletest//:gtest_main",
    ]
)

cc_test  (
    name = "iwdigits_test",
    srcs = [
        "iwdigits_test.cc",
    ],
    deps = [
        ":iwmisc",
	"@com_google_protobuf//:protobuf",
        "@googletest//:gtest_main",
    ],
    linkopts = [
        "-lm",
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
        "-lm",
    ],
)

cc_test(
    name = "logical_expression_test",
    srcs = [
        "logical_expression_test.cc",
    ],
    deps = [
        ":iwmisc",
        "//Foundational/iwstring:iwstring",
        "@googletest//:gtest_main",
    ]
)

cc_test(
    name = "matcher_test",
    srcs = [
        "matcher_test.cc",
    ],
    deps = [
        ":iwmisc",
        "//Foundational/iwstring:iwstring",
        "@googletest//:gtest_main",
    ]
)

cc_test  (
    name = "memoized_floats_test",
    srcs = [
        "memoized_floats_test.cc",
    ],
    deps = [
        ":memoized_floats",
        "@benchmark//:benchmark",
	"@com_google_protobuf//:protobuf",
        "@googletest//:gtest_main",
    ],
    linkopts = [
        "-lm",
    ],
)

cc_test  (
    name = "memoized_floats_benchmark",
    srcs = [
        "memoized_floats_benchmark.cc",
    ],
    deps = [
        ":memoized_floats",
        "@benchmark//:benchmark_main",
    ],
    linkopts = [
    ],
)
