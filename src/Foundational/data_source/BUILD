package(default_visibility = ["//visibility:public"])
cc_library(
    name = "binary_data_file",
    visibility = ["//visibility:public"],
    srcs = [
        "binary_data_file.cc",
    ],
    hdrs = [
        "binary_data_file.h",
    ],
    deps = [
        "//Foundational/iwstring:iwstring",
    ],

)
cc_library(
    name = "iwstring_data_source",
    visibility = ["//visibility:public"],
    srcs = [
        "iwmmap.cc",
        "iwstring_data_source.cc",
        "string_data_source.cc",
    ],
    hdrs = [
        "iwmmap.h",
        "iwstring_data_source.h",
        "string_data_source.h",
    ],
    deps = [
        "//Foundational/iwmisc:iwconfig",
        "//Foundational/iwmisc:iwre2",
        "//Foundational/iwaray:iwaray",
        "@re2//:re2",
        "//Foundational/iwstring:iwstring",
    ],
)

cc_library(
    name = "iwtfdata_record",
    visibility = ["//visibility:public"],
    srcs = [
        "tfdatarecord.cc",
    ],
    hdrs = [
        "tfdatarecord.h",
    ],
    deps = [
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwstring:iwstring",
        "@crc32//:crc32",
    ],
)

proto_library(
    name = "proto_for_testing",
    srcs = [
        "proto_for_testing.proto",
    ],
)

cc_proto_library(
    name = "proto_for_testing_cc_proto",
    deps = [
        ":proto_for_testing",
    ],
)

cc_test(
    name = "binary_data_file_test",
    timeout="short",
    srcs = [
       "binary_data_file_test.cc",
    ],
    deps = [
        ":binary_data_file",
        ":proto_for_testing_cc_proto",
        "@googletest//:gtest_main",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "iwstring_data_source_test",
    srcs = [
       "iwstring_data_source_test.cc",
    ],
    deps = [
      ":iwstring_data_source",
        "@googletest//:gtest_main",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "test_tfdata_record",
    srcs = [
        "test_tfdata_record.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        ":iwtfdata_record",
    ],
)

cc_binary(
    name = "tiwds",
    srcs = [
        "tiwds.cc",
    ],
    deps = [
        ":iwstring_data_source",
        "//Foundational/cmdline:iwcmdline",
    ],
    linkopts = [
        "-lz",
    ],
)
