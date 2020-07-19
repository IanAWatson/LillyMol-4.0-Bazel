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
        "//Foundational/iwaray:iwaray",
        "@re2//:re2",
        "//Foundational/iwstring:iwstring",
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
    name = "tiwds",
    srcs = [
        "tiwds.cc",
    ],
    deps = [
        ":iwstring_data_source",
        "//Foundational/cmdline:iwcmdline",
    ],
)
