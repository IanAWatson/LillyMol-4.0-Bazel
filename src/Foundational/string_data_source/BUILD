cc_library(
    name = "string_data_source",
    visibility = ["//visibility:public"],
    srcs = [
        "iwstring_string_data_source.cc",
    ],
    hdrs = [
        "iwstring_string_data_source.h",
    ],
    deps = [
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwstring:iwstring",
        "//Foundational/iwmisc:iwre2",
        "@re2//:re2",
    ],
)

cc_test (
    name = "iwstring_string_data_source_test",
    srcs = [
        "iwstring_string_data_source_test.cc",
    ],
    deps = [
        ":string_data_source",
        "//Foundational/iwstring:iwstring",
        "//Foundational/iwmisc:iwre2",
        "@re2//:re2",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ],
)
