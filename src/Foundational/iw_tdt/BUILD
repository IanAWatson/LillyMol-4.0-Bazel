cc_library(
    name = "iw_tdt",
    visibility = ["//visibility:public"],
    srcs = [
        "iw_tdt.cc",
        "_iw_tdt_float.cc",
        "_iw_tdt_double.cc",
        "_iw_tdt_int.cc",
        "_iw_tdt_iwstring.cc",
    ],
    hdrs = [
        "iw_tdt.h",
    ],
    deps = [
        "//Foundational/iwmisc:iwconfig",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwaray:iwaray",
        "//Foundational/iwstring:iwstring",
    ],
)

cc_binary (
    name = "tdt_sort",
    srcs = [
        "tdt_sort.cc",
    ],
    deps = [
        ":iw_tdt",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline_v2:cmdline_v2",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwqsort:iwqsort",
    ],
)


cc_binary (
    name = "tdt_stats",
    srcs = [
        "tdt_stats.cc",
    ],
    deps = [
        ":iw_tdt",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
    ],
)

