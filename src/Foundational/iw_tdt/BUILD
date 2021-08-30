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
        "//Foundational/mtrand:iwrandom",
    ],
)
