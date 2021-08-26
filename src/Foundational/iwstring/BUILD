cc_library(
     name = "iwstring",
     visibility = ["//visibility:public"],
     srcs = [
        "append_number.cc",
        "hexprint.cc",
        "iwgetline.cc",
        "iwstring.cc",
        "IWString_class.cc",
        "iwstring_and_file_descriptor.cc",
        "iwwrite.cc",
        "iwwrite_block.cc",
        "iwstrncasecmp.cc",
        "iwstring_strstream.cc",
        "iwzlib.cc",
        "iw_stl_hash_map.cc",
        "remove_suffix.cc",
        "string_change.cc",
        "string_relationals.cc",
        "unhtml.cc",
     ],
     hdrs = [
         "iwhash.h",
         "iwstring.h",
         "iw_stl_hash_set.h",
         "iw_stl_hash_map.h",
         "iwzlib.h",
     ],
     deps = [
         "//Foundational/iwaray:iwaray",
         "//Foundational/iwmisc:iwconfig",
     ],
     linkopts = [
         "-lm",
     ]
)

cc_binary(
    name = "tsclass",
    srcs = [
        "tsclass.cc",
        "should_match.cc",
        "should_match.h",
    ],
    deps = [
        ":iwstring",
        "//Foundational/cmdline:iwcmdline",
    ],
    linkopts = [
        "-lm",
    ],
)

cc_binary(
    name = "tspassbyref",
    srcs = [
        "tspassbyref.cc",
    ],
    deps = [
        ":iwstring",
    ],
)
