cc_binary(
    name = "leveldb_list",
    srcs = [
        "leveldb_list.cc",
    ],
    deps = [
        "//Foundational/cmdline:iwcmdline",
        "@leveldb//:leveldb",
    ],
    linkopts = [
        "-lz",
        "-lpthread",
        "-lre2",
    ],
)
