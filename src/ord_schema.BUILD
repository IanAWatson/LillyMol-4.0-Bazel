proto_library(
    name = "reaction_proto",
    srcs = [
        "proto/reaction.proto"
    ],
)

cc_proto_library(
    name = "reaction_cc_proto",
    deps = [
        ":reaction_proto",
    ],
    visibility = ["//visibility:public"],
)

proto_library(
    name = "dataset_proto",
    srcs = [
        "proto/dataset.proto"
    ],
    deps = [
        ":reaction_proto",
    ],
)

cc_proto_library(
    name = "dataset_cc_proto",
    deps = [
        ":dataset_proto",
    ],
    visibility = ["//visibility:public"],
)

