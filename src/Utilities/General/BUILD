load("@com_google_protobuf//:protobuf.bzl", "py_proto_library")
load("@rules_proto_grpc//ruby:defs.bzl", "ruby_proto_compile", "ruby_proto_library")

proto_library(
    name = "class_label_translation_proto",
    srcs = [
        "class_label_translation.proto",
    ],
)

cc_proto_library(
    name = "class_label_translation_cc_proto",
    deps = [
        "class_label_translation_proto",
    ],
    visibility = [
        "//visibility:public",
    ],
)

py_proto_library(
    name = "class_label_translation_py_proto",
    srcs = [
        "class_label_translation.proto",
    ]
)

ruby_proto_compile(
    name = "class_label_translation_rb_proto",
    protos = [
        "class_label_translation_proto",
    ],
)

proto_library(
    name = "descriptor_model_proto",
    srcs = [
        "descriptor_model.proto",
    ],
)

cc_proto_library(
    name = "descriptor_model_cc_proto",
    deps = [
        "descriptor_model_proto",
    ],
)

ruby_proto_compile(
    name = "descriptor_model_rb_proto",
    protos = [
        "descriptor_model_proto",
    ],
)

proto_library(
    name = "feature_scaling_proto",
    srcs = [
        "feature_scaling.proto",
    ],
)

cc_proto_library(
    name = "feature_scaling_cc_proto",
    deps = [
        "feature_scaling_proto",
    ],
    visibility = [
        "//visibility:public",
    ],
)

py_proto_library(
    name = "feature_scaling_py_proto",
    srcs = [
        "feature_scaling.proto",
    ]
)

cc_library(
    name = "scaler",
    srcs = [
    ],
    hdrs = [
        "scaler.h",
    ],
    deps = [
        "feature_scaling_cc_proto",
    ],
    visibility = [
        "//visibility:public",
    ],
)

cc_binary(
    name = "average",
    srcs = [
        "average.cc",
        "iwtokeniser.cc",
        "iwtokeniser.h",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwstring:iwstring",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:iwre2",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "class_label_translation",
    srcs = [
        "class_label_translation.cc",
    ],
    deps = [
        ":class_label_translation_cc_proto",
        "//Foundational/cmdline_v2:cmdline_v2",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:proto_support",
        "//Foundational/iwstring:iwstring",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "concat_files",
    srcs = [
        "concat_files.cc",
        "iwtokeniser.cc",
        "iwtokeniser.h",
    ],
    deps = [
        "//Foundational/iwstring:iwstring",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "correlate",
    srcs = [
       "correlate.cc",
    ],
    deps = [
        "//Foundational/iwstring:iwstring",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:iw_tabular_data",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "descriptor_file_to_svml",
    srcs = [
        "descriptor_file_to_svml.cc",
    ],
    deps = [
        "descriptor_model_cc_proto",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:proto_support",
        "//Foundational/iwstring:iwstring",
    ],
)

cc_binary(
    name = "dfilefilter",
    srcs = [
       "dfilefilter.cc",
    ],
    deps = [
        "//Foundational/iwstring:iwstring",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:iw_tabular_data",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "feature_scaling",
    srcs = [
        "feature_scaling.cc",
    ],
    deps = [
        ":scaler",
        ":feature_scaling_cc_proto",
        "//Foundational/cmdline_v2:cmdline_v2",
        "//Foundational/iwmisc:proto_support",
        "//Foundational/iwstring:iwstring",
        "//Foundational/data_source:iwstring_data_source",
    ],
    linkopts = [
        "-lz",
    ],
)

py_binary(
    name = "feature_scaling_create",
    srcs = [
        "feature_scaling_create.py",
    ],
    deps = [
        "feature_scaling_py_proto",
    ],
)


cc_binary(
    name = "fetch_sdf_quick",
    srcs = [
        "fetch_sdf_quick.cc",
    ],
    deps = [
        "//Foundational/iwstring:iwstring",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "fetch_smiles",
    srcs = [
        "fetch_smiles.cc",
        "iwtokeniser.cc",
        "iwtokeniser.h",
    ],
    deps = [
        "//Foundational/iwstring:iwstring",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "fetch_smiles_quick",
    srcs = [
        "fetch_smiles_quick.cc",
        "iwtokeniser.cc",
        "iwtokeniser.h",
    ],
    deps = [
        "//Foundational/iwstring:iwstring",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "iwsplit",
    srcs = [
        "iwsplit.cc",
    ],
    deps = [
        "//Foundational/cmdline_v2:cmdline_v2",
        "//Foundational/iwstring:iwstring",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "iwcut",
    srcs = [
        "iwcut.cc",
        "iwtokeniser.cc",
        "iwtokeniser.h",
    ],
    deps = [
        "//Foundational/iwstring:iwstring",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "normalise",
    srcs = [
        "normalise.cc",
        "iwtokeniser.cc",
        "iwtokeniser.h",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwstring:iwstring",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "notenoughvariance",
    srcs = [
        "notenoughvariance.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwstring:iwstring",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "random_records",
    srcs = [
        "random_records.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwstring:iwstring",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "tcount",
    srcs = [
        "tcount.cc",
        "tokens_in_quoted_string.cc",
        "tokens_in_quoted_string.h",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iwstring:iwstring",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "unique_rows",
    srcs = [
        "unique_rows.cc",
        "iwtokeniser.cc",
        "iwtokeniser.h",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline_v2:cmdline_v2",
        "//Foundational/iwstring:iwstring",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_test(
    name = "scaler_test",
    srcs = [
        "scaler_test.cc",
    ],
    deps = [
        ":feature_scaling_cc_proto",
        ":scaler",
        "@googletest//:gtest_main",
        "@googletest//:gtest",
    ],
)
