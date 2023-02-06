load("@com_google_protobuf//:protobuf.bzl", "py_proto_library")
load("@rules_proto_grpc//ruby:defs.bzl", "ruby_proto_compile", "ruby_proto_library")

proto_library(
    name = "gfp_proto",
    srcs = [
        "gfp.proto",
    ],
)

cc_proto_library(
    name = "gfp_cc_proto",
    deps = [
        ":gfp_proto",
    ],
)

py_proto_library(
    name = "gfp_py_proto",
    srcs = [
        "gfp.proto",
    ],
)

proto_library(
    name = "gfp_model_proto",
    srcs = [
        "gfp_model.proto",
    ],
)

cc_proto_library(
    name = "gfp_model_cc_proto",
    deps = [
        ":gfp_model_proto",
    ],
)

py_proto_library(
    name = "gfp_model_py_proto",
    srcs = [
        "gfp_model.proto",
    ],
)

ruby_proto_compile(
    name = "gfp_model_rb_proto",
    protos = [
        "gfp_model_proto",
    ],
)

proto_library(
    name = "gfp_to_svm_lite_proto",
    srcs = [
        "gfp_to_svm_lite.proto",
    ],
)

cc_proto_library(
    name = "gfp_to_svm_lite_cc_proto",
    deps = [
        ":gfp_to_svm_lite_proto",
    ],
)

py_proto_library(
    name = "gfp_to_svm_lite_py_proto",
    srcs = [
        "gfp_to_svm_lite.proto",
    ],
)

proto_library(
    name = "iwstats_proto",
    srcs = [
        "iwstats.proto",
    ],
)

cc_proto_library( 
    name = "iwstats_cc_proto",
    deps = [
        "iwstats_proto",
    ],
)

py_proto_library(
    name = "iwstats_py_proto",
    srcs = [
        "iwstats.proto",
    ],
)

proto_library(
    name = "nearneighbors",
    srcs = [
      "nearneighbors.proto",
    ],
)

cc_proto_library(
    name = "nearneighbors_cc_proto",
    deps = [
        "nearneighbors",
    ],
)

py_proto_library(
    name = "nearneighbors_py_proto",
    srcs = [
        "nearneighbors.proto",
    ],
)

ruby_proto_compile(
    name = "iwstats_rb_proto",
    protos = [
        "iwstats_proto",
    ],
)

cc_binary(
  name = "properties_benchmark",
  srcs = [
      "properties_benchmark.cc",
  ],
  deps = [
      ":gfp",
      ":tanimoto_float",
        "@benchmark//:benchmark",
  ],
)

cc_binary(
    name = "descriptor_file_to_01_fingerprints",
    srcs = [
        "descriptor_file_to_01_fingerprints.cc",
    ],
    deps = [
        "gfp",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "gfp_catboost_score",
    srcs = [
        "gfp_catboost_score.cc",
        "Catboost.model.bin.cpp",
    ],
    deps = [
        ":gfp",
        ":gfp_to_svm_lite_cc_proto",
        ":gfp_model_cc_proto",
        "//Foundational/cmdline_v2:cmdline_v2",
        "//Foundational/data_source:iwstring_data_source",
        "//Utilities/General:class_label_translation_cc_proto",
        "//Utilities/General:scaler",
    ],
    tags = [
        "catboost",
    ],
)

cc_binary(
    name = "gfp_leader",
    srcs = [
        "gfp_leader_v2.cc",
        "leader.h",
    ],
    deps = [
        ":gfp",
        "//Foundational/accumulator:accumulator",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "gfp_leader_standard",
    srcs = [
        "gfp_leader_standard.cc",
    ],
    copts = [
#       "-fopenmp=libiomp5",  llvm
        "-fopenmp",
    ],
    deps = [
        ":gfp",
        ":gfp_standard",
        ":nearneighbors_cc_proto",
        "//Foundational/data_source:iwtfdata_record",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
    ],
    linkopts = [
        "-lz",
#       "-L/usr/lib/llvm-11/lib/",  # clang
        "-lgomp",
    ],
)

cc_binary(
    name = "gfp_naive_bayesian",
    srcs = [
      "gfp_naive_bayesian.cc",
      "bit_subset_v2.h",
    ],
    deps = [
        ":gfp",
        ":spread_v2",
        ":smiles_id_dist",
        "//Foundational/accumulator:accumulator",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "gfp_spread",
    srcs = [
        "gfp_spread_v2.cc",
    ],
    deps = [
        ":gfp",
        ":spread_v2",
        ":smiles_id_dist",
        "//Foundational/accumulator:accumulator",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "gfp_spread_buckets_v2",
    srcs = [
        "gfp_spread_buckets_v2.cc",
    ],
    deps = [
        ":gfp",
        ":spread_v2",
        ":smiles_id_dist",
        "//Foundational/accumulator:accumulator",
        "//Foundational/iwmisc:iwmisc",
    ],
)

cc_binary(
    name = "gfp_spread_standard",
    srcs = [
        "gfp_spread_standard.cc",
    ],
    copts = [
#       "-fopenmp=libiomp5", llvm
        "-fopenmp",
    ],
    deps = [
        ":gfp",
        ":gfp_standard",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
    ],
    linkopts = [
        "-lz",
#       "-L/usr/lib/llvm-11/lib/",  # clang
        "-lgomp",
    ],
)

cc_binary(
    name = "gfp_spread_standard_mpi",
    srcs = [
        "gfp_spread_standard_mpi.cc",
    ],
    deps = [
        ":gfp",
        ":gfp_standard",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iw_tdt:iw_tdt",
    ],
    copts = [
        "-I /usr/include/x86_64-linux-gnu/mpich",
    ],
    tags = [
        "mpi",
    ],
    linkopts = [
        "-lfabric",
        "-lmpich",
        "-lhwloc",
    ],
)

cc_binary(
    name = "gfp_nearneighbours_single_file",
    srcs = [
        "gfp_nearneighbours_single_file.cc",
        "fp_and_smiles.cc",
    ],
    deps = [
        ":gfp",
        "//Foundational/histogram:iwhistogram",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "gfp_lnearneighbours",
    srcs = [
        "gfp_lnearneighbours.cc",
        "fp_and_smiles.cc",
    ],
    deps = [
        ":gfp",
        ":neighbour_list",
    ],
    linkopts = [
        "-lz",
    ],
)

# Runs faster single threaded with no omp options specified.
cc_binary(
    name = "gfp_lnearneighbours_standard",
    srcs = [
        "gfp_lnearneighbours_standard.cc",
        "fp_and_smiles.cc",
    ],
    copts = [
#       "-fopenmp=libiomp5", llvm
        "-fopenmp",
    ],
    deps = [
        ":gfp",
        ":gfp_standard",
        ":neighbour_list",
    ],
    linkopts = [
#       "-L/usr/lib/llvm-11/lib/",  # clang
        "-lgomp",
    ],
)

cc_binary(
    name = "gfp_pairwise_distances",
    srcs = [
        "gfp_pairwise_distances.cc",
    ],
    deps = [
        "//Foundational/histogram:iwhistogram",
        ":gfp",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "gfp_profile_activity_by_bits",
    srcs = [
        "gfp_profile_activity_by_bits.cc",
    ],
    deps = [
        ":gfp",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "gfp_to_descriptors",
    srcs = [
        "gfp_to_descriptors.cc",
    ],
    deps = [
        ":gfp",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "gfp_to_proto",
    srcs = [
        "gfp_to_proto.cc",
    ],
    deps = [
        ":gfp",
        ":gfp_cc_proto",
        "//Foundational/data_source:iwtfdata_record",
        "//Foundational/iwmisc:iwmisc",
        "@snappy//:snappy",
    ],
)

cc_binary(
    name = "gfp_to_svm_lite",
    srcs = [
        "gfp_to_svm_lite.cc",
    ],
    deps = [
        ":gfp",
        ":gfp_to_svm_lite_cc_proto",
        "//Foundational/iwmisc:activity_from_file",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "nplotnn",
    srcs = [
        "nplotnn.cc",
    ],
    deps = [
        ":distance_scaling",
        ":smiles_id_dist",
        "//Foundational/histogram:iwhistogram",
        "//Foundational/iwmisc:iwmisc",
        "//Molecule_Lib:iwmolecule",
        "//Molecule_Lib:moleculeio",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "gfp_distance_filter_standard",
    srcs = [
        "gfp_distance_filter_standard.cc",
    ],
    deps = [
        ":gfp_standard",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt:iw_tdt",
    ]
)

cc_binary(
    name = "gfp_distance_matrix",
    srcs = [
        "gfp_distance_matrix.cc",
    ],
    deps = [
        ":build_pool",
        ":gfp",
        "//Foundational/histogram:iwhistogram",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "gfp_svmfp_evaluate",
    srcs = [
        "gfp_svmfp_evaluate.cc",
    ],
    deps = [
        ":gfp",
        ":gfp_to_svm_lite_cc_proto",
        ":gfp_model_cc_proto",
        "//Foundational/cmdline_v2:cmdline_v2",
        "//Foundational/data_source:iwstring_data_source",
        "//Utilities/General:class_label_translation_cc_proto",
        "//Utilities/General:scaler",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "test_sparse_bitvector_performance",
    srcs = [
        "test_sparse_bitvector_performance.cc",
    ],
    deps = [
        ":build_pool",
        ":gfp",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/iw_tdt:iw_tdt",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
    ],
    linkopts = [
        "-lz",
    ],
)


cc_binary(
    name = "iwstats",
    srcs = [
        "iwstats.cc",
    ],
    deps = [
        ":bsquared",
        ":metric",
        ":iwstats_cc_proto",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:proto_support",
        "//Foundational/iwqsort:iwqsort",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_binary(
    name = "tdt_sort",
    srcs = [
        "tdt_sort.cc",
    ],
    deps = [
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline_v2:cmdline_v2",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwmisc:iwre2",
        "//Foundational/iwqsort:iwqsort",
        "//Foundational/iw_tdt:iw_tdt",
    ],
    linkopts = [
        "-lz",
    ],
)

cc_library(
    name = "spread_v2",
    srcs = [
        "spread_v2.cc",
    ],
    hdrs = [
        "spread_v2.h",
    ],
    deps = [
        ":gfp",
        ":smiles_id_dist",
    ],
)

cc_library(
    name = "build_pool",
    srcs = [
        "build_pool.cc",
    ],
    hdrs = [
    ],
    deps = [
        ":gfp",
        ":smiles_id_dist",
    ],
)

cc_library(
    name = "smiles_id_dist",
    srcs = [
        "smiles_id_dist.cc",
    ],
    hdrs = [
        "smiles_id_dist.h",
    ],
    deps = [
        ":gfp",
    ],
)


cc_library(
    name = "gfp",
    srcs = [
        "bic_table.cc",
        "bit_subset.cc",
        "gfp.cc",
        "dyfp.cc",
        "fixed_size_counted_fingerprint_uchar.cc",
        "_fixed_size_counted_fingerprint_uchar.cc",
        "_fixed_size_counted_fingerprint_uint.cc",
        "multi_conformer.cc",
        "sparsefp.cc",
        "sparse_collection.cc",
        "tversky.cc",
        "various_distance_metrics.cc",
    ],
    hdrs = [
        "bit_subset.h",
        "gfp.h",
        "dyfp.h",
        "fixed_size_counted_fingerprint.h",
        "multi_conformer.h",
        "sparsefp.h",
        "sparse_collection.h",
        "tversky.h",
        "various_distance_metrics.h",
    ],
    deps = [
        "gfp_to_svm_lite_cc_proto",
        ":tanimoto_float",
        "//Foundational/accumulator:accumulator",
        "//Foundational/cmdline:iwcmdline",
        "//Foundational/data_source:iwstring_data_source",
        "//Foundational/iwbits:iwbits",
        "//Foundational/iwmisc:iwmisc",
        "//Foundational/iwmisc:proto_support",
        "//Foundational/iwmisc:sparse_fp_creator",
        "//Foundational/iwqsort:iwqsort",
        "//Foundational/iw_tdt:iw_tdt",
        "@re2//:re2",
    ],
)

cc_library(
    name = "gfp_standard",
    visibility = [
        "//visibility:public",
    ],
    srcs = [
        "gfp_standard.cc",
    ],
    hdrs = [
        "gfp_standard.h",
    ],
    deps = [
        ":gfp",
        ":smiles_id_dist",
        "//Foundational/iwmisc:iwmisc",
    ],
)

cc_library(
    name = "neighbour_list",
    srcs = [
        "neighbour_list.cc",
        "_neighbour_list_fp_and_smiles.cc",
    ],
    hdrs = [
        "neighbour_list.h",
    ],
    deps = [
        ":smiles_id_dist",
    ],
)

cc_library(
    name = "sparse_collection",
    visibility = [
        "//visibility:public",
    ],
    srcs = [
        "sparse_collection.cc",
    ],
    hdrs = [
        "sparse_collection.h",
    ],
    deps = [
        ":gfp",
        "//Foundational/cmdline:iwcmdline",
    ],
)

cc_library(
    name = "distance_scaling",
    srcs = [
        "distance_scaling.cc",
    ],
    hdrs = [
        "distance_scaling.h",
    ],
    deps = [
        "//Foundational/iwmisc:iwmisc",
    ],
)

cc_library(
    name = "bsquared",
    srcs = [
        "bsquared.cc"
    ],
    hdrs = [
        "bsquared.h"
    ],
    deps = [
        "//Foundational/accumulator:accumulator"
    ],
)

cc_library(
    name = "metric",
    srcs = [
        "Enrichment.cc",
        "Metric.cc"
    ],
    hdrs = [
        "Enrichment.h",
        "Metric.h"
    ],
    deps = [
        "//Foundational/accumulator:accumulator"
    ],
)

cc_library(
    name = "tanimoto_float",
    srcs = [
        "tanimoto_float.cc",
    ],
    hdrs = [
        "tanimoto_float.h",
    ],
    deps = [
    ],
    copts = [
        "-mavx512vl -mavx -mssse3"
    ]
)

cc_binary(
    name = "sparsefp_benchmark",
    srcs = [
        "sparsefp_benchmark.cc",
    ],
    deps = [
        "@benchmark//:benchmark",
        ":gfp",
        "//Foundational/iwmisc:sparse_fp_creator",
    ],
)


cc_test(
    name = "bit_subset_test",
    srcs = [
        "bit_subset_test.cc",
    ],
    timeout="short",
    deps = [
        ":gfp",
        "@googletest//:gtest_main",
    ],
)

cc_test (
    name = "tanimoto_float_test",
    srcs = [
        "tanimoto_float_test.cc",
    ],
    deps = [
        ":tanimoto_float",
        "@googletest//:gtest_main",
    ],
)

cc_test(
    name = "sparsefp_test",
    srcs = [
        "sparsefp_test.cc",
    ],
    deps = [
        ":gfp",
        "@googletest//:gtest_main",
    ],
    timeout="short",
)
