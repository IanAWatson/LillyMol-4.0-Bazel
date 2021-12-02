project(gfp)
find_package(Protobuf REQUIRED)
include_directories(SYSTEM ${Protobuf_INCLUDE_DIRS})

set(proto_src
  ${CMAKE_CURRENT_SOURCE_DIR}/gfp_to_svm_lite.proto
)
protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS ${proto_src})

set(gfp_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/bic_table.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/bit_subset.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/gfp.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/dyfp.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/fixed_size_counted_fingerprint_uchar.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_fixed_size_counted_fingerprint_uchar.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_fixed_size_counted_fingerprint_uint.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/multi_conformer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/sparsefp.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/sparse_collection.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/tversky.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/various_distance_metrics.cc
  ${PROTO_SRCS}
  ${PROTO_HDRS}
)

add_library(gfp ${gfp_srcs})
target_include_directories(gfp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(smiles_id_dist_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/smiles_id_dist.cc
)
add_library(smiles_id_dist ${smiles_id_dist_srcs})
target_include_directories(smiles_id_dist PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(spread_v2_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/spread_v2.cc
)
add_library(spread_v2 ${spread_v2_srcs})
target_include_directories(spread_v2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(gfp_standard_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/gfp_standard.cc
)
add_library(gfp_standard ${gfp_standard_srcs})
target_include_directories(gfp_standard PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(neighbour_list_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/neighbour_list.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_neighbour_list_fp_and_smiles.cc
)
add_library(neighbour_list ${neighbour_list_srcs})
target_include_directories(neighbour_list PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(metric_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/Enrichment.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/Metric.cc
)
add_library(metric ${metric_srcs})
target_include_directories(metric PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(bsquared_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/bsquared.cc
)
add_library(bsquared ${bsquared_srcs})
target_include_directories(bsquared PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(common_libs
  accumulator
  cmdline
  iwstring_data_source
  iwre2
  gfp
  sparse_fp_creator
  iwbits
  iw_tdt
  iwstring
  iwmisc
  ${IW_EXTRA_LIBS}
)

set(gfp_spread_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/gfp_spread_v2.cc
)
set(gfp_spread_libs
  spread_v2
  smiles_id_dist
  ${common_libs}
  iwrandom
)
add_executable(gfp_spread ${gfp_spread_srcs})
target_link_libraries(gfp_spread LINK_PUBLIC ${gfp_spread_libs})


set(gfp_leader_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/gfp_leader_v2.cc
)
set(gfp_leader_libs
  ${common_libs}
  smiles_id_dist
)
add_executable(gfp_leader ${gfp_leader_srcs})
target_link_libraries(gfp_leader LINK_PUBLIC ${gfp_leader_libs})


set(gfp_leader_standard_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/gfp_leader_standard.cc
)
set(gfp_leader_standard_libs
  ${common_libs}
  smiles_id_dist
  gfp_standard
)
add_executable(gfp_leader_standard ${gfp_leader_standard_srcs})
target_link_libraries(gfp_leader_standard LINK_PUBLIC ${gfp_leader_standard_libs} -fopenmp)


set(gfp_naive_bayesian_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/gfp_naive_bayesian.cc
)
set(gfp_naive_bayesian_libs
  ${common_libs}
)
add_executable(gfp_naive_bayesian ${gfp_naive_bayesian_srcs})
target_link_libraries(gfp_naive_bayesian LINK_PUBLIC ${gfp_naive_bayesian_libs} -fopenmp)


set(gfp_spread_standard_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/gfp_spread_standard.cc
)
set(gfp_spread_standard_libs
  ${common_libs}
  gfp_standard
  iwrandom
)
add_executable(gfp_spread_standard ${gfp_spread_standard_srcs})
target_link_libraries(gfp_spread_standard LINK_PUBLIC ${gfp_spread_standard_libs} -fopenmp)


set(gfp_nearneighbours_single_file_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/gfp_nearneighbours_single_file.cc
)
set(gfp_nearneighbours_single_file_libs
  ${common_libs}
  gfp_standard
  iwrandom
  iwhistogram
)
add_executable(gfp_nearneighbours_single_file ${gfp_nearneighbours_single_file_srcs})
target_link_libraries(gfp_nearneighbours_single_file LINK_PUBLIC ${gfp_nearneighbours_single_file_libs} -fopenmp)


set(gfp_lnearneighbours_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/gfp_lnearneighbours.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/fp_and_smiles.cc
)
set(gfp_lnearneighbours_libs
  smiles_id_dist
  neighbour_list
  smiles_id_dist
  ${common_libs}
)
add_executable(gfp_lnearneighbours ${gfp_lnearneighbours_srcs})
target_link_libraries(gfp_lnearneighbours LINK_PUBLIC ${gfp_lnearneighbours_libs})


set(iwstats_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/iwstats.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/fp_and_smiles.cc
)
set(iwstats_libs
  metric
  bsquared
  ${common_libs}
  iwrandom
)
add_executable(iwstats ${iwstats_srcs})
target_link_libraries(iwstats LINK_PUBLIC ${iwstats_libs})
