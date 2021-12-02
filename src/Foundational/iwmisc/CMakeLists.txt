project(iwmisc)
set(source_files
  ${CMAKE_CURRENT_SOURCE_DIR}/append_first_token_of_name.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/dash_d.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/dash_f.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/dash_s.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/dash_x.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/endian.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/fetch_via_regexp.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/fraction_as_string.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/int_comparator.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwdigits.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwdirname.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/logical_expression.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/msi_object.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/new_double.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/new_float.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/new_int.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/proto_support.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/report_progress.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/write_space_suppressed_string.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_numeric_data_from_file_int.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_numeric_data_from_file_float.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_set_or_unset_float.cc
)

set(header_files
  ${CMAKE_CURRENT_SOURCE_DIR}/fetch_via_regexp.h
  ${CMAKE_CURRENT_SOURCE_DIR}/iwarchive.h
  ${CMAKE_CURRENT_SOURCE_DIR}/iwdigits.h
  ${CMAKE_CURRENT_SOURCE_DIR}/iwminmax.h
  ${CMAKE_CURRENT_SOURCE_DIR}/iw_time.h
  ${CMAKE_CURRENT_SOURCE_DIR}/logical_expression.h
  ${CMAKE_CURRENT_SOURCE_DIR}/minmaxspc.h
  ${CMAKE_CURRENT_SOURCE_DIR}/misc.h
  ${CMAKE_CURRENT_SOURCE_DIR}/msi_object.h
  ${CMAKE_CURRENT_SOURCE_DIR}/numeric_data_from_file.h
  ${CMAKE_CURRENT_SOURCE_DIR}/primes.h
  ${CMAKE_CURRENT_SOURCE_DIR}/proto_support.h
  ${CMAKE_CURRENT_SOURCE_DIR}/report_progress.h
  ${CMAKE_CURRENT_SOURCE_DIR}/set_or_unset.h
  ${CMAKE_CURRENT_SOURCE_DIR}/timsort.hpp
)

set(iwconfig_header ${CMAKE_CURRENT_SOURCE_DIR}/iwconfig.h)

add_library(iwconfig INTERFACE)
target_sources(iwconfig INTERFACE "$<BUILD_INTERFACE:${iwconfig_header}>")
target_include_directories(iwconfig INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

add_library(iwmisc ${source_files} )

set(sparse_fp_creator_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/sparse_fp_creator.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/sparse_fp_creator_support.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/md5c.c
)

add_library(sparse_fp_creator ${sparse_fp_creator_srcs})

set(iwre2_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/iwre2.cc
)
add_library(iwre2 ${iwre2_srcs})

set(tabular_data_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/iw_tabular_data.cc
)
add_library(tabular_data ${tabular_data_srcs})
