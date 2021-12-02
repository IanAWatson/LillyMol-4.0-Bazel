project(utilities)
find_package(Protobuf REQUIRED)

set(tcount_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/tcount.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/tokens_in_quoted_string.cc
)
set(tcount_libs
  cmdline
  iwstring_data_source
  iwre2
  iwstring
  ${IW_EXTRA_LIBS}
)

add_executable(tcount ${tcount_srcs})
target_link_libraries(tcount LINK_PUBLIC ${tcount_libs})

set(average_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/average.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwtokeniser.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwtokeniser.h
)
set(average_libs
  accumulator
  cmdline
  iwstring_data_source
  iwstring
  iwre2
  ${IW_EXTRA_LIBS}
)

add_executable(average ${average_srcs})
target_link_libraries(average LINK_PUBLIC ${average_libs})

set(concat_files_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/concat_files.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwtokeniser.cc
)
set(concat_files_libs
  accumulator
  cmdline
  iwmisc
  iwstring_data_source
  iwstring
  iwre2
  ${IW_EXTRA_LIBS}
)
add_executable(concat_files ${concat_files_srcs})
target_link_libraries(concat_files LINK_PUBLIC ${concat_files_libs})


set(correlate_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/correlate.cc
)
set(correlate_libs
  accumulator
  cmdline
  tabular_data
  iwbits
  iwstring_data_source
  iwstring
  iwre2
  iwmisc
  ${IW_EXTRA_LIBS}
)
add_executable(correlate ${correlate_srcs})
target_link_libraries(correlate LINK_PUBLIC ${correlate_libs})


set(dfilefilter_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/dfilefilter.cc
)
set(dfilefilter_libs
  accumulator
  cmdline
  iwstring_data_source
  iwstring
  iwmisc
  iwre2
  ${IW_EXTRA_LIBS}
)
add_executable(dfilefilter ${dfilefilter_srcs})
target_link_libraries(dfilefilter LINK_PUBLIC ${dfilefilter_libs})


set(fetch_sdf_quick_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/fetch_sdf_quick.cc
)
set(fetch_sdf_quick_libs
  accumulator
  cmdline
  iwstring_data_source
  iwstring
  iwmisc
  iwre2
  ${IW_EXTRA_LIBS}
)
add_executable(fetch_sdf_quick ${fetch_sdf_quick_srcs})
target_link_libraries(fetch_sdf_quick LINK_PUBLIC ${fetch_sdf_quick_libs})


set(fetch_smiles_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/fetch_smiles.cc
)
set(fetch_smiles_libs
  accumulator
  cmdline
  iwstring_data_source
  iwstring
  iwre2
  ${IW_EXTRA_LIBS}
)
add_executable(fetch_smiles ${fetch_smiles_srcs})
target_link_libraries(fetch_smiles LINK_PUBLIC ${fetch_smiles_libs})


set(fetch_smiles_quick_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/fetch_smiles_quick.cc
)
set(fetch_smiles_quick_libs
  accumulator
  cmdline
  iwstring_data_source
  iwstring
  iwre2
  ${IW_EXTRA_LIBS}
)
add_executable(fetch_smiles_quick ${fetch_smiles_quick_srcs})
target_link_libraries(fetch_smiles_quick LINK_PUBLIC ${fetch_smiles_quick_libs})


set(iwsplit_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/iwsplit.cc
)
set(iwsplit_libs
  accumulator
  cmdline
  iwstring_data_source
  iwstring
  iwre2
  ${IW_EXTRA_LIBS}
)
add_executable(iwsplit ${iwsplit_srcs})
target_link_libraries(iwsplit LINK_PUBLIC ${iwsplit_libs})


set(iwcut_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/iwcut.cc
)
set(iwcut_libs
  accumulator
  cmdline
  iwmisc
  iwstring_data_source
  iwstring
  iwre2
  ${IW_EXTRA_LIBS}
)
add_executable(iwcut ${iwcut_srcs})
target_link_libraries(iwcut LINK_PUBLIC ${iwcut_libs})


set(normalise_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/normalise.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwtokeniser.cc
)
set(normalise_libs
  accumulator
  cmdline
  iwstring_data_source
  iwstring
  iwre2
  ${IW_EXTRA_LIBS}
)
add_executable(normalise ${normalise_srcs})
target_link_libraries(normalise LINK_PUBLIC ${normalise_libs})


set(random_records_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/random_records.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwtokeniser.cc
)
set(random_records_libs
  accumulator
  cmdline
  iwstring_data_source
  iwstring
  iwre2
  ${IW_EXTRA_LIBS}
)
add_executable(random_records ${random_records_srcs})
target_link_libraries(random_records LINK_PUBLIC ${random_records_libs})


set(unique_rows_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/unique_rows.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwtokeniser.cc
)
set(unique_rows_libs
  accumulator
  cmdline
  iwstring_data_source
  iwstring
  iwre2
  ${IW_EXTRA_LIBS}
)
add_executable(unique_rows ${unique_rows_srcs})
target_link_libraries(unique_rows LINK_PUBLIC ${unique_rows_libs})

set(proto_src
  ${CMAKE_CURRENT_SOURCE_DIR}/class_label_translation.proto
)
protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS ${proto_src})


set(class_label_translation_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/class_label_translation.cc
  ${PROTO_SRCS}
  ${PROTO_HDRS}
)
set(class_label_translation_libs
  accumulator
  cmdline
  iwstring_data_source
  iwstring
  iwre2
  ${IW_EXTRA_LIBS}
)
add_executable(class_label_translation ${class_label_translation_srcs})
target_link_libraries(class_label_translation LINK_PUBLIC ${class_label_translation_libs} ${Protobuf_LIBRARIES})
