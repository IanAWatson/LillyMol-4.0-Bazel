project(data_source)

set(iwstring_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/iwstring_data_source.cc
)
add_library(iwstring_data_source ${iwstring_srcs})

set(binary_file_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/binary_data_file.cc
)
add_library(binary_data_file ${binary_file_srcs})
target_include_directories(iwstring_data_source PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
