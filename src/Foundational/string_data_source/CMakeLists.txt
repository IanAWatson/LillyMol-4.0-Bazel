project(string_data_source)

set(srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/iwstring_string_data_source.cc
)

add_library(string_data_source ${srcs})
target_include_directories(string_data_source PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
