project(iwstring_data_source)

set(srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/iwstring_string_data_source.cc
)

add_library(iwstring_string_data_source ${srcs})
