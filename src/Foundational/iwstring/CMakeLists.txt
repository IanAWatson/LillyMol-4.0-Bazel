project(iwstring)

set(srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/append_number.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/hexprint.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwgetline.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwstring.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/IWString_class.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwstring_and_file_descriptor.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwwrite.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwwrite_block.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwstrncasecmp.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwstring_strstream.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwzlib.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iw_stl_hash_map.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/remove_suffix.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/string_change.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/string_relationals.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/unhtml.cc
)

add_library(iwstring ${srcs})
target_include_directories (iwstring PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
