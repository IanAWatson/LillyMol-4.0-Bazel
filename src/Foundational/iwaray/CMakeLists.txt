cmake_minimum_required(VERSION 3.16.3)
project(iwaray LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

set(srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/_extending_resizable_array_float.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_extending_resizable_array_int.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_extending_resizable_array_long.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_extending_resizable_array_long_long.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_extending_resizable_array_uchar.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_extending_resizable_array_ushort.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_iwaray_float.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_iwaray_int.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_resizable_array_char2.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_resizable_array_char.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_resizable_array_double.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_resizable_array_float.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_resizable_array_int.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_resizable_array_int_star.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_resizable_array_long.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_resizable_array_long_long.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_resizable_array_p_int.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_resizable_array_uchar.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_resizable_array_uint.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_resizable_array_unsigned_long.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_resizable_array_ushort.cc
)

set(hdrs
  ${CMAKE_CURRENT_SOURCE_DIR}/iwaray.h
  ${CMAKE_CURRENT_SOURCE_DIR}/iwaray_np_.h
  ${CMAKE_CURRENT_SOURCE_DIR}/iwaray_op.h
)

file(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/*.cc)
#file(GLOB SRC_FILES  "${CMAKE_CURRENT_SOURCE_DIR}/*.cc")
add_library(iwaray ${srcs})
target_link_libraries(iwaray INTERFACE iwmisc)
target_include_directories (iwaray PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
