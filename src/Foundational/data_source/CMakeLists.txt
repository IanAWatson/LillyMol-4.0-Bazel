project(data_source)

set(iwstring_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/iwstring_data_source.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwmmap.cc
)
add_library(iwstring_data_source ${iwstring_srcs})

target_include_directories(iwstring_data_source PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

IF(IWTEST)
include(GoogleTest)
enable_testing()
find_package(Protobuf REQUIRED)

protobuf_generate_cpp(for_testing_proto_srcs for_testing_proto_hdrs proto_for_testing.proto)

set(iwstring_data_source_test_srcs
  iwstring_data_source_test.cc
)
set(iwstring_data_source_test_libs
  iwstring_data_source
  iwstring
  iwre2
  gtest
  gmock
  gtest_main
  -lre2 -lz
)
add_executable(iwstring_data_source_test ${iwstring_data_source_test_srcs})
target_link_libraries(iwstring_data_source_test LINK_PUBLIC ${iwstring_data_source_test_libs})

add_test(
  NAME iwstring_data_source_test
  COMMAND iwstring_data_source_test
)

gtest_discover_tests(iwstring_data_source_test)

endif(IWTEST)
