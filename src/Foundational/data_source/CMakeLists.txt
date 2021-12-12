project(data_source)

set(iwstring_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/iwstring_data_source.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwmmap.cc
)
add_library(iwstring_data_source ${iwstring_srcs})

set(binary_file_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/binary_data_file.cc
)
add_library(binary_data_file ${binary_file_srcs})
target_include_directories(iwstring_data_source PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

IF(IWTEST)
include(GoogleTest)
enable_testing()
find_package(Protobuf REQUIRED)

protobuf_generate_cpp(for_testing_proto_srcs for_testing_proto_hdrs proto_for_testing.proto)

set(binary_data_file_test_srcs
  binary_data_file_test.cc
  ${for_testing_proto_hdrs} 
  ${for_testing_proto_srcs}
)
set(binary_data_file_test_libs
  binary_data_file
  iwstring
  iwstring_data_source
  iwbits
  gtest
  gmock
  gtest_main
  -lprotobuf -lsnappy -lz
)
add_executable(binary_data_file_test ${binary_data_file_test_srcs})
target_link_libraries(binary_data_file_test LINK_PUBLIC ${binary_data_file_test_libs})

add_test(
  NAME binary_data_file_test
  COMMAND binary_data_file_test
)
gtest_discover_tests(binary_data_file_test)

set(iwstring_data_source_test_srcs
  iwstring_data_source_test.cc
)
set(iwstring_data_source_test_libs
  binary_data_file
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
