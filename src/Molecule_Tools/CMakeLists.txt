project(molecule_tools)

set(third_party_libs
  -lz -lre2 -lprotobuf -lsnappy
)

set(fileconv_srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/fileconv.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/fix_structures.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/fix_structures.h
  ${CMAKE_CURRENT_SOURCE_DIR}/do_remove_duplicate_fragments.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/do_remove_duplicate_fragments.h
  ${CMAKE_CURRENT_SOURCE_DIR}/known_fragment_data.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/known_fragment_data.h
)

set(fileconv_libs
  cmdline
  accumulator
  moleculeio
  iwmolecule
  iwmisc
  iwbits
  binary_data_file
  iwstring_data_source
  string_data_source
  iwre2
  iwstring
  iwrandom
  xmlParser
  ${third_party_libs}
)
add_executable(fileconv ${fileconv_srcs})
target_link_libraries(fileconv LINK_PUBLIC ${fileconv_libs})
