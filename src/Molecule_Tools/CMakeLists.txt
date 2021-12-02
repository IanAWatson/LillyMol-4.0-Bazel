project(molecule_tools)

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
  iwmolecule
)
add_executable(fileconv ${fileconv_srcs})
target_link_libraries(fileconv LINK_PUBLIC ${fileconv_libs})
