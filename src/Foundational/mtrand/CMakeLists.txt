project(iwrandom)
set(srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/iwrandom.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/mtrand.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/normal_distribution.cc
)

add_library(iwrandom ${srcs})
