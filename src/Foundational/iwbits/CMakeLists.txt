project(iwbits)

set(srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/bits_in_common.c
  ${CMAKE_CURRENT_SOURCE_DIR}/du_bin2ascii.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/fastprint.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/fixed_size_counted_fingerprint.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/hexbits.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwbits.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/swap_byte_order.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/fixed_bit_vector.cc
)

add_library(iwbits ${srcs})
