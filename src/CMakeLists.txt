cmake_minimum_required(VERSION 3.16.1)
project(LillyMol)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/build)
set(IW_EXTRA_LIBS -lz -lre2)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

enable_testing()

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
#find_package(GTest REQUIRED CONFIG)


#FetchContent_Declare(
#  clif
#  GIT_REPOSITORY https://github.com/google/clif
#  GIT_TAG 25a438fad1534f70c758c3d3c30ec361aeca63e2
#)
#FetchContent_MakeAvailable(clif)

add_compile_options(-march=native)
add_link_options(-lz -lre2)

add_subdirectory(Foundational/iwaray)
add_subdirectory(Foundational/accumulator)
add_subdirectory(Foundational/iwstring)
add_subdirectory(Foundational/cmdline)
add_subdirectory(Foundational/cmdline_v2)
add_subdirectory(Foundational/data_source)
add_subdirectory(Foundational/histogram)
add_subdirectory(Foundational/iwbits)
add_subdirectory(Foundational/iwmisc)
add_subdirectory(Foundational/iwqsort)
add_subdirectory(Foundational/iw_tdt)
add_subdirectory(Foundational/jama)
add_subdirectory(Foundational/string_data_source)
add_subdirectory(Foundational/tnt)
add_subdirectory(Foundational/xmlParser)
add_subdirectory(Utilities/General)
add_subdirectory(Utilities/GFP_Tools)
add_subdirectory(Molecule_Lib)
add_subdirectory(Molecule_Tools)

#configure_file(Foundational/iwaray/iwaray.h iwaray.h)
