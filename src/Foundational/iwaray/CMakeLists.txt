cmake_minimum_required(VERSION 3.16.3)
project(iwaray)
set(CMAKE_CXX_STANDARD 17)

file(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/*.cc)
#file(GLOB SRC_FILES "*.cc")
add_library(iwaray ${SRC_FILES})
target_link_libraries(iwaray PUBLIC iwconfig)
