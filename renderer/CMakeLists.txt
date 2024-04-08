cmake_minimum_required(VERSION 3.17)

project(SoRender)

file(GLOB_RECURSE HEADER_FILES include "*.hpp" "*.cuh")
file(GLOB_RECURSE SRC_FILES src "*.cpp" "*.cu")

include_directories(include ${SOSIM_SIMULATOR_ROOT}/include)

set(TARGET sorender)
add_library(${TARGET} ${HEADER_FILES} ${SRC_FILES})
target_link_libraries(${TARGET} PUBLIC sosim-rt)