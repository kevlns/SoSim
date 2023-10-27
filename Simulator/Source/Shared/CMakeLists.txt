# @author        : Long Shen
# @date          : 2023/10/11
# @description   : CMakeList.txt for ...
# @version       : 1.0

cmake_minimum_required(VERSION 3.22)

project(SoSim-Shared)

enable_language(CUDA)
find_package(CUDA REQUIRED)

# ns-lib
FILE(GLOB NeighborSearchUGB_SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/NeighborSearchUGB/*.cpp)
add_library(sosim-shared-ns-lib ${NeighborSearchUGB_SRC_FILES})
target_include_directories(sosim-shared-ns-lib PRIVATE ${CUDA_INCLUDE_DIRS})
target_link_libraries(sosim-shared-ns-lib PRIVATE ${CUDA_LIBRARIES})
LIST(APPEND shared-libs sosim-shared-ns-lib)

# sphKernel-lib
FILE(GLOB_RECURSE SPHKernel_SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/SPHKernel/*.cu)
add_library(sosim-shared-sphKernel-lib ${SPHKernel_SRC_FILES})
target_include_directories(sosim-shared-sphKernel-lib PRIVATE ${CUDA_INCLUDE_DIRS})
target_link_libraries(sosim-shared-sphKernel-lib PRIVATE ${CUDA_LIBRARIES})
LIST(APPEND shared-libs sosim-shared-sphKernel-lib)


set(sosim-shared-libs ${shared-libs} PARENT_SCOPE)


