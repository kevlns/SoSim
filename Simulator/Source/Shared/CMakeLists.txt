# @author        : Long Shen
# @date          : 2023/10/11
# @description   : CMakeList.txt for ...
# @version       : 1.0

cmake_minimum_required(VERSION 3.22)

project(SoSim-Shared)

enable_language(CUDA)
find_package(CUDA REQUIRED)

set(THIS_TARGET sosim_Shared)
FILE(GLOB_RECURSE SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR} "*.c*" "*.h*")
add_library(${THIS_TARGET} SHARED ${SRC_FILES})
target_include_directories(${THIS_TARGET} PRIVATE ${CUDA_INCLUDE_DIRS})
target_link_libraries(${THIS_TARGET} PRIVATE ${CUDA_LIBRARIES})

target_compile_definitions(${THIS_TARGET} PRIVATE SOSIM_EXPORTS)

set(SOSIM_SHARED_LIB ${THIS_TARGET} CACHE STRING "SoSim add shared lib." FORCE)