# @author        : Long Shen
# @date          : 2023/10/11
# @description   : CMakeList.txt for ...
# @version       : 1.0

cmake_minimum_required(VERSION 3.11.4)

project(SoSim-Framework)

find_package(CUDA REQUIRED)

set(THIS_TARGET sosim_Framework)
FILE(GLOB SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.c*)
add_library(${THIS_TARGET} ${SRC_FILES})
target_include_directories(${THIS_TARGET} PRIVATE ${CUDA_INCLUDE_DIRS})
target_link_libraries(${THIS_TARGET} PRIVATE ${CUDA_LIBRARIES} ${SOSIM_SHARED_LIB})

set(SOSIM_FRAMEWORK_LIB ${THIS_TARGET} CACHE STRING "SoSim add framework lib." FORCE)