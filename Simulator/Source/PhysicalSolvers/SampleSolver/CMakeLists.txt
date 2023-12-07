# @author        : Long Shen
# @date          : 2023/11/26
# @description   : CMakeList.txt for 
# @version       : 1.0

cmake_minimum_required(VERSION 3.11.4)

project(sosim-sampleSolver)

find_package(CUDA REQUIRED)

FILE(GLOB SRC_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cu)

set(TARGET_NAME sosim-solver-sample)
add_library(${TARGET_NAME} ${SRC_FILES})
include_directories(${SOSIM_INCLUDE_DIR} ${CUDA_INCLUDE_DIRS})
target_link_libraries(${TARGET_NAME} PRIVATE ${CUDA_LIBRARIES} ${SOSIM_FRAMEWORK_LIBS} ${SOSIM_SHARED_LIBS})

set(current_solver ${TARGET_NAME} PARENT_SCOPE)
