# @author        : Long Shen
# @date          : 2023/10/11
# @description   : CMakeList.txt for ...
# @version       : 1.0

cmake_minimum_required(VERSION 3.22)

project(SoSim-Framework)

set(TARGET_NAME sosim-framework-lib)

find_package(CUDA REQUIRED)
FILE(GLOB SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_library(${TARGET_NAME} ${SRC_FILES})
target_include_directories(${TARGET_NAME} PRIVATE ${CUDA_INCLUDE_DIRS})
target_link_libraries(${TARGET_NAME} PRIVATE ${CUDA_LIBRARIES})

#find_package(CUDA)