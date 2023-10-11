# @author        : Long Shen
# @date          : 2023/10/11
# @description   : CMakeList.txt for ...
# @version       : 1.0

cmake_minimum_required(VERSION 3.22)

project(SoSim-DFSPH)

set(TARGET_NAME sosim-dfsph-lib)

FILE(GLOB SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/*.cu)

add_library(${TARGET_NAME} ${SRC_FILES})

target_link_libraries(${TARGET_NAME} PRIVATE sosim-framework-lib sosim-shared-lib)

#find_package(CUDA)