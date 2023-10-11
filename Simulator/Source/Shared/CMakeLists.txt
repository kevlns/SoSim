# @author        : Long Shen
# @date          : 2023/10/11
# @description   : CMakeList.txt for ...
# @version       : 1.0

cmake_minimum_required(VERSION 3.22)

project(SoSim-Shared)

set(TARGET_NAME sosim-shared-lib)

FILE(GLOB_RECURSE SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR} "*.cpp" "*.cu")

add_library(${TARGET_NAME} ${SRC_FILES})