# @author        : Long Shen
# @date          : 2023/10/24
# @description   : CMakeList.txt for Mixture-SuspensionLiquid model, so called MSMSPH
# @version       : 1.0

cmake_minimum_required(VERSION 3.11.4)

project(SoSim-MSMSPH)

set(TARGET_NAME sosim-msmsph-lib)

FILE(GLOB SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/*.cu)

add_library(${TARGET_NAME} ${SRC_FILES})

find_package(CUDA REQUIRED)
target_include_directories(${TARGET_NAME} PRIVATE ${CUDA_INCLUDE_DIRS})
target_link_libraries(${TARGET_NAME} PRIVATE ${CUDA_LIBRARIES} sosim-framework-lib ${sosim-shared-libs})