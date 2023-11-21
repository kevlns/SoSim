# @author        : Long Shen
# @date          : 2023/10/11
# @description   : CMakeList.txt for ...
# @version       : 1.0

cmake_minimum_required(VERSION 3.22)

project(SoSim-Framework)

enable_language(CUDA)
find_package(CUDA REQUIRED)

FILE(GLOB SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# object-lib
set(THIS_TARGET sosim-framework-obj)
add_library(${THIS_TARGET} ${CMAKE_CURRENT_SOURCE_DIR}/object.cpp)
target_include_directories(${THIS_TARGET} PRIVATE ${CUDA_INCLUDE_DIRS})
target_link_libraries(${THIS_TARGET} PRIVATE ${CUDA_LIBRARIES})
LIST(APPEND framework ${THIS_TARGET})

# object-lib
set(THIS_TARGET sosim-framework-solver)
add_library(${THIS_TARGET} ${CMAKE_CURRENT_SOURCE_DIR}/solver.cpp)
target_include_directories(${THIS_TARGET} PRIVATE ${CUDA_INCLUDE_DIRS})
target_link_libraries(${THIS_TARGET} PRIVATE ${CUDA_LIBRARIES} sosim-framework-obj)
LIST(APPEND framework ${THIS_TARGET})

#scene-lib
set(THIS_TARGET sosim-framework-scene)
add_library(${THIS_TARGET} ${CMAKE_CURRENT_SOURCE_DIR}/scene.cpp)
target_include_directories(${THIS_TARGET} PRIVATE ${CUDA_INCLUDE_DIRS})
target_link_libraries(${THIS_TARGET} PRIVATE ${CUDA_LIBRARIES} sosim-framework-solver)
LIST(APPEND framework ${THIS_TARGET})

#simulator-lib
set(THIS_TARGET sosim-framework-simulator)
add_library(${THIS_TARGET} ${CMAKE_CURRENT_SOURCE_DIR}/simulator.cpp)
target_include_directories(${THIS_TARGET} PRIVATE ${CUDA_INCLUDE_DIRS})
target_link_libraries(${THIS_TARGET} PRIVATE ${CUDA_LIBRARIES} sosim-framework-scene)
LIST(APPEND framework ${THIS_TARGET})

set(sosim-framework ${framework} PARENT_SCOPE)
