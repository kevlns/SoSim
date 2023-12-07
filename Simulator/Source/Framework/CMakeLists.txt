# @author        : Long Shen
# @date          : 2023/10/11
# @description   : CMakeList.txt for ...
# @version       : 1.0

cmake_minimum_required(VERSION 3.11.4)

project(SoSim-Framework)

find_package(CUDA REQUIRED)


# object-lib
set(THIS_TARGET sosim-framework-obj)
FILE(GLOB OBJ_SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/SubObjects/*.cpp)
add_library(${THIS_TARGET} ${OBJ_SRC_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/object.cpp)
target_include_directories(${THIS_TARGET} PRIVATE ${CUDA_INCLUDE_DIRS})
target_link_libraries(${THIS_TARGET} PRIVATE ${CUDA_LIBRARIES} ${SOSIM_SHARED_LIBS})
LIST(APPEND framework ${THIS_TARGET})

# solver-lib
set(THIS_TARGET sosim-framework-solver)
add_library(${THIS_TARGET} ${CMAKE_CURRENT_SOURCE_DIR}/solver.cpp)
target_include_directories(${THIS_TARGET} PRIVATE ${CUDA_INCLUDE_DIRS})
target_link_libraries(${THIS_TARGET} PRIVATE ${CUDA_LIBRARIES} ${SOSIM_SHARED_LIBS} ${framework})
LIST(APPEND framework ${THIS_TARGET})

#scene-lib
set(THIS_TARGET sosim-framework-scene)
add_library(${THIS_TARGET} ${CMAKE_CURRENT_SOURCE_DIR}/scene.cpp)
target_include_directories(${THIS_TARGET} PRIVATE ${CUDA_INCLUDE_DIRS})
target_link_libraries(${THIS_TARGET} PRIVATE ${CUDA_LIBRARIES} ${SOSIM_SHARED_LIBS} ${framework})
LIST(APPEND framework ${THIS_TARGET})

#simulator-lib
set(THIS_TARGET sosim-framework-simulator)
add_library(${THIS_TARGET} ${CMAKE_CURRENT_SOURCE_DIR}/simulator.cpp)
target_include_directories(${THIS_TARGET} PRIVATE ${CUDA_INCLUDE_DIRS})
target_link_libraries(${THIS_TARGET} PRIVATE ${CUDA_LIBRARIES} ${SOSIM_SHARED_LIBS} ${framework})
LIST(APPEND framework ${THIS_TARGET})

set(sosim-framework ${framework} PARENT_SCOPE)
