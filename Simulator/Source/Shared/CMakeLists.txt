# @author        : Long Shen
# @date          : 2023/10/11
# @description   : CMakeList.txt for ...
# @version       : 1.0

cmake_minimum_required(VERSION 3.22)

project(SoSim-Shared)

enable_language(CUDA)
find_package(CUDA REQUIRED)

# model-tool-lib
set(THIS_TARGET sosim-shared-modeltool-lib)
FILE(GLOB ModelUtils_SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/ModelUtils/*.cpp)
add_library(${THIS_TARGET} ${ModelUtils_SRC_FILES})
target_include_directories(${THIS_TARGET} PRIVATE ${CUDA_INCLUDE_DIRS})
target_link_libraries(${THIS_TARGET} PRIVATE ${CUDA_LIBRARIES})
LIST(APPEND shared-libs ${THIS_TARGET})

# cuda-tool-lib
set(THIS_TARGET sosim-shared-cutool-lib)
FILE(GLOB CudaUtils_SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/CudaUtils/*.cu)
add_library(${THIS_TARGET} ${CudaUtils_SRC_FILES})
target_include_directories(${THIS_TARGET} PRIVATE ${CUDA_INCLUDE_DIRS})
target_link_libraries(${THIS_TARGET} PRIVATE ${CUDA_LIBRARIES})
LIST(APPEND shared-libs ${THIS_TARGET})

# ns-lib
set(THIS_TARGET sosim-shared-ns-lib)
FILE(GLOB NeighborSearchUGB_SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/NeighborSearchUGB/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/NeighborSearchUGB/*.cu)
add_library(${THIS_TARGET} ${NeighborSearchUGB_SRC_FILES})
target_include_directories(${THIS_TARGET} PRIVATE ${CUDA_INCLUDE_DIRS})
target_link_libraries(${THIS_TARGET} PRIVATE ${CUDA_LIBRARIES} sosim-shared-cutool-lib)
LIST(APPEND shared-libs ${THIS_TARGET})

# sphKernel-lib
set(THIS_TARGET sosim-shared-sphKernel-lib)
FILE(GLOB_RECURSE SPHKernel_SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/SPHKernel/*.cu)
add_library(${THIS_TARGET} ${SPHKernel_SRC_FILES})
target_include_directories(${THIS_TARGET} PRIVATE ${CUDA_INCLUDE_DIRS})
target_link_libraries(${THIS_TARGET} PRIVATE ${CUDA_LIBRARIES})
LIST(APPEND shared-libs ${THIS_TARGET})


set(sosim-shared-libs ${shared-libs} PARENT_SCOPE)


