cmake_minimum_required(VERSION 3.22)

project(SoSim)

enable_language(CUDA)

set(CMAKE_CXX_STANDARD 17)

set(SOSIM_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
set(SOSIM_INCLUDE_DIR ${SOSIM_ROOT}/Simulator/Interface)
set(SOSIM_SRC_ROOT ${SOSIM_ROOT}/Simulator/Source)
set(SOSIM_SHARED_LIBS "")

add_subdirectory(Simulator/Source)

find_package(CUDA REQUIRED)

set(TARGET_NAME SoSim)
add_executable(${TARGET_NAME} main.cpp)
include_directories(${SOSIM_INCLUDE_DIR} ${CUDA_INCLUDE_DIRS})
target_link_libraries(${TARGET_NAME} PRIVATE sosim-framework
        ${SOSIM_SHARED_LIBS}
        sosim-mmsph)

