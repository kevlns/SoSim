cmake_minimum_required(VERSION 3.22)

project(SoSim)

#enable_language(CUDA)

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_OSX_DEPLOYMENT_TARGET 13.0)

set(SOSIM_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
set(SOSIM_INCLUDE_DIR ${SOSIM_ROOT}/Simulator/Interface)
set(SOSIM_SRC_ROOT ${SOSIM_ROOT}/Simulator/Source)

add_subdirectory(Simulator/Source)

set(TARGET_NAME SoSim)
add_executable(${TARGET_NAME} main.cpp)
include_directories(${SOSIM_INCLUDE_DIR})
target_link_libraries(${TARGET_NAME} PRIVATE sosim-framework-lib
        sosim-dfsph-lib)
