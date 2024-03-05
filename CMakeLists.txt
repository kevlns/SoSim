cmake_minimum_required(VERSION 3.22)

project(SoSim)

enable_language(CUDA)

set(CMAKE_CXX_STANDARD 17)
set(SOSIM_SIMULATOR_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/simulator)
set(CMAKE_TOOLCHAIN_FILE ${SOSIM_SIMULATOR_ROOT}/thirdparty/vcpkg/scripts/buildsystems/vcpkg.cmake)

add_subdirectory(simulator/src)





