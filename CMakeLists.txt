cmake_minimum_required(VERSION 3.22)

project(SoSim)
include(Simulator/Interface/Private/BuildSystem/sosim_install.cmake)

enable_language(CUDA)
find_package(CUDA REQUIRED)
find_package(Git QUIET)

set(CMAKE_CXX_STANDARD 17)

set(SOSIM_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
set(SOSIM_THIRDPARTY_ROOT ${SOSIM_ROOT}/Simulator/ThirdParty)
set(SOSIM_SRC_ROOT ${SOSIM_ROOT}/Simulator/Source)
set(SOSIM_INCLUDE_DIR
        ${SOSIM_ROOT}/Simulator/Interface
        ${SOSIM_THIRDPARTY_ROOT}/Public
        ${SOSIM_THIRDPARTY_ROOT}/Public/imgui
        ${SOSIM_THIRDPARTY_ROOT}/Public/glad)
set(SOSIM_INSTALL_DIR ${CMAKE_BINARY_DIR} CACHE STRING "SoSim install dir.")

set(SOSIM_FRAMEWORK_LIB "" CACHE STRING "SoSim framework lib(static)." FORCE)
set(SOSIM_SHARED_LIB "" CACHE STRING "SoSim shared libs." FORCE)
set(SOSIM_PHYSICAL_SOLVERS "" CACHE STRING "SoSim solver lib(static)." FORCE)
set(SOSIM_GUI_LIB "" CACHE STRING "SoSim gui lib(static)." FORCE)
set(SOSIM_VERSION_MAJOR 1 CACHE STRING "SoSim version major.")
set(SOSIM_VERSION_MINOR 1 CACHE STRING "SoSim version minor.")
set(SOSIM_VERSION_PATCH 1 CACHE STRING "SoSim version patch.")

execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        RESULT_VARIABLE GIT_SUBMOD_RESULT)
if (NOT GIT_SUBMOD_RESULT EQUAL "0")
    message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
endif ()

add_subdirectory(Simulator/ThirdParty)

add_subdirectory(Simulator/Source)

add_subdirectory(Simulator/Test)

add_executable(SoSim main.cpp)
include_directories(${SOSIM_INCLUDE_DIR} ${CUDA_INCLUDE_DIRS})
target_link_libraries(SoSim PRIVATE
        ${CUDA_LIBRARIES}
        ${SOSIM_FRAMEWORK_LIB}
        ${SOSIM_SHARED_LIB}
        ${SOSIM_PHYSICAL_SOLVERS}
        ${SOSIM_GUI_LIB})

if (MSVC)
    add_custom_command(TARGET SoSim
            POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:${SOSIM_SHARED_LIB}>
            $<TARGET_FILE_DIR:SoSim>)
endif ()


