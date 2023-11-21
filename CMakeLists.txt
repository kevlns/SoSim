cmake_minimum_required(VERSION 3.22)

project(SoSim)

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
        ${SOSIM_THIRDPARTY_ROOT}/Public/imgui)
set(SOSIM_SHARED_LIBS "")
set(SOSIM_FRAMEWORK_LIBS "")
set(SOSIM_PHYSICAL_SOLVERS "")

execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        RESULT_VARIABLE GIT_SUBMOD_RESULT)
if(NOT GIT_SUBMOD_RESULT EQUAL "0")
    message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
endif()

add_subdirectory(Simulator/ThirdParty)

add_subdirectory(Simulator/Source)

add_subdirectory(Simulator/Test)


