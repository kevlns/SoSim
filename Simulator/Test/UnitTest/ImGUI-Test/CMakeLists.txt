# @author        : Long Shen
# @date          : 2023/11/21
# @description   : CMakeList.txt for 
# @version       : 1.0

cmake_minimum_required(VERSION 3.11.4)

project(SoSim-ImGUI-Test)

find_package(CUDA REQUIRED)
find_package(OpenGL REQUIRED)

file(GLOB imgui_srcs ${SOSIM_THIRDPARTY_ROOT}/Public/imgui/src/*.cpp)

set(TARGET_NAME SoSim-test-imgui.exe)
add_executable(${TARGET_NAME} imgui_test.cpp ${imgui_srcs})
include_directories(${SOSIM_INCLUDE_DIR} ${CUDA_INCLUDE_DIRS})
target_link_libraries(${TARGET_NAME} PRIVATE
        ${CUDA_LIBRARIES}
        ${SOSIM_FRAMEWORK_LIBS}
        ${SOSIM_SHARED_LIBS}
        ${SOSIM_PHYSICAL_SOLVERS}
        glfw
        OpenGL::GL)