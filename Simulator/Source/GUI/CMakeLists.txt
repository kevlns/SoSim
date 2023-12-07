# @author        : Long Shen
# @date          : 2023/11/21
# @description   : CMakeList.txt for 
# @version       : 1.0

cmake_minimum_required(VERSION 3.11.4)

project(SoSim-GUI)

set(TARGET_NAME sosim-gui)

FILE(GLOB SRC_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cu)
list(APPEND SRC_FILES ${SOSIM_THIRDPARTY_ROOT}/Public/glad/src/glad.c)

FILE(GLOB imgui_srcs ${SOSIM_THIRDPARTY_ROOT}/Public/imgui/src/*.cpp)
list(APPEND SRC_FILES ${imgui_srcs})

add_library(${TARGET_NAME} ${SRC_FILES})

find_package(CUDA REQUIRED)
find_package(OpenGL REQUIRED)
target_include_directories(${TARGET_NAME} PRIVATE ${CUDA_INCLUDE_DIRS})
include_directories(${SOSIM_INCLUDE_DIR} ${CUDA_INCLUDE_DIRS})
target_link_libraries(${TARGET_NAME} PRIVATE
        ${CUDA_LIBRARIES}
        ${SOSIM_FRAMEWORK_LIBS}
        ${SOSIM_SHARED_LIBS}
        ${SOSIM_PHYSICAL_SOLVERS}
        OpenGL::GL
        glfw
        glm)

set(sosim-gui ${TARGET_NAME} PARENT_SCOPE)