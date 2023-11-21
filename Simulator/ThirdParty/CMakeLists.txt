# @author        : Long Shen
# @date          : 2023/11/21
# @description   : CMakeList.txt for 
# @version       : 1.0

cmake_minimum_required(VERSION 3.11.4)

include_directories(${SOSIM_INCLUDE_DIR})

project(SoSim-ThirdParty)

# build glfw-3.3.8 and install the includes
add_subdirectory(Private/glfw)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Private/glfw/include/GLFW/ DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/Public/GLFW)

# install imgui includes and src
file(GLOB imgui_headers ${CMAKE_CURRENT_SOURCE_DIR}/Private/imgui/*.h)
list(APPEND imgui_headers
        ${CMAKE_CURRENT_SOURCE_DIR}/Private/imgui/backends/imgui_impl_glfw.h
        ${CMAKE_CURRENT_SOURCE_DIR}/Private/imgui/backends/imgui_impl_opengl3.h
        ${CMAKE_CURRENT_SOURCE_DIR}/Private/imgui/backends/imgui_impl_opengl3_loader.h)
file(COPY ${imgui_headers} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/Public/imgui)

file(GLOB imgui_srcs ${CMAKE_CURRENT_SOURCE_DIR}/Private/imgui/*.cpp)
list(APPEND imgui_srcs
        ${CMAKE_CURRENT_SOURCE_DIR}/Private/imgui/backends/imgui_impl_glfw.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Private/imgui/backends/imgui_impl_opengl3.cpp)
file(COPY ${imgui_srcs} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/Public/imgui/src)

# install json includes
file(GLOB json_headers ${CMAKE_CURRENT_SOURCE_DIR}/Private/json/*.hpp)
file(COPY ${json_headers} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/Public/json)