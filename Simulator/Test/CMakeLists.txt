# @author        : Long Shen
# @date          : 2023/11/21
# @description   : CMakeList.txt for 
# @version       : 1.0

cmake_minimum_required(VERSION 3.11.4)

project(SoSim-Test)

add_subdirectory(UnitTest/ImGUI-Test)

#find_package(CUDA REQUIRED)
#
#file(GLOB SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
#
#set(TARGET_NAME SoSim.test.exe)
#add_executable(${TARGET_NAME} ${SRC_FILES})
#
#include_directories(${SOSIM_INCLUDE_DIR} ${CUDA_INCLUDE_DIRS})
#target_link_libraries(${TARGET_NAME} PRIVATE
#        ${SOSIM_FRAMEWORK_LIBS}
#        ${SOSIM_SHARED_LIBS}
#        ${SOSIM_PHYSICAL_SOLVERS})
#
#target_link_libraries(${TARGET_NAME} PRIVATE
#        glfw)
