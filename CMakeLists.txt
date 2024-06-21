cmake_minimum_required(VERSION 3.17)

project(SoSim)

enable_language(CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_ARCHITECTURES 75)
set(SOSIM_SIMULATOR_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/simulator)
set(SORENDER_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/renderer)

set(CMAKE_TOOLCHAIN_FILE ${SOSIM_SIMULATOR_ROOT}/thirdparty/vcpkg/scripts/buildsystems/vcpkg.cmake CACHE FILEPATH "" FORCE)

if (WIN32)
    set(CMAKE_PREFIX_PATH ${SOSIM_SIMULATOR_ROOT}/thirdparty/vcpkg/installed/x64-windows/share CACHE PATH "" FORCE)
endif ()

if (WIN32)
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        list(APPEND CMAKE_PREFIX_PATH "${SOSIM_SIMULATOR_ROOT}/thirdparty/o3d/open3d-lib-win-amd64-0.18.0-debug")
        set(OPEN3D_DLL_PATH ${SOSIM_SIMULATOR_ROOT}/thirdparty/o3d/open3d-lib-win-amd64-0.18.0-debug/bin/Open3D.dll)
    else()
        list(APPEND CMAKE_PREFIX_PATH "${SOSIM_SIMULATOR_ROOT}/thirdparty/o3d/open3d-lib-win-amd64-0.18.0-release")
        set(OPEN3D_DLL_PATH ${SOSIM_SIMULATOR_ROOT}/thirdparty/o3d/open3d-lib-win-amd64-0.18.0-release/bin/Open3D.dll)
    endif()
elseif (UNIX)
    list(APPEND CMAKE_PREFIX_PATH "${SOSIM_SIMULATOR_ROOT}/thirdparty/o3d/open3d-lib-linux-0.18.0")
endif()

add_subdirectory(simulator)

#add_subdirectory(renderer)

add_subdirectory(scenes)






