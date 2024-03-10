project(sosim-libs)

find_package(CUDAToolkit)
find_package(assimp CONFIG REQUIRED)

include_directories(${SOSIM_SIMULATOR_ROOT}/include ${CUDAToolkit_INCLUDE_DIRS})

# build sosim-libs[static]
file(GLOB_RECURSE LIB_HEADERS ${SOSIM_SIMULATOR_ROOT}/include/libs "*.hpp" "*.cuh")
file(GLOB_RECURSE LIB_SRC libs "*.cpp" "*.cu")
add_library(sosim-libs ${LIB_HEADERS} ${LIB_SRC})
target_link_libraries(sosim-libs PUBLIC assimp::assimp CUDA::cudart)

# build base framework
file(GLOB FRAMEWORK_HEADERS
        ${SOSIM_SIMULATOR_ROOT}/include/framework/object.hpp
        ${SOSIM_SIMULATOR_ROOT}/include/framework/solver.hpp)
file(GLOB FRAMEWORK_SRC
        framework/object.cpp
        framework/solver.cpp)
add_library(sosim-framework ${FRAMEWORK_HEADERS} ${FRAMEWORK_SRC})
target_link_libraries(sosim-framework PUBLIC sosim-libs)

# TODO build solvers dll
file(GLOB DFSPH_HEADERS
        ${SOSIM_SIMULATOR_ROOT}/include/solvers/DFSPH/*.hpp)
file(GLOB DFSPH_SRC
        solvers/DFSPH/*.cpp
        solvers/DFSPH/*.cuh
        solvers/DFSPH/*.cu)
add_library(sosim-dfsph ${DFSPH_HEADERS} ${DFSPH_SRC})
target_link_libraries(sosim-dfsph PUBLIC sosim-framework)

# TODO build runtime lib/dll
file(GLOB RT_HEADERS
        ${SOSIM_SIMULATOR_ROOT}/include/framework/object_manager.hpp
        ${SOSIM_SIMULATOR_ROOT}/include/framework/solver_manager.hpp)
file(GLOB RT_SRC
        framework/object_manager.cpp
        framework/solver_manager.cpp)
add_library(sosim-rt ${RT_HEADERS} ${RT_SRC})
target_link_libraries(sosim-rt PRIVATE sosim-dfsph)