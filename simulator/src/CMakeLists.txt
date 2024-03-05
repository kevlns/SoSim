project(sosim-libs)

find_package(CUDA REQUIRED)
find_package(assimp CONFIG REQUIRED)

include_directories(${SOSIM_SIMULATOR_ROOT}/include ${CUDA_INCLUDE_DIRS})

# build sosim-libs[static]
file(GLOB_RECURSE LIB_HEADERS ${SOSIM_SIMULATOR_ROOT}/include/libs "*.hpp" "*.cuh")
file(GLOB_RECURSE LIB_SRC libs "*.cpp" "*.cu")
add_library(sosim-libs ${LIB_HEADERS} ${LIB_SRC})
target_link_libraries(sosim-libs PUBLIC ${CUDA_LIBRARIES} assimp::assimp)

# build framework
file(GLOB FRAMEWORK_HEADERS ${SOSIM_SIMULATOR_ROOT}/include/framework/*.hpp)
file(GLOB FRAMEWORK_SRC framework/*.cpp)
add_library(sosim-framework ${FRAMEWORK_HEADERS} ${FRAMEWORK_SRC})
target_link_libraries(sosim-framework PUBLIC sosim-libs)