# TODO build solvers dll

include_directories(${SOSIM_SIMULATOR_ROOT}/include ${CUDAToolkit_INCLUDE_DIRS})

# dfsph
file(GLOB DFSPH_HEADERS
        ${SOSIM_SIMULATOR_ROOT}/include/solvers/DFSPH/*.hpp)
file(GLOB DFSPH_SRC
        DFSPH/*.cpp
        DFSPH/*.cuh
        DFSPH/*.cu)
add_library(sosim-dfsph ${DFSPH_HEADERS} ${DFSPH_SRC})
target_link_libraries(sosim-dfsph PUBLIC
        sosim-framework)
list(APPEND solver-libs sosim-dfsph)

# wcsph
file(GLOB WCSPH_HEADERS
        ${SOSIM_SIMULATOR_ROOT}/include/solvers/WCSPH/*.hpp)
file(GLOB WCSPH_SRC
        WCSPH/*.cpp
        WCSPH/*.cuh
        WCSPH/*.cu)
add_library(sosim-wcsph ${WCSPH_HEADERS} ${WCSPH_SRC})
target_link_libraries(sosim-wcsph PUBLIC
        sosim-framework)
list(APPEND solver-libs sosim-wcsph)

set(SOSIM_SOLVER_LIBS ${solver-libs} PARENT_SCOPE)