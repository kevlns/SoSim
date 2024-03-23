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

# JL21-CT
file(GLOB JL21CT_HEADERS
        ${SOSIM_SIMULATOR_ROOT}/include/solvers/JL21-ConformationTensor/*.hpp)
file(GLOB JL21CT_SRC
        JL21-ConformationTensor/*.cpp
        JL21-ConformationTensor/*.cuh
        JL21-ConformationTensor/*.cu)
add_library(sosim-jl21ct ${JL21CT_HEADERS} ${JL21CT_SRC})
target_link_libraries(sosim-jl21ct PUBLIC
        sosim-framework)
list(APPEND solver-libs sosim-jl21ct)

set(SOSIM_SOLVER_LIBS ${solver-libs} PARENT_SCOPE)