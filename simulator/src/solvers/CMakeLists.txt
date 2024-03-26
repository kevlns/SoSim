# TODO build solvers dll

include_directories(${SOSIM_SIMULATOR_ROOT}/include ${CUDAToolkit_INCLUDE_DIRS})

# dfsph
set(TARGET sosim-dfsph)
file(GLOB DFSPH_HEADERS
        ${SOSIM_SIMULATOR_ROOT}/include/solvers/DFSPH/*.hpp)
file(GLOB DFSPH_SRC
        DFSPH/*.cpp
        DFSPH/*.cuh
        DFSPH/*.cu)
add_library(${TARGET} ${DFSPH_HEADERS} ${DFSPH_SRC})
target_link_libraries(${TARGET} PUBLIC
        sosim-framework)
list(APPEND solver-libs ${TARGET})

# wcsph
set(TARGET sosim-wcsph)
file(GLOB WCSPH_HEADERS
        ${SOSIM_SIMULATOR_ROOT}/include/solvers/WCSPH/*.hpp)
file(GLOB WCSPH_SRC
        WCSPH/*.cpp
        WCSPH/*.cuh
        WCSPH/*.cu)
add_library(${TARGET} ${WCSPH_HEADERS} ${WCSPH_SRC})
target_link_libraries(${TARGET} PUBLIC
        sosim-framework)
list(APPEND solver-libs ${TARGET})

# JL21-CT
set(TARGET sosim-jl21ct)
file(GLOB JL21CT_HEADERS
        ${SOSIM_SIMULATOR_ROOT}/include/solvers/JL21-ConformationTensor/*.hpp)
file(GLOB JL21CT_SRC
        JL21-ConformationTensor/*.cpp
        JL21-ConformationTensor/*.cuh
        JL21-ConformationTensor/*.cu)
add_library(${TARGET} ${JL21CT_HEADERS} ${JL21CT_SRC})
target_link_libraries(${TARGET} PUBLIC
        sosim-framework)
list(APPEND solver-libs ${TARGET})

# IM-CT
set(TARGET sosim-imct)
file(GLOB IMCT_HEADERS
        ${SOSIM_SIMULATOR_ROOT}/include/solvers/IM-ConformationTensor/*.hpp)
file(GLOB IMCT_SRC
        IMS-ConformationTensor/*.cpp
        IMS-ConformationTensor/*.cuh
        IMS-ConformationTensor/*.cu)
add_library(${TARGET} ${IMCT_HEADERS} ${IMCT_SRC})
target_link_libraries(${TARGET} PUBLIC
        sosim-framework)
list(APPEND solver-libs ${TARGET})

set(SOSIM_SOLVER_LIBS ${solver-libs} PARENT_SCOPE)