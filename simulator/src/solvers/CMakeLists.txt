# TODO build solvers dll

include_directories(${SOSIM_SIMULATOR_ROOT}/include ${CUDAToolkit_INCLUDE_DIRS})

find_package(Bullet CONFIG REQUIRED)

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

# IMM
set(TARGET sosim-imm)
file(GLOB IMM_HEADERS
        ${SOSIM_SIMULATOR_ROOT}/include/solvers/IMM/*.hpp)
file(GLOB IMM_SRC
        IMM/*.cpp
        IMM/*.cuh
        IMM/*.cu)
add_library(${TARGET} ${IMM_HEADERS} ${IMM_SRC})
target_link_libraries(${TARGET} PUBLIC
        sosim-framework
        ${BULLET_LIBRARIES})
list(APPEND solver-libs ${TARGET})

# IMM-CT
set(TARGET sosim-immct)
file(GLOB IMMCT_HEADERS
        ${SOSIM_SIMULATOR_ROOT}/include/solvers/IMM-CT/*.hpp)
file(GLOB IMMCT_SRC
        IMM-CT/*.cpp
        IMM-CT/*.cuh
        IMM-CT/*.cu)
add_library(${TARGET} ${IMMCT_HEADERS} ${IMMCT_SRC})
target_link_libraries(${TARGET} PUBLIC
        sosim-framework
        ${BULLET_LIBRARIES})
list(APPEND solver-libs ${TARGET})

# IMM-v2
set(TARGET sosim-imm-v2)
file(GLOB IMM_HEADERS
        ${SOSIM_SIMULATOR_ROOT}/include/solvers/IMM-v2/*.hpp)
file(GLOB IMM_SRC
        IMM-v2/*.cpp
        IMM-v2/*.cuh
        IMM-v2/*.cu)
add_library(${TARGET} ${IMM_HEADERS} ${IMM_SRC})
target_link_libraries(${TARGET} PUBLIC
        sosim-framework
        ${BULLET_LIBRARIES})
list(APPEND solver-libs ${TARGET})

# PBF
set(TARGET sosim-pbf)
file(GLOB PBF_HEADERS
        ${SOSIM_SIMULATOR_ROOT}/include/solvers/PBF/*.hpp)
file(GLOB PBF_SRC
        PBF/*.cpp
        PBF/*.cuh
        PBF/*.cu)
add_library(${TARGET} ${PBF_HEADERS} ${PBF_SRC})
target_link_libraries(${TARGET} PUBLIC
        sosim-framework
        ${BULLET_LIBRARIES})
list(APPEND solver-libs ${TARGET})

# VisPBF
set(TARGET sosim-vis_pbf)
file(GLOB VisPBF_HEADERS
        ${SOSIM_SIMULATOR_ROOT}/include/solvers/VisPBF/*.hpp)
file(GLOB VisPBF_SRC
        VisPBF/*.cpp
        VisPBF/*.cuh
        VisPBF/*.cu)
add_library(${TARGET} ${VisPBF_HEADERS} ${VisPBF_SRC})
target_link_libraries(${TARGET} PUBLIC
        sosim-framework
        ${BULLET_LIBRARIES})
list(APPEND solver-libs ${TARGET})

## wcsph
#set(TARGET sosim-wcsph)
#file(GLOB WCSPH_HEADERS
#        ${SOSIM_SIMULATOR_ROOT}/include/solvers/WCSPH/*.hpp)
#file(GLOB WCSPH_SRC
#        WCSPH/*.cpp
#        WCSPH/*.cuh
#        WCSPH/*.cu)
#add_library(${TARGET} ${WCSPH_HEADERS} ${WCSPH_SRC})
#target_link_libraries(${TARGET} PUBLIC
#        sosim-framework)
#list(APPEND solver-libs ${TARGET})
#
## JL21-CT
#set(TARGET sosim-jl21ct)
#file(GLOB JL21CT_HEADERS
#        ${SOSIM_SIMULATOR_ROOT}/include/solvers/JL21-ConformationTensor/*.hpp)
#file(GLOB JL21CT_SRC
#        JL21-ConformationTensor/*.cpp
#        JL21-ConformationTensor/*.cuh
#        JL21-ConformationTensor/*.cu)
#add_library(${TARGET} ${JL21CT_HEADERS} ${JL21CT_SRC})
#target_link_libraries(${TARGET} PUBLIC
#        sosim-framework)
#list(APPEND solver-libs ${TARGET})

set(SOSIM_SOLVER_LIBS ${solver-libs} PARENT_SCOPE)