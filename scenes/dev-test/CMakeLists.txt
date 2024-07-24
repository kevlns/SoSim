project(SoSim-DevTest-viewer)

include_directories(${SOSIM_SIMULATOR_ROOT}/include)

set(Target SoSim-DevTest-Matrix_v2-Viewer)
add_executable(${Target} matrix_v2_test.cpp)
target_link_libraries(${Target} PUBLIC
        sosim-rt)