project(SoSim-TVCG-viewer)

include_directories(${SOSIM_SIMULATOR_ROOT}/include)

set(Target SoSim-TVCG-Test-Viewer)
add_executable(${Target} imm_v2_test.cpp)
target_link_libraries(${Target} PUBLIC
        sosim-rt)

INSTALL(TARGETS SoSim-TVCG-Test-Viewer DESTINATION TVCG)