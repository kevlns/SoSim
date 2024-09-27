project(SoSim-CMPB-viewer)

include_directories(${SOSIM_SIMULATOR_ROOT}/include)

add_executable(SoSim-CMPB-Inject-Viewer emitter_test.cpp)
target_link_libraries(SoSim-CMPB-Inject-Viewer PUBLIC
        sosim-rt)

INSTALL(TARGETS SoSim-CMPB-Inject-Viewer DESTINATION BIBM)