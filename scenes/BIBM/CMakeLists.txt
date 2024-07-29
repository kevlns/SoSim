project(SoSim-BIBM-viewer)

include_directories(${SOSIM_SIMULATOR_ROOT}/include)

add_executable(SoSim-BIBM-Inject-Viewer inject.cpp)
target_link_libraries(SoSim-BIBM-Inject-Viewer PUBLIC
        sosim-rt)

add_executable(SoSim-BIBM-EmitterTest-Viewer emitter_test.cpp)
target_link_libraries(SoSim-BIBM-EmitterTest-Viewer PUBLIC
        sosim-rt)

add_executable(SoSim-BIBM-Mix-Viewer mix.cpp)
target_link_libraries(SoSim-BIBM-Mix-Viewer PUBLIC
        sosim-rt)

INSTALL(TARGETS SoSim-BIBM-Inject-Viewer DESTINATION BIBM)
INSTALL(TARGETS SoSim-BIBM-EmitterTest-Viewer DESTINATION BIBM)
INSTALL(TARGETS SoSim-BIBM-Mix-Viewer DESTINATION BIBM)