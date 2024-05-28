project(SoSim-BIBM2024Scene-viewer)

include_directories(${SOSIM_SIMULATOR_ROOT}/include)


add_executable(SoSim-SceneBIBM-MixIMSCT-Viewer mix_imsct.cpp)
target_link_libraries(SoSim-SceneBIBM-MixIMSCT-Viewer PUBLIC
        sosim-rt)

add_executable(SoSim-SceneBIBM-MixIMS-Viewer mix_ims.cpp)
target_link_libraries(SoSim-SceneBIBM-MixIMS-Viewer PUBLIC
        sosim-rt)

add_executable(SoSim-SceneBIBM-MixDFSPH-Viewer mix_dfsph.cpp)
target_link_libraries(SoSim-SceneBIBM-MixDFSPH-Viewer PUBLIC
        sosim-rt)