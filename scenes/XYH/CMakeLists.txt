project(SoSim-BIBM2024Scene-viewer)

include_directories(${SOSIM_SIMULATOR_ROOT}/include)


add_executable(SoSim-SceneXYH-DFSPH-Viewer dfsph.cpp)
target_link_libraries(SoSim-SceneXYH-DFSPH-Viewer PUBLIC
        sosim-rt)