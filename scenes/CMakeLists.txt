project(SoSim-scene-viewer)

include_directories(${SOSIM_SIMULATOR_ROOT}/include)

# base scene
add_executable(SoSim-SceneBase-Viewer scene_base.cpp)
target_link_libraries(SoSim-SceneBase-Viewer PUBLIC
        sosim-rt)

# wcsph scene
add_executable(SoSim-SceneWCSPH-Viewer scene_wcsph.cpp)
target_link_libraries(SoSim-SceneWCSPH-Viewer PUBLIC
        sosim-rt)

# mixct scene
add_executable(SoSim-SceneMixCT-Viewer scene_mixct.cpp)
target_link_libraries(SoSim-SceneMixCT-Viewer PUBLIC
        sosim-rt)