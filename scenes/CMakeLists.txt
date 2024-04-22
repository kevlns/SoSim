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
add_executable(SoSim-SceneJL21CT-Viewer scene_jl21ct.cpp)
target_link_libraries(SoSim-SceneJL21CT-Viewer PUBLIC
        sosim-rt)

# ismct scene
add_executable(SoSim-SceneISMCT-Viewer scene_ismct.cpp)
target_link_libraries(SoSim-SceneISMCT-Viewer PUBLIC
        sosim-rt)

add_subdirectory(SCA2024)