project(SoSim-SCA2024Scene-viewer)

include_directories(${SOSIM_SIMULATOR_ROOT}/include)

# ismct scene
add_executable(SoSim-SceneISMCT-Buckling-Viewer scene_buckling.cpp)
target_link_libraries(SoSim-SceneISMCT-Buckling-Viewer PUBLIC
        sosim-rt)

add_executable(SoSim-SceneISMCT-Diamonds-Viewer scene_diamonds.cpp)
target_link_libraries(SoSim-SceneISMCT-Diamonds-Viewer PUBLIC
        sosim-rt)

add_executable(SoSim-SceneISMCT-Stir-Viewer scene_stir.cpp)
target_link_libraries(SoSim-SceneISMCT-Stir-Viewer PUBLIC
        sosim-rt)

add_executable(SoSim-SceneISMCT-Diff-Viewer scene_diff.cpp)
target_link_libraries(SoSim-SceneISMCT-Diff-Viewer PUBLIC
        sosim-rt)

add_executable(SoSim-SceneISMCT-Diff2-Viewer scene_diff_2.cpp)
target_link_libraries(SoSim-SceneISMCT-Diff2-Viewer PUBLIC
        sosim-rt)

add_executable(SoSim-SceneISMCT-Shear-Viewer scene_shear.cpp)
target_link_libraries(SoSim-SceneISMCT-Shear-Viewer PUBLIC
        sosim-rt)