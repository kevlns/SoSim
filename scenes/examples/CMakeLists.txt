project(SoSim-Examples-viewer)

include_directories(${SOSIM_SIMULATOR_ROOT}/include)

add_executable(SoSim-Examples-GenObj-Viewer example_create_object.cpp)
target_link_libraries(SoSim-Examples-GenObj-Viewer PUBLIC
        sosim-rt)

add_executable(SoSim-Examples-DFSPH-Viewer example_dfsph.cpp)
target_link_libraries(SoSim-Examples-DFSPH-Viewer PUBLIC
        sosim-rt)

add_executable(SoSim-Examples-IMM-Viewer example_imm.cpp)
target_link_libraries(SoSim-Examples-IMM-Viewer PUBLIC
        sosim-rt)

add_executable(SoSim-Examples-IMMCT-Viewer example_immct.cpp)
target_link_libraries(SoSim-Examples-IMMCT-Viewer PUBLIC
        sosim-rt)

if(WIN32)
    message(STATUS ${OPEN3D_DLL_PATH})
    add_custom_command(TARGET SoSim-Examples-GenObj-Viewer POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${OPEN3D_DLL_PATH}
            $<TARGET_FILE_DIR:SoSim-Examples-GenObj-Viewer>
    )
endif ()
