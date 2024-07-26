project(SoSim-SCA-viewer)

include_directories(${SOSIM_SIMULATOR_ROOT}/include)

set(Target SoSim-SCA-Diamonds-Viewer)
add_executable(${Target} diamonds.cpp)
target_link_libraries(${Target} PUBLIC
        sosim-rt)

set(Target SoSim-SCA-Mix-Viewer)
add_executable(${Target} mix.cpp)
target_link_libraries(${Target} PUBLIC
        sosim-rt)

set(Target SoSim-SCA-BucklingParticle-Viewer)
add_executable(${Target} buckling_particle.cpp)
target_link_libraries(${Target} PUBLIC
        sosim-rt)

set(Target SoSim-SCA-BucklingDonut-Viewer)
add_executable(${Target} buckling_donut.cpp)
target_link_libraries(${Target} PUBLIC
        sosim-rt)

set(Target SoSim-SCA-Stir-Viewer)
add_executable(${Target} stir.cpp)
target_link_libraries(${Target} PUBLIC
        sosim-rt)

set(Target SoSim-SCA-ShearBlocks-Viewer)
add_executable(${Target} shear_blocks.cpp)
target_link_libraries(${Target} PUBLIC
        sosim-rt)

set(Target SoSim-SCA-StirCompare-Viewer)
add_executable(${Target} stir_compare.cpp)
target_link_libraries(${Target} PUBLIC
        sosim-rt)

INSTALL(TARGETS SoSim-SCA-Diamonds-Viewer DESTINATION SCA)
INSTALL(TARGETS SoSim-SCA-Mix-Viewer DESTINATION SCA)
INSTALL(TARGETS SoSim-SCA-BucklingParticle-Viewer DESTINATION SCA)
INSTALL(TARGETS SoSim-SCA-BucklingDonut-Viewer DESTINATION SCA)
INSTALL(TARGETS SoSim-SCA-Stir-Viewer DESTINATION SCA)
INSTALL(TARGETS SoSim-SCA-ShearBlocks-Viewer DESTINATION SCA)
INSTALL(TARGETS SoSim-SCA-StirCompare-Viewer DESTINATION SCA)