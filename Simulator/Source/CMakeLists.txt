# @author        : Long Shen
# @date          : 2023/10/11
# @description   : CMakeList.txt for ...
# @version       : 1.0

include_directories(${SOSIM_INCLUDE_DIR})

# build framework
add_subdirectory(Framework)

# build shared-libs
add_subdirectory(Shared)

# build solvers
add_subdirectory(PhysicalSolvers)

# build gui
add_subdirectory(GUI)

target_link_libraries(${SOSIM_FRAMEWORK_LIB} PRIVATE ${sosim-solvers})
target_link_libraries(${SOSIM_GUI_LIB} PRIVATE ${SOSIM_FRAMEWORK_LIB} ${sosim-solvers})