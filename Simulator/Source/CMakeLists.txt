# @author        : Long Shen
# @date          : 2023/10/11
# @description   : CMakeList.txt for ...
# @version       : 1.0

include_directories(${SOSIM_INCLUDE_DIR})

set(sosim-framework "")
add_subdirectory(Framework)
set(SOSIM_FRAMEWORK_LIBS ${sosim-framework} PARENT_SCOPE)

set(sosim-shared "")
add_subdirectory(Shared)
set(SOSIM_SHARED_LIBS ${sosim-shared} PARENT_SCOPE)

set(sosim-solvers "")
add_subdirectory(PhysicalSolvers)
set(SOSIM_PHYSICAL_SOLVERS ${sosim-solvers} PARENT_SCOPE)


target_link_libraries(sosim-framework-scene PUBLIC ${sosim-solvers})