# @author        : Long Shen
# @date          : 2023/11/8
# @description   : CMakeList.txt for 
# @version       : 1.0

cmake_minimum_required(VERSION 3.11.4)

set(solvers "")
set(current_solver "")

add_subdirectory(MixtureModelForNonNewtonFlow)
list(APPEND solvers ${current_solver})

add_subdirectory(SampleSolver)
list(APPEND solvers ${current_solver})

#add_subdirectory(NEXT_SOLVER)
#list(APPEND solvers current_solver)


set(sosim-solvers ${solvers} PARENT_SCOPE)