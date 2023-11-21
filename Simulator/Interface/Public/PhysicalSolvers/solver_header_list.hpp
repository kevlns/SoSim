//@author        : Long Shen
//@date          : 2023/11/12
//@description   : include all solver headers
//@version       : 1.0

#ifndef SOSIM_SOLVER_HEADER_LIST_HPP
#define SOSIM_SOLVER_HEADER_LIST_HPP

#include "Public/PhysicalSolvers/MixtureModelForNonNewtonFlow/mmsph_solver.hpp"

namespace SoSim {

    enum SolverInstanceType : uint8_t {
        MMSPH_MMSPHSolver
    };

}

#endif //SOSIM_SOLVER_HEADER_LIST_HPP
