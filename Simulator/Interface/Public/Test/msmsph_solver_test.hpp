//@author        : Long Shen
//@date          : 2023/10/26
//@description   :
//@version       : 1.0

#ifndef SOSIM_MSMSPH_SOLVER_TEST_HPP
#define SOSIM_MSMSPH_SOLVER_TEST_HPP

#include <memory>

#include "Public/Framework/solver.hpp"
#include "Public/MSMSPH/msmsph_solver.hpp"

using namespace SoSim;

int MSMSPH_Test() {

    MSMSPH::MSMSPHSolver solver;

    solver.addParts("E:\\Projects\\SoSim\\Simulator\\Source\\MSMSPH\\ObjectConfig\\overall_objs.json");

    MSMSPH::MSMSPHSolver::SolverConfig config;
    config.dt = 0.01;
    config.gravity = {0, -9.8, 0};
    solver.setConfig(config);

    solver.initialize();

    solver.run();

    solver.destroy();

    return 0;
}


#endif //SOSIM_MSMSPH_SOLVER_TEST_HPP
