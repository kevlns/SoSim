//@author        : Long Shen
//@date          : 2023/11/9
//@description   :
//@version       : 1.0

#ifndef SOSIM_MMSPH_SOLVER_TEST_HPP
#define SOSIM_MMSPH_SOLVER_TEST_HPP

#include "Public/Framework/solver.hpp"
#include "Public/PhysicalSolvers/MixtureModelForNonNewtonFlow/mmsph_solver.hpp"

using namespace SoSim;

int mmsph_test() {
    MMSPH::MMSPHSolver solver;
    SolverConfig config{};
    config.block_num = 1;
    config.thread_num = 1;
    config.scene_lb = {-20, -20, -20};
    config.scene_size = {40, 40, 40};
    config.dt = 0.001;

    solver.setSolverConfig(&config);

    solver.addParticles(
            "E:\\Projects\\SoSim\\Simulator\\Source\\PhysicalSolvers\\MixtureModelForNonNewtonFlow\\overall_objs.json");

    return 0;
}

#endif //SOSIM_MMSPH_SOLVER_TEST_HPP
