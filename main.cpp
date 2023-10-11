#include <memory>

#include "Public/DFSPH/dfsph_solver.hpp"

using namespace SoSim;

int main() {
    std::shared_ptr<Solver> solver = std::make_shared<DFSPH::DFSPHSolver>();

    Solver::SolverConfig config{0.001, 0};
    solver->setConfig(config);
}