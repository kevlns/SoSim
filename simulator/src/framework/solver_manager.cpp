//
// Created by ADMIN on 2024/3/9.
//

#include "framework/solver_manager.hpp"

namespace SoSim {

    template<typename SolverType>
    Solver *SolverManager::createSolver() {
        Solver *solver = new SolverType;
        m_solvers.insert(solver);
        return solver;
    }

    void SolverManager::removeSolver(Solver *solver) {
        if (m_solvers.count(solver) > 0) {
            solver->destroy();
            m_solvers.erase(solver);
            delete solver;
        }
    }

    void SolverManager::destroy() {
        for (auto solver: m_solvers) {
            solver->destroy();
            delete solver;
        }
        m_solvers.clear();
    }

    // explicit template instance
    template Solver *SolverManager::createSolver<DFSPHSolver>();
}
