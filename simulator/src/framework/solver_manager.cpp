//
// Created by ADMIN on 2024/3/9.
//

#include <memory>

#include "framework/solver_manager.hpp"

namespace SoSim {

    SolverManager::~SolverManager() {
        for (auto solver: m_solvers) {
            solver.reset();
        }
        m_solvers.clear();
    }

    template<class SolverType>
    std::shared_ptr<Solver> SolverManager::createSolver() {
        auto solver = std::make_shared<SolverType>();
        m_solvers.insert(solver);
        return solver;
    }

    void SolverManager::removeSolver(std::shared_ptr<Solver> solver) {
        if (m_solvers.count(solver) > 0) {
            m_solvers.erase(solver);
            solver.reset();
        }
    }

    // explicit template instance
    template std::shared_ptr<Solver> SolverManager::createSolver<DFSPHSolver>();

    template std::shared_ptr<Solver> SolverManager::createSolver<WCSPHSolver>();

    template std::shared_ptr<Solver> SolverManager::createSolver<JL21CTSolver>();

    template std::shared_ptr<Solver> SolverManager::createSolver<IMSCTSolver>();
}
