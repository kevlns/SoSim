//
// Created by ADMIN on 2024/3/9.
//

#ifndef SOSIM_SOLVER_MANAGER_HPP
#define SOSIM_SOLVER_MANAGER_HPP

#include <set>

// solvers includes
#include "solvers/DFSPH/dfsph_solver.hpp"

namespace SoSim {

    class SolverManager {
    public:
        template<class SolverType>
        Solver *createSolver();

        void removeSolver(Solver *solver);

        void destroy();

    private:
        std::set<Solver *> m_solvers;
    };

}

#endif //SOSIM_SOLVER_MANAGER_HPP
