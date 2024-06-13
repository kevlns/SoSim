//
// Created by ADMIN on 2024/3/9.
//

#ifndef SOSIM_SOLVER_MANAGER_HPP
#define SOSIM_SOLVER_MANAGER_HPP

#include <set>

// solvers includes
#include "solvers/DFSPH/dfsph_solver.hpp"
#include "solvers/IMM/imm_solver.hpp"
#include "solvers/IMM-CT/immct_solver.hpp"
#include "solvers/PBF/pbf_solver.hpp"

namespace SoSim {

    class SolverManager {
    public:
        ~SolverManager();

        template<class SolverType>
        std::shared_ptr<Solver> createSolver();

        void removeSolver(std::shared_ptr<Solver> solver);

    private:
        std::set<std::shared_ptr<Solver>> m_solvers;
    };

}

#endif //SOSIM_SOLVER_MANAGER_HPP
