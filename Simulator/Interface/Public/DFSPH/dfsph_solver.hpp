//@author        : Long Shen
//@date          : 2023/10/11
//@description   : 
//@version       : 1.0

#ifndef SOSIM_DFSPH_SOLVER_HPP
#define SOSIM_DFSPH_SOLVER_HPP

#include <iostream>

#include "Public/Framework/solver.hpp"
#include "Private/DFSPH/data_pack.hpp"
#include "Public/Shared/NeighborSearchUGB/neighbor_search_ugb.hpp"

namespace SoSim {

    using NeighborSearcher = NSUGB::NeighborSearchUGB;

    namespace DFSPH {

        class DFSPHSolver : public Solver {
        public:
            using Solver::SolverConfig;

        public:

            DFSPHSolver();

            ~DFSPHSolver() override;

            void initialize() override;

            bool isInitialized() const override;

            void run() override;

            void destroy() override;

            void setConfig(const SolverConfig &config) override;

        protected:
            void step() override;

        private:
            SolverConfig m_solverConfig{};
            bool m_isInit{false};
            NeighborSearcher *m_neighborSearcher{nullptr};

        private:
            ConstParams *m_device_constParams{nullptr};
            DynamicParams *m_device_simParams{nullptr};
        };
    }
}

#endif //SOSIM_DFSPH_SOLVER_HPP
