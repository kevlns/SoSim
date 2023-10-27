//@author        : Long Shen
//@date          : 2023/10/11
//@description   : 
//@version       : 1.0

#ifndef SOSIM_DFSPH_SOLVER_HPP
#define SOSIM_DFSPH_SOLVER_HPP

#include <iostream>

#include "Public/Framework/solver.hpp"
#include "Public/Shared/NeighborSearchUGB/neighbor_search_ugb.hpp"
#include "Private/DFSPH/data_pack.hpp"

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

            void setParticleRadius(float particle_radius);

            void setSceneInfo(float3 scene_lb, float3 scene_size);

        protected:
            void step() override;

        private:
            SolverConfig m_solverConfig{};
            bool m_isInit{false};
            NeighborSearcher *m_neighborSearcher{nullptr};
            ConstParams m_host_cp{};
            float m_particleRadius{0};

        private:
            ConstParams *m_device_cp{nullptr};
            DynamicParams *m_device_dp{nullptr};
        };
    }
}

#endif //SOSIM_DFSPH_SOLVER_HPP
