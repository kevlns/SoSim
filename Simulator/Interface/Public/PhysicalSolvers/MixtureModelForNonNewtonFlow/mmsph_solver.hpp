//@author        : Long Shen
//@date          : 2023/11/8
//@description   :
//@version       : 1.0

#ifndef SOSIM_MMSPH_SOLVER_HPP
#define SOSIM_MMSPH_SOLVER_HPP

#include <vector>
#include <iostream>

#include "Public/Framework/solver.hpp"
#include "Public/Shared/NeighborSearchUGB/neighbor_search_ugb.hpp"
#include "Private/PhysicalSolvers/MixtureModelForNonNewtonFlow/data_pack.hpp"
#include "Public/PhysicalSolvers/MixtureModelForNonNewtonFlow/enum_type.hpp"

namespace SoSim {

    using NeighborSearcher = NSUGB::NeighborSearcher;

    namespace MMSPH {

        // if you need more config, create a new config data struct then inherit the SolverConfig
        /*
             struct MMSPHSolverConfig : public SolverConfig {
                ...
             }
         */
        // ...

        class MMSPHSolver : public Solver {

        public:
            MMSPHSolver() = default;

            ~MMSPHSolver() override = default;

            void initialize() override;

            bool isInitialized() const override;

            void run() override;

            void destroy() override;

            void setSolverConfig(const SolverConfig *config) override;

            void attachObject(Object *obj) override;

            void addParticles(const std::string &obj_json_path) override;

        protected:
            void step() override;

        private:
            bool m_isInit{false};
            double m_mem{0};

            SolverConfig *m_config{nullptr};
            NeighborSearcher *m_neighborSearcher{nullptr};
            ConstParams m_host_cp;
            DynamicParams m_host_dp;
            ConstParams *m_device_cp;
            DynamicParams *m_device_dp;

            std::vector<float3> m_host_pos;
            std::vector<float3> m_host_vel;
            std::vector<float2> m_host_phase;
            std::vector<Material> m_host_mat;
            std::vector<bool> m_host_isDynamic;


        };

    }

}

#endif //SOSIM_MMSPH_SOLVER_HPP
