//@author        : Long Shen
//@date          : 2023/10/24
//@description   :
//@version       : 1.0

#ifndef SOSIM_MSMSPH_SOLVER_HPP
#define SOSIM_MSMSPH_SOLVER_HPP

#include <vector>

#include "Public/Framework/solver.hpp"
#include "Public/Shared/NeighborSearchUGB/neighbor_search_ugb.hpp"
#include "Private/MSMSPH/data_pack.hpp"

namespace SoSim {

    using NeighborSearcher = NSUGB::NeighborSearcher;

    namespace MSMSPH {

        class MSMSPHSolver : public Solver {
        public:
            using Solver::SolverConfig;

        public:

            MSMSPHSolver() = default;

            ~MSMSPHSolver() override = default;

            void initialize() override;

            bool isInitialized() const override;

            void run() override;

            void destroy() override;

            void setConfig(const SolverConfig &config) override;

            void attachObject() override;

//            void addParticles(const std::string &obj_json) override;

            void addParts(const std::string &obj_json) override;

            void setPhaseDensity(float den1, float den2);

            void setParticleRadius(float particle_radius);

            void setSceneInfo(float3 scene_lb, float3 scene_size);

            void setBNeighborSearcherPtr(NeighborSearcher *bNS);

            NeighborSearcher *getNeighborSearcherPtr() const;

        protected:
            void step() override;

        private:
            void initCP();

            void initDP();

        private:
            ConstParams m_host_cp{};
            DynamicParams m_host_dp{};
            std::vector<float3> m_host_pos;
            std::vector<float3> m_host_vel;
            std::vector<float> m_host_den;
            std::vector<Material> m_host_mat;
            std::vector<Phase> m_host_oPhase;

        private:
            SolverConfig m_solverConfig{};
            bool m_isInit{false};
            NeighborSearcher *m_neighborSearcher{nullptr};
            NeighborSearcher *m_bNeighborSearcher{nullptr};
            float m_particleRadius{0};
            float3 m_sceneLB;
            float3 m_sceneSize;
            bool m_rIsInit{false};
            bool m_sIsInit{false};
            bool m_isStart{true};
            double m_mem{0};
            uint32_t m_blockNum;
            uint32_t m_threadNum;

        private:
            ConstParams *m_device_cp{nullptr};
            DynamicParams *m_device_dp{nullptr};
        };

    }
}


#endif //SOSIM_MSMSPH_SOLVER_HPP
