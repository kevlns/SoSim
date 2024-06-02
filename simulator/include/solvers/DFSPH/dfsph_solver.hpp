//
// Created by ADMIN on 2024/3/7.
//

#ifndef SOSIM_DFSPH_SOLVER_HPP
#define SOSIM_DFSPH_SOLVER_HPP

#include <optional>
#include <set>

#include "framework/MetaFramework/solver.hpp"
#include "framework/object_manager.hpp"
#include "libs/NeighborSearchL/unified_grid_ns.hpp"
#include "dfsph_parameters.hpp"

namespace SoSim {

    struct DFSPHSolverConfig : public SolverConfig {
        float dt{1 / 60.f};
        float cur_sim_time{0};
        Vec3f gravity{0, -9.8, 0};
        float rest_vis{0.001};
        float rest_density{1000};
        float rest_rigid_density{1200};
        Vec3f scene_lb{-15, -15, -15};
        Vec3f scene_size{30, 30, 30};
        unsigned max_neighborNum{35};
        unsigned kernel_blocks{0};
        unsigned kernel_threads{0};

        unsigned particle_num_copy{0};
        bool export_data{false};
    };

    class DFSPHSolver : public Solver {
    public:
        DFSPHSolver();

        ~DFSPHSolver() override;

        std::shared_ptr<SolverConfig> getConfig() override;

        void attachObject(std::shared_ptr<Object> object) override;

        void attachParticleEmitter(std::shared_ptr<ParticleEmitter> emitter) override {};

        void detachObject(std::shared_ptr<Object> object) override;

        void detachParticleEmitter(std::shared_ptr<ParticleEmitter> emitter) override {};

        bool initialize() override;


        void run(float total_time) override;

    protected:
        void step() override;

    private:
        void mergeObjects();

        void deviceMalloc();

        void destroy();

        void exportData(DFSPHSolverConfig *solver_config);

    private:
        bool m_change_occur{false}; // for re-config solver
        bool m_is_init{false};
        bool m_is_start{true};
        std::shared_ptr<SolverConfig> m_config;
        std::set<std::shared_ptr<Object>> m_objects;
        DFSPHConstantParams m_host_const;
        DFSPHConstantParams *m_device_const{nullptr};
        DFSPHDynamicParams m_host_data;
        DFSPHDynamicParams *m_device_data{nullptr};
        NeighborSearchUG m_neighborSearch;
    };

}

#endif //SOSIM_DFSPH_SOLVER_HPP
