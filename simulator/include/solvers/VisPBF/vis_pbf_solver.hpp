//
// Created by ADMIN on 2024/6/13.
//

#ifndef SOSIM_VIS_PBF_SOLVER_HPP
#define SOSIM_VIS_PBF_SOLVER_HPP

#include <optional>
#include <set>

#include "framework/MetaFramework/solver.hpp"
#include "framework/object_manager.hpp"
#include "libs/NeighborSearchL/unified_grid_ns.hpp"
#include "vis_pbf_parameters.hpp"

namespace SoSim {

    struct VisPBFSolverConfig : public SolverConfig {
        // common
        float dt{1 / 60.f};
        float cur_sim_time{0};
        Vec3f gravity{0, -9.8, 0};
        unsigned kernel_threads{0};
        unsigned kernel_blocks{0};
        Vec3f scene_lb{-15, -15, -15};
        Vec3f scene_size{30, 30, 30};
        unsigned max_neighborNum{60};

        // custom
        float rest_density;
        float rest_rigid_density{1000};
        float rest_bound_density{1000};
        unsigned pbf_iter_num{3};

        // export setting
        bool export_data{false};
        std::optional<std::string> export_path;
        float export_fps{30};
        std::string export_partial;

    };

    class VisPBFSolver : public Solver {
    public:
        VisPBFSolver();

        ~VisPBFSolver() override;

        std::shared_ptr<SolverConfig> getConfig() override;

        void attachObject(std::shared_ptr<Object> object) override;

        void attachParticleEmitter(std::shared_ptr<ParticleEmitter> emitter) override;

        void detachObject(std::shared_ptr<Object> object) override;

        void detachParticleEmitter(std::shared_ptr<ParticleEmitter> emitter) override;

        bool initialize() override;
        
        void run(float total_time) override;

    protected:
        void step() override;

    private:
        void mergeObjects();

        void destroy();

        void exportAsPly();

    private:
        bool m_change_occur{false}; // for re-config solver
        bool m_is_init{false};
        bool m_is_start{true};
        bool m_is_crash{false};
        std::shared_ptr<SolverConfig> m_config;
        std::set<std::shared_ptr<Object>> m_objects;
        std::set<std::shared_ptr<ParticleEmitter>> m_emitters;
        VisPBFConstantParams m_host_const;
        VisPBFConstantParams *m_device_const{nullptr};
        VisPBFDynamicParams m_host_data;
        VisPBFDynamicParams *m_device_data{nullptr};
        NeighborSearchUG m_neighborSearch;
        Vec3ui m_unified_part_type_start_index; // start index of three sim-mat: ISMCT_NONNEWTON, DYNAIMC_RIGID, FIXED_BOUND
        std::vector<unsigned> m_obj_start_index; // start index of each obj

        // custom
        std::vector<Vec3f> pos_all;
        std::vector<Vec3f> vel_all;
        std::vector<Material> mat_all;
    };

}

#endif //SOSIM_VIS_PBF_SOLVER_HPP
