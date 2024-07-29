//
// Created by ADMIN on 2024/3/7.
//

#ifndef SOSIM_IMM_SOLVER_v2_HPP
#define SOSIM_IMM_SOLVER_v2_HPP

#include <optional>
#include <set>

#include "framework/MetaFramework/solver.hpp"
#include "framework/object_manager.hpp"
#include "libs/NeighborSearchL/unified_grid_ns.hpp"
#include "imm_parameters_v2.hpp"

namespace SoSim {

    struct IMMSolverConfig_v2 : public SolverConfig {
        // common
        float dt{1 / 60.f};
        float cur_sim_time{0};
        Vec3f gravity{0, -9.8, 0};
        unsigned kernel_threads{0};
        unsigned kernel_blocks{0};

        // export setting
        bool export_data{false};
        std::optional<std::string> export_path;
        float export_fps{30};
        std::string export_partial;
        bool export_phase{false};

        // custom
        Vec3f scene_lb{-15, -15, -15};
        Vec3f scene_size{30, 30, 30};
        unsigned max_neighborNum{60};
        float rest_viscosity{0.001};
        float rest_rigid_density{1000.};
        float rest_bound_density{1000.};
        float div_free_threshold;
        float incompressible_threshold;

        // imm
        float Cf;
        float Cd;
        std::vector<float> phase_rest_density; // only support 9 phases
        std::vector<Vec3f> phase_color;
    };

    class IMMSolver_v2 : public Solver {
    public:
        IMMSolver_v2();

        ~IMMSolver_v2() override;

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
        IMMConstantParams_v2 m_host_const;
        IMMConstantParams_v2 *m_device_const{nullptr};
        IMMDynamicParams_v2 m_host_data;
        IMMDynamicParams_v2 *m_device_data{nullptr};
        NeighborSearchUG m_neighborSearch;
        Vec3ui m_unified_part_type_start_index; // start index of three sim-mat: ISMCT_NONNEWTON, DYNAIMC_RIGID, FIXED_BOUND
        std::vector<unsigned> m_obj_start_index; // start index of each obj

        float *m_device_phase_density{nullptr};
        Vec3f *m_device_phase_color{nullptr};

        // custom
        std::vector<Vec3f> pos_all;
        std::vector<Vec3f> vel_all;
        std::vector<Material> mat_all;
        std::vector<float> vol_frac_all;
    };

}

#endif //SOSIM_IMM_SOLVER_v2_HPP
