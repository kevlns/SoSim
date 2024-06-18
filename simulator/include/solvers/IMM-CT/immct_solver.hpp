//
// Created by ADMIN on 2024/3/26.
//

#ifndef SOSIM_IMMCT_SOLVER_HPP
#define SOSIM_IMMCT_SOLVER_HPP

#include <optional>
#include <set>

#include "framework/MetaFramework/solver.hpp"
#include "framework/MetaFramework/object.hpp"

#include "libs/NeighborSearchL/unified_grid_ns.hpp"
#include "immct_parameters.hpp"

namespace SoSim {
    struct IMMCTSolverConfig : public SolverConfig {
        // common
        float dt{1 / 60.f};
        float cur_sim_time{0};
        Vec3f gravity{0, -9.8, 0};
        unsigned kernel_threads{0};
        unsigned kernel_blocks{0};
        bool export_data{false};
        float export_fps{35};

        // custom
        Vec3f scene_lb{-15, -15, -15};
        Vec3f scene_size{30, 30, 30};
        unsigned max_neighborNum{35};
        float rest_viscosity{0.001};
        Vec2f rest_density;
        float rest_rigid_density;
        float rest_bound_density;
        float div_free_threshold;
        float incompressible_threshold;
        float Cf;
        float Cd;
        Vec3f phase1_color;
        Vec3f phase2_color;
        float Cd0;
        float ct_thinning_exp0{0};
        float solution_vis_base{0.02};
        float solution_vis_max{8000};
        float ct_relaxation_time{0.1};
        float polymer_vol_frac0{0.6};
        float vis_bound_damp_factor{0.1};

        // export setting
        std::optional<std::string> export_path;
        int export_gap{1};
        std::string export_partial;
        bool export_phase{false};

        // method compare
        float phase1_vis{0.01};
        float phase2_vis{0.01};
    };

    class IMMCTSolver : public Solver {
    public:
        IMMCTSolver();

        ~IMMCTSolver() override;

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

        void syncObjectDeviceJitData();

    private:
        bool m_change_occur{false}; // for re-config solver
        bool m_is_init{false};
        bool m_is_start{true};
        bool m_is_crash{false};
        std::shared_ptr<SolverConfig> m_config;
        std::set<std::shared_ptr<Object>> m_objects;
        std::set<std::shared_ptr<ParticleEmitter>> m_emitters;
        IMMCTConstantParams m_host_const;
        IMMCTConstantParams *m_device_const{nullptr};
        IMMCTDynamicParams m_host_data;
        IMMCTDynamicParams *m_device_data{nullptr};
        NeighborSearchUG m_neighborSearch;
        Vec3ui m_unified_part_type_start_index; // start index of three sim-mat: ISMCT_NONNEWTON, DYNAIMC_RIGID, FIXED_BOUND
        std::vector<unsigned> m_obj_start_index; // start index of each obj

        // custom
        std::vector<Vec3f> pos_all;
        std::vector<Vec3f> vel_all;
        std::vector<Material> mat_all;
        std::vector<Vec2f> vol_frac_all;
    };
}

#endif //SOSIM_IMMCT_SOLVER_HPP
