//
// Created by ADMIN on 2024/3/22.
//

#ifndef SOSIM_JL21_CT_SOLVER_HPP
#define SOSIM_JL21_CT_SOLVER_HPP

#include <optional>
#include <set>

#include "framework/solver.hpp"
#include "framework/object_manager.hpp"

#include "libs/NeighborSearchL/unified_grid_ns.hpp"
#include "jl21_ct_parameters.hpp"

namespace SoSim {
    struct JL21CTSolverConfig : public SolverConfig {
        float dt{1 / 60.f};
        float cur_sim_time{0};
        Vec3f gravity{0, -9.8, 0};
        float wc_stiff{1000};
        float rest_vis{0.001};
        Vec2f rest_density{1000, 1000};
        Vec3f scene_lb{-15, -15, -15};
        Vec3f scene_size{30, 30, 30};
        unsigned max_neighborNum{35};
        unsigned kernel_blocks{0};
        unsigned kernel_threads{0};
        Vec3f phase1_color;
        Vec3f phase2_color;
        float kd;
        float Cd;
        float relation_time;
        float compliance_0;
        float compliance_inf;
        float alpha_up_threshold;
        float alpha_low_threshold;

        bool export_data{false};
    };

    class JL21CTSolver : public Solver {
    public:
        JL21CTSolver();

        ~JL21CTSolver() override;

        std::shared_ptr<SolverConfig> getConfig() override;

        void attachObject(std::shared_ptr<Object> object) override;

        void detachObject(std::shared_ptr<Object> object) override;

        bool initialize() override;

        void run(float total_time) override;

    protected:
        void step() override;

    private:
        void mergeObjects();

        void deviceMalloc();

        void destroy();

        void exportAsPly();

    private:
        bool m_change_occur{false}; // for re-config solver
        bool m_is_init{false};
        bool m_is_start{true};
        std::shared_ptr<SolverConfig> m_config;
        std::set<std::shared_ptr<Object>> m_objects;
        JL21CTConstantParams m_host_const;
        JL21CTConstantParams *m_device_const{nullptr};
        JL21CTDynamicParams m_host_data;
        JL21CTDynamicParams *m_device_data{nullptr};
        NeighborSearchUG m_neighborSearch;

        // special
        std::vector<Vec3f> pos_all;
        std::vector<Vec3f> vel_all;
        std::vector<Material> mat_all;
        std::vector<Vec2f> alpha_all;
    };
}

#endif //SOSIM_JL21_CT_SOLVER_HPP
