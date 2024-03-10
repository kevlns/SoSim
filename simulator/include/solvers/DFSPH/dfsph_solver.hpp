//
// Created by ADMIN on 2024/3/7.
//

#ifndef SOSIM_DFSPH_SOLVER_HPP
#define SOSIM_DFSPH_SOLVER_HPP

#include <optional>
#include <set>

#include "framework/solver.hpp"
#include "framework/object_manager.hpp"
#include "libs/NeighborSearchL/unified_grid_ns.hpp"
#include "dfsph_parameters.hpp"

namespace SoSim {

    struct DFSPHSolverConfig {
        float dt{1 / 60.f};
        float cur_sim_time{0};
        Vec3f gravity{0, -9.8, 0};
        float rest_vis{0.001};
        float rest_density{1000};
        Vec3f scene_lb{-15, -15, -15};
        Vec3f scene_size{30, 30, 30};
        unsigned max_neighborNum{35};
        unsigned kernel_blocks{0};
        unsigned kernel_threads{0};
    };

    class DFSPHSolver : public Solver {
    public:
        ~DFSPHSolver() override = default;

        void setConfig(DFSPHSolverConfig *config);

        DFSPHSolverConfig *getConfig();

        void attachObject(Object *object) override;

        void detachObject(Object *object) override;

        bool initialize() override;

        void destroy() override;

        void run(float total_time) override;

    protected:
        void step() override;

    private:
        void mergeObjects();

        void deviceMalloc();

    private:
        bool m_changeOccur{false}; // for re-config solver
        bool m_isInit{false};
        std::optional<DFSPHSolverConfig *> m_config;
        std::set<Object *> m_objects;
        DFSPHConstantParams m_host_const;
        DFSPHConstantParams *m_device_const;
        DFSPHDynamicParams m_host_data;
        DFSPHDynamicParams *m_device_data;
        NeighborSearchUG m_neighborSearch;
        std::vector<Vec3f> m_host_pos_reorder;
    };

}

#endif //SOSIM_DFSPH_SOLVER_HPP
