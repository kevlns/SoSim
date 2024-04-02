//
// Created by ADMIN on 2024/3/5.
//

#ifndef SOSIM_SOLVER_HPP
#define SOSIM_SOLVER_HPP

#include <optional>
#include <set>

#include "core/math/matrix.hpp"
#include "object.hpp"

namespace SoSim {

    struct SolverConfig {
        virtual ~SolverConfig() = default;

        // common
        float dt{1 / 60.f};
        float cur_sim_time{0};
        Vec3f gravity{0, -9.8, 0};
        unsigned kernel_threads{0};
        unsigned kernel_blocks{0};
        bool export_data{false};
    };

    class Solver {
    public:
        virtual ~Solver() = default;

        virtual std::shared_ptr<SolverConfig> getConfig() = 0;

        virtual void attachObject(std::shared_ptr<Object> object) = 0;

        virtual void detachObject(std::shared_ptr<Object> object) = 0;

        virtual bool initialize() = 0;

        virtual void run(float total_time) = 0;

    protected:
        virtual void step() = 0;

    private:
        std::set<std::shared_ptr<Object>> m_objects;
    };


}

#endif //SOSIM_SOLVER_HPP
