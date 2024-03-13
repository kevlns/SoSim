//
// Created by ADMIN on 2024/3/5.
//

#ifndef SOSIM_SOLVER_HPP
#define SOSIM_SOLVER_HPP

#include <optional>
#include <set>

#include "core/data_type.hpp"
#include "framework/object.hpp"

namespace SoSim {

    struct SolverConfig {
        virtual ~SolverConfig() = default;

        float dt;
        float cur_sim_time;
        Vec3f gravity;
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
