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

    class Solver {
    public:
        virtual ~Solver() = default;

        virtual void attachObject(Object *object) = 0;

        virtual void detachObject(Object *object) = 0;

        virtual bool initialize() = 0;

        virtual void destroy() = 0;

        virtual void run(float total_time) = 0;

    protected:
        virtual void step() = 0;

    private:
        std::set<Object *> m_objects;
    };


}

#endif //SOSIM_SOLVER_HPP
