//@author        : Long Shen
//@date          : 2023/10/11
//@description   : 
//@version       : 1.0

#ifndef SOSIM_SOLVER_HPP
#define SOSIM_SOLVER_HPP

#include <vector_types.h>
#include <iostream>

#include "Public/Framework/object.hpp"
#include "Public/Framework/framework_config.hpp"

namespace SoSim {
    class Solver {
    public:
        Solver() = default;

        virtual ~Solver() = default;

        virtual void refresh() = 0;

        virtual void runTimeRange(double t) = 0;

        virtual void runSingleStep() = 0;

        virtual void destroy() = 0;

        virtual void attachObject(Object *obj) = 0;

        virtual void removeObject(Object *obj) = 0;

        virtual SolverConfig *getConfig() = 0;
    };
}

#endif //SOSIM_SOLVER_HPP
