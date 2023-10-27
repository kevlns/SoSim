//@author        : Long Shen
//@date          : 2023/10/11
//@description   : 
//@version       : 1.0

#ifndef SOSIM_SOLVER_HPP
#define SOSIM_SOLVER_HPP

#include <vector_types.h>

namespace SoSim {

    class Solver {
    public:
        struct SolverConfig {
            double dt{0};
            double cur_time{0};
            float3 gravity{0.f, -9.8f, 0.0};
        };

    public:
        Solver() = default;

        virtual ~Solver() = default;

        virtual void initialize() = 0;

        virtual bool isInitialized() const = 0;

        virtual void run() = 0;

        virtual void destroy() = 0;

        virtual void setConfig(const SolverConfig &config) = 0;

        virtual void addObject() = 0;

    protected:
        virtual void step() = 0;
    };

}

#endif //SOSIM_SOLVER_HPP
