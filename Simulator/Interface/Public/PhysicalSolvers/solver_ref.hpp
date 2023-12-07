//@author        : Long Shen
//@date          : 2023/11/12
//@description   : include all solver headers
//@version       : 1.0

#ifndef SOSIM_SOLVER_REF_HPP
#define SOSIM_SOLVER_REF_HPP

namespace SoSim {

    enum SolverInstanceType : uint8_t {
        SOLVER_NONE = 0,

        SAMPLE_SOLVER,
        /**
         *  add your own solver blow
         */

        SOLVER_BOTTOM
    };
    static inline const char *solverItems[] = {"None",
                                               "SampleSolver"};

}

#endif //SOSIM_SOLVER_REF_HPP
