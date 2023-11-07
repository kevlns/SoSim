//@author        : Long Shen
//@date          : 2023/10/30
//@description   :
//@version       : 1.0

#ifndef SOSIM_MAT_HPP
#define SOSIM_MAT_HPP

namespace SoSim {
    enum Material : uint8_t {
        FLUID,
        RIGID,
        ELASTIC,
        BOUND,
    };

    enum Phase : uint8_t {
        PHASE1,
        PHASE2,
        PHASE3,
    };
}

#endif //SOSIM_MAT_HPP
