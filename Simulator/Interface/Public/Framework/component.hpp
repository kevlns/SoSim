//@author        : Long Shen
//@date          : 2023/10/27
//@description   :
//@version       : 1.0

#ifndef SOSIM_COMPONENT_HPP
#define SOSIM_COMPONENT_HPP

#include <vector_types.h>
#include <iostream>

namespace SoSim {

    struct Component {
    };

    struct BaseMoveComponent : Component {
        float3 pos;
        float3 vel;
    };

}

#endif //SOSIM_COMPONENT_HPP
