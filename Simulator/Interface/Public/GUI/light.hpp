//@author        : Long Shen
//@date          : 2023/11/21
//@description   :
//@version       : 1.0

#ifndef SOSIM_LIGHT_HPP
#define SOSIM_LIGHT_HPP

#include <iostream>

#include "glm/glm.hpp"

namespace SoSim {

    struct PointLight {
        glm::vec3 pos;
        glm::vec4 color;
        float intensity;
    };

    struct AreaLight {

    };

    struct EnvLight {
        std::string map_path;
    };

}

#endif //SOSIM_LIGHT_HPP
