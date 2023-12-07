//@author        : Long Shen
//@date          : 2023/11/26
//@description   :
//@version       : 1.0

#ifndef MODEL_BUILDER_HPP
#define MODEL_BUILDER_HPP

#include <tuple>
#include <vector>
#include <vector_types.h>

#include "Public/Framework/framework_config.hpp"

namespace SoSim {

    extern std::vector<float3> genFromObjectConfig(const ObjectConfig *objectConfig);

    extern std::vector<float3> genParticleCube(const ObjectConfig *object_config);

    extern std::vector<float3> genParticleBox(const ObjectConfig *object_config);

    extern std::vector<float3> genParticlePlaneX(const ObjectConfig *object_config);

    extern std::vector<float3> genParticlePlaneZ(const ObjectConfig *object_config);

    extern std::vector<float3> genParticleCylinder(const ObjectConfig *object_config);

}

#endif //MODEL_BUILDER_HPP
