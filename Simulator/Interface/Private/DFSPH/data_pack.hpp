//@author        : Long Shen
//@date          : 2023/10/11
//@description   : 
//@version       : 1.0

#ifndef SOSIM_DFSPH_DATA_PACK_HPP
#define SOSIM_DFSPH_DATA_PACK_HPP

#include <vector_types.h>

namespace SoSim::DFSPH {

    struct ConstParams {
        uint32_t totalParticleNum;
        float3 sceneLB;
        float3 sceneSize;

        float dt;
        float3 gravity;
        float sph_h;
        float cross_vis0;
        float cross_visInf;
        float cross_k;
        float cross_a;
    };

    struct DynamicParams {

    };

}

#endif //SOSIM_DFSPH_DATA_PACK_HPP
