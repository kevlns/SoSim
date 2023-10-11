//@author        : Long Shen
//@date          : 2023/10/11
//@description   : 
//@version       : 1.0

#ifndef SOSIM_DATA_PACK_HPP
#define SOSIM_DATA_PACK_HPP

namespace SoSim::NSUGB {

    struct ConstParams {
        unsigned maxNeighborNum;
        unsigned totalParticleNum;
        unsigned cellNum;
        float cellLength;
//        float3 sceneLB;
//        float3 sceneSize;
//        uint3 gridSize;
//        int3 cellOffsets[];
    };

    struct DynamicParams {

    };

}

#endif //SOSIM_DATA_PACK_HPP
