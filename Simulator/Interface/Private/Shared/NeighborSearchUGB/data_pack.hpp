//@author        : Long Shen
//@date          : 2023/10/11
//@description   :
//@version       : 1.0

#ifndef SOSIM_NSUGB_DATA_PACK_HPP
#define SOSIM_NSUGB_DATA_PACK_HPP

#include <vector_types.h>
#include <vector>

namespace SoSim::NSUGB {

    struct ConstParams {
        uint32_t maxNeighborNum;
        uint32_t totalParticleNum;
        uint32_t cellNum;
        float cellLength;
        float3 sceneLB;
        float3 sceneSize;
        uint3 gridSize;
        int3 cellOffsets[27];
    };

    struct DynamicParams {
        uint32_t *particleIndices;
        uint32_t *cellIndices;
        uint32_t *cellStart;
        uint32_t *cellEnd;
        uint32_t *neighborNum;
        uint32_t *neighbors;
    };

}  // namespace SoSim::NSUGB

#endif  // SOSIM_NSUGB_DATA_PACK_HPP
