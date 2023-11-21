//@author        : Long Shen
//@date          : 2023/10/11
//@description   :
//@version       : 1.0

#ifndef SOSIM_NSUGB_DATA_PACK_HPP
#define SOSIM_NSUGB_DATA_PACK_HPP

#include <cuda_runtime.h>
#include <vector_types.h>
#include <vector>

#include "Public/Shared/CudaUtils/cuda_tool.hpp"

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

        void destroy() const {
            cudaFree(particleIndices);
            cudaFree(cellIndices);
            cudaFree(cellStart);
            cudaFree(cellEnd);
            cudaFree(neighborNum);
            cudaFree(neighbors);

            cudaGetLastError_t("ERROR::NeighborSearcher::DynamicParams destroy() failed.");
        }
    };

}  // namespace SoSim::NSUGB

#endif  // SOSIM_NSUGB_DATA_PACK_HPP
