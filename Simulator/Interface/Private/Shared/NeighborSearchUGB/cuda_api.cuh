//@author        : Long Shen
//@date          : 2023/10/27
//@description   :
//@version       : 1.0

#ifndef SOSIM_CUDA_API_CUH
#define SOSIM_CUDA_API_CUH

#include "Private/Shared/NeighborSearchUGB/data_pack.hpp"
#include "Private/BuildSystem/macro_helper.hpp"

namespace SoSim::NSUGB {

    SOSIM_API __device__ inline int3 floor_to_int3(const float3 &v);

    SOSIM_API __device__ inline int3 getCellPos(const float3 &pos, const float3 &sceneLB, float cellLength);

    SOSIM_API __device__ inline uint32_t getCellId(const int3 &cellPos, const uint3 &gridSize);

    SOSIM_API __device__ inline bool cellIsAvailable(const int3 &cellPos, const uint3 &gridSize);

    SOSIM_API __device__ inline bool cellIsActivated(uint32_t cellId, const uint32_t *cellStart);

    SOSIM_API __host__ void
    ns_resetDevPtr(ConstParams &h_cp, DynamicParams &h_dp);

    SOSIM_API __global__ void
    ns_calcParticleHashValue(ConstParams *d_cp, DynamicParams *d_dp, float3 *d_pos);

    SOSIM_API __host__ void
    ns_sortByHashValue(ConstParams &h_cp, DynamicParams &h_dp);

    SOSIM_API __global__ void
    ns_findCellRange(ConstParams *d_cp, DynamicParams *d_dp);

    SOSIM_API __global__ void
    ns_findNeighbors(ConstParams *d_cp, DynamicParams *d_dp, float3 *d_pos);

    SOSIM_API __host__ void
    update_cuda(ConstParams &h_cp, DynamicParams &h_dp, ConstParams *d_cp, DynamicParams *d_dp, uint32_t blockNum,
                uint32_t threadNum, float3 *d_pos);

}

#endif //SOSIM_CUDA_API_CUH
