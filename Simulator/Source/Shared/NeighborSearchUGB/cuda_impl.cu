//@author        : Long Shen
//@date          : 2023/10/27
//@description   :
//@version       : 1.0

#include <thrust/device_vector.h>
#include <thrust/sort.h>

#include "Private/Shared/NeighborSearchUGB/cuda_api.cuh"
#include "Public/Shared/Math/helper_math.hpp"

namespace SoSim::NSUGB {

    __device__ inline int3
    floor_to_int3(const float3 &v) {
        return make_int3(static_cast<int>(floor(v.x)),
                         static_cast<int>(floor(v.y)),
                         static_cast<int>(floor(v.z)));
    }

    __device__ inline int3
    getCellPos(const float3 &pos, const float3 &sceneLB, float cellLength) {
        int3 cellPos = floor_to_int3((pos - sceneLB) / cellLength);
        return cellPos;
    }

    __device__ inline uint32_t
    getCellId(const int3 &cellPos, const uint3 &gridSize) {
        uint32_t cellId =
                cellPos.z * (gridSize.y * gridSize.x) + cellPos.y * gridSize.x + cellPos.x;
        return cellId;
    }

    __device__ inline bool
    cellIsAvailable(const int3 &cellPos, const uint3 &gridSize) {
        auto cellStartPos = make_int3(0, 0, 0);
        auto cellEndPos = make_int3(gridSize);

        return (cellPos.x >= cellStartPos.x && cellPos.y >= cellStartPos.y && cellPos.z >= cellStartPos.z &&
                cellPos.x <= cellEndPos.x && cellPos.y <= cellEndPos.y && cellPos.z <= cellEndPos.z);
    }

    __device__ inline bool
    cellIsActivated(uint32_t cellId, const uint32_t *cellStart) {
        return (cellStart[cellId] != UINT_MAX);
    }

    __host__ void
    ns_resetDevPtr(ConstParams &h_cp, DynamicParams &h_dp) {
        static size_t size1 = h_cp.cellNum * sizeof(uint32_t);
        static size_t size2 = h_cp.totalParticleNum * sizeof(uint32_t);
        static size_t size3 = h_cp.totalParticleNum * h_cp.maxNeighborNum * sizeof(uint32_t);

        cudaMemset(h_dp.cellStart, UINT_MAX, size1);
        cudaMemset(h_dp.cellEnd, UINT_MAX, size1);
        cudaMemset(h_dp.neighborNum, 0, size2);
        cudaMemset(h_dp.neighbors, UINT_MAX, size3);
    }

    __global__ void
    ns_calcParticleHashValue(ConstParams *d_cp, DynamicParams *d_dp, float3 *d_pos) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_cp->totalParticleNum)
            return;

        auto cellPos = getCellPos(d_pos[i], d_cp->sceneLB, d_cp->cellLength);
        if (cellIsAvailable(cellPos, d_cp->gridSize)) {
            uint32_t cellId = getCellId(cellPos, d_cp->gridSize);
            d_dp->cellIndices[i] = cellId;
            d_dp->particleIndices[i] = i;
        }
    }

    __host__ void
    ns_sortByHashValue(ConstParams &h_cp, DynamicParams &h_dp) {
        thrust::device_ptr<uint32_t> keys_dev_ptr(h_dp.cellIndices);
        thrust::device_ptr<uint32_t> values_dev_ptr(h_dp.particleIndices);

        // 使用thrust::sort_by_key函数根据键进行排序
        thrust::sort_by_key(keys_dev_ptr, keys_dev_ptr + h_cp.totalParticleNum, values_dev_ptr);
    }

    __global__ void
    ns_findCellRange(ConstParams *d_cp, DynamicParams *d_dp) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        uint32_t pre_i = i - 1;
        if (i >= d_cp->totalParticleNum)
            return;

        uint32_t curCellId = d_dp->cellIndices[i];
        if (i == 0)
            d_dp->cellStart[curCellId] = 0;
        else {
            uint32_t preCellId = d_dp->cellIndices[pre_i];
            if (curCellId != preCellId) {
                d_dp->cellStart[curCellId] = i;
                d_dp->cellEnd[preCellId] = pre_i;
            }

            if (i == d_cp->totalParticleNum - 1)
                d_dp->cellEnd[curCellId] = i;
        }
    }

    __global__ void
    ns_findNeighbors(ConstParams *d_cp, DynamicParams *d_dp, float3 *d_pos) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_cp->totalParticleNum)
            return;

        auto p_i = d_dp->particleIndices[i];
        auto pos_i = d_pos[p_i];
        auto pn_index = p_i * d_cp->maxNeighborNum;
        int3 curCellPos = getCellPos(d_pos[p_i], d_cp->sceneLB, d_cp->cellLength);
        for (auto offset: d_cp->cellOffsets) {
            int3 cellPos = curCellPos + offset;
            auto cellId = getCellId(cellPos, d_cp->gridSize);
            if (cellIsAvailable(cellPos, d_cp->gridSize) && cellIsActivated(cellId, d_dp->cellStart)) {
                for (uint32_t j = d_dp->cellStart[cellId]; j <= d_dp->cellEnd[cellId]; ++j) {
                    auto p_j = d_dp->particleIndices[j];
                    auto pos_j = d_pos[p_j];
                    if (length(pos_i - pos_j) <= d_cp->cellLength) {
                        if (d_dp->neighborNum[p_i] < d_cp->maxNeighborNum) {
                            auto ind_offset = d_dp->neighborNum[p_i]++;
                            d_dp->neighbors[pn_index + ind_offset] = p_j;
                        }
                    }
                }
            }
        }
    }

    __host__ void
    update_cuda(ConstParams &h_cp, DynamicParams &h_dp, ConstParams *d_cp, DynamicParams *d_dp, uint32_t blockNum,
                uint32_t threadNum, float3 *d_pos) {

        ns_resetDevPtr(h_cp, h_dp);

        ns_calcParticleHashValue<<< blockNum, threadNum >>>(d_cp, d_dp, d_pos);

        ns_sortByHashValue(h_cp, h_dp);

        ns_findCellRange<<< blockNum, threadNum >>>(d_cp, d_dp);

        ns_findNeighbors<<< blockNum, threadNum >>>(d_cp, d_dp, d_pos);

    }

}