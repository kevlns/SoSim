//
// Created by ADMIN on 2024/3/5.
//
#include "unified_grid_ns_kernel_api.cuh"

#include <thrust/device_vector.h>
#include <thrust/sort.h>

#include "libs/NeighborSearchL/unified_grid_ns.hpp"

namespace SoSim {

    /**
 * @brief math helper, floor [float3] to [int3]
 *
 * @param[in] v 3d-vector of float
 *
 * @return    3d-vector of int
 */
    __device__ inline Vec3i
    floor_to_int3(const Vec3f &v) {
        return {static_cast<int>(floor(v.x)),
                static_cast<int>(floor(v.y)),
                static_cast<int>(floor(v.z))};
    }

/**
 * @brief compute cell pos by particle pos
 *
 * @param[in] pos particle pos
 * @param[in] sceneLB left-bottom of the scene
 * @param[in] cellLength cell length
 *
 * @return    cell pos, 3d-vector of int
 */
    __device__ inline Vec3i
    getCellPos(Vec3f &pos, Vec3f &sceneLB, float cellLength) {
        Vec3i cellPos = floor_to_int3((pos - sceneLB) / cellLength);
        return cellPos;
    }

/**
 * @brief compute cell id by cell pos
 *
 * @param[in] cellPos cell pos
 * @param[in] gridSize size of background grid
 *
 * @return    cell id, uint32_t
 */
    __device__ inline unsigned
    getCellId(const Vec3i &cellPos, const Vec3ui &gridSize) {
        unsigned cellId =
                cellPos.z * (gridSize.y * gridSize.x) + cellPos.y * gridSize.x + cellPos.x;
        return cellId;
    }

/**
 * @brief check if cur cell is available
 *
 * @param[in] cellPos cell pos
 * @param[in] gridSize size of background grid
 *
 * @return    true if available, false otherwise
 */
    __device__ inline bool
    cellIsAvailable(const Vec3i &cellPos, const Vec3ui &gridSize) {
        Vec3i cellStartPos = {0, 0, 0};
        Vec3i cellEndPos = {static_cast<int>(gridSize.x), static_cast<int>(gridSize.y), static_cast<int>(gridSize.z)};

        return (cellPos.x >= cellStartPos.x && cellPos.y >= cellStartPos.y && cellPos.z >= cellStartPos.z &&
                cellPos.x <= cellEndPos.x && cellPos.y <= cellEndPos.y && cellPos.z <= cellEndPos.z);
    }

/**
 * @brief check if cur cell is activated
 *
 * @param[in] cellId cell id
 * @param[in] cellStart device pointer of the cell start array
 *
 * @return    true if cell is not empty, false otherwise
 */
    __device__ inline bool
    cellIsActivated(const unsigned cellId, const unsigned *cellStart) {
        return (cellStart[cellId] != UINT_MAX);
    }

    extern __host__ void
    resetDevPtr(NeighborSearchUGConfig &h_config, NeighborSearchUGParams &h_params) {
        static size_t size1 = h_config.cellNum;
        static size_t size2 = h_config.particle_num;
        static size_t size3 = h_config.particle_num * h_config.maxNeighborNum;

        cudaMemset(h_params.cellStart_cuData, UINT_MAX, size1 * sizeof(uint32_t));
        cudaMemset(h_params.cellEnd_cuData, UINT_MAX, size1 * sizeof(uint32_t));
        cudaMemset(h_params.neighborNum_cuData, 0, size2 * sizeof(uint32_t));
        cudaMemset(h_params.neighbors_cuData, UINT_MAX, size3 * sizeof(uint32_t));
    }

    extern __global__ void
    calcParticleHashValue(NeighborSearchUGConfig *d_config, NeighborSearchUGParams *d_params, Vec3f *pos) {
        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_config->particle_num)
            return;

        auto cellPos = getCellPos(pos[i], d_config->sceneLB, d_config->cellLength);
        if (cellIsAvailable(cellPos, d_config->gridSize)) {
            uint32_t cellId = getCellId(cellPos, d_config->gridSize);
            d_params->cellIndices_cuData[i] = cellId;
            d_params->particleIndices_cuData[i] = i;
        }
    }

    extern __host__ void
    sortByHashValue(NeighborSearchUGConfig &h_config, NeighborSearchUGParams &h_params) {
        thrust::device_ptr<uint32_t> keys_dev_ptr(h_params.cellIndices_cuData);
        thrust::device_ptr<uint32_t> values_dev_ptr(h_params.particleIndices_cuData);

        // use thrust::sort_by_key to order by key
        thrust::sort_by_key(keys_dev_ptr, keys_dev_ptr + h_config.particle_num, values_dev_ptr);
    }

    extern __global__ void
    findCellRange(NeighborSearchUGConfig *d_config, NeighborSearchUGParams *d_params) {
        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
        unsigned pre_i = i - 1;
        if (i >= d_config->particle_num)
            return;

        uint32_t curCellId = d_params->cellIndices_cuData[i];
        if (i == 0)
            d_params->cellStart_cuData[curCellId] = 0;
        else {
            uint32_t preCellId = d_params->cellIndices_cuData[pre_i];
            if (curCellId != preCellId) {
                d_params->cellStart_cuData[curCellId] = i;
                d_params->cellEnd_cuData[preCellId] = pre_i;
            }

            if (i == d_config->particle_num - 1)
                d_params->cellEnd_cuData[curCellId] = i;
        }
    }

    extern __global__ void
    findNeighbors(NeighborSearchUGConfig *d_config, NeighborSearchUGParams *d_params, Vec3f *pos) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_config->particle_num)
            return;

        auto p_i = d_params->particleIndices_cuData[i];
        auto pos_i = pos[p_i];
        auto pn_index = p_i * d_config->maxNeighborNum;
        d_params->neighborNum_cuData[p_i] = 1;
        d_params->neighbors_cuData[pn_index] = p_i;
        Vec3i curCellPos = getCellPos(pos[p_i], d_config->sceneLB, d_config->cellLength);
        for (int t = 0; t < 27; ++t) {
            auto offset = d_params->cellOffsets_cuData[t];
            Vec3i cellPos = curCellPos + offset;
            auto cellId = getCellId(cellPos, d_config->gridSize);
            if (cellIsAvailable(cellPos, d_config->gridSize) && cellIsActivated(cellId, d_params->cellStart_cuData)) {
                for (unsigned j = d_params->cellStart_cuData[cellId]; j <= d_params->cellEnd_cuData[cellId]; ++j) {
                    auto p_j = d_params->particleIndices_cuData[j];
                    if (p_j == p_i)
                        continue;
                    auto pos_j = pos[p_j];
                    if ((pos_i - pos_j).length() > 1e-6 && (pos_i - pos_j).length() <= d_config->cellLength) {
                        if (d_params->neighborNum_cuData[p_i] < d_config->maxNeighborNum) {
                            auto ind_offset = d_params->neighborNum_cuData[p_i]++;
                            d_params->neighbors_cuData[pn_index + ind_offset] = p_j;
                        }
                    }
                }
            }
        }
    }

    extern __host__ void
    cu_update(NeighborSearchUGConfig &h_config, NeighborSearchUGParams &h_params, NeighborSearchUGConfig *d_config,
              NeighborSearchUGParams *d_params, Vec3f *pos) {
        auto block_num = h_config.kernel_blocks;
        auto thread_num = h_config.kernel_threads;

        resetDevPtr(h_config, h_params);

        calcParticleHashValue<<<block_num, thread_num>>>(d_config, d_params, pos);

        sortByHashValue(h_config, h_params);

        findCellRange<<<block_num, thread_num>>>(d_config, d_params);

        findNeighbors<<<block_num, thread_num>>>(d_config, d_params, pos);

        cudaGetLastError();
    }

}