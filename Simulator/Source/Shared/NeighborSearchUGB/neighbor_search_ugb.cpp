//@author        : Long Shen
//@date          : 2023/10/11
//@description   : 
//@version       : 1.0

#include <cuda_runtime.h>
#include <iostream>
#include <iomanip>
#include <string>

#include "Public/Shared/NeighborSearchUGB/neighbor_search_ugb.hpp"
#include "Public/Shared/CudaUtils/cuda_tool.hpp"
#include "Private/Shared/NeighborSearchUGB/cuda_api.cuh"

namespace SoSim::NSUGB {

    void NeighborSearchUGB::initialize(float3 scene_lb, float3 scene_size, unsigned int total_particle_num,
                                       float sph_support_radius) {
        m_isInit = true;

        m_host_cp.maxNeighborNum = 35;
        m_host_cp.totalParticleNum = total_particle_num;
        m_host_cp.sceneLB = scene_lb;
        m_host_cp.sceneSize = scene_size;
        m_host_cp.cellLength = sph_support_radius;

        m_host_cp.gridSize = make_uint3(
                static_cast<uint32_t>(std::ceil(m_host_cp.sceneSize.x / m_host_cp.cellLength)),
                static_cast<uint32_t>(std::ceil(m_host_cp.sceneSize.y / m_host_cp.cellLength)),
                static_cast<uint32_t>(std::ceil(m_host_cp.sceneSize.z / m_host_cp.cellLength)));
        m_host_cp.cellNum = m_host_cp.gridSize.x * m_host_cp.gridSize.y * m_host_cp.gridSize.z;

        int3 _cellOffsets_[27] = {
                make_int3(-1, -1, -1),
                make_int3(-1, -1, 0),
                make_int3(-1, -1, 1),
                make_int3(-1, 0, -1),
                make_int3(-1, 0, 0),
                make_int3(-1, 0, 1),
                make_int3(-1, 1, -1),
                make_int3(-1, 1, 0),
                make_int3(-1, 1, 1),
                make_int3(0, -1, -1),
                make_int3(0, -1, 0),
                make_int3(0, -1, 1),
                make_int3(0, 0, -1),
                make_int3(0, 0, 0),
                make_int3(0, 0, 1),
                make_int3(0, 1, -1),
                make_int3(0, 1, 0),
                make_int3(0, 1, 1),
                make_int3(1, -1, -1),
                make_int3(1, -1, 0),
                make_int3(1, -1, 1),
                make_int3(1, 0, -1),
                make_int3(1, 0, 0),
                make_int3(1, 0, 1),
                make_int3(1, 1, -1),
                make_int3(1, 1, 0),
                make_int3(1, 1, 1),
        };

        memcpy(m_host_cp.cellOffsets, _cellOffsets_, 27 * sizeof(int3));

        cudaMalloc_t((void **) &m_device_cp, sizeof(ConstParams), m_mem);
        cudaMemcpy(m_device_cp, &m_host_cp, sizeof(ConstParams), cudaMemcpyHostToDevice);
        m_isInit &= cudaGetLastError_t("NeighborSearchUGB::initialize():: init m_device_cp failed!");

        size_t size1 = m_host_cp.totalParticleNum * sizeof(uint32_t);
        size_t size2 = m_host_cp.cellNum * sizeof(uint32_t);
        size_t size3 = m_host_cp.totalParticleNum * m_host_cp.maxNeighborNum * sizeof(uint32_t);

        cudaMalloc_t((void **) &m_host_dp.particleIndices, size1, m_mem);
        cudaMalloc_t((void **) &m_host_dp.cellIndices, size1, m_mem);
        cudaMalloc_t((void **) &m_host_dp.neighborNum, size1, m_mem);
        cudaMalloc_t((void **) &m_host_dp.cellStart, size2, m_mem);
        cudaMalloc_t((void **) &m_host_dp.cellEnd, size2, m_mem);
        cudaMalloc_t((void **) &m_host_dp.neighbors, size3, m_mem);

        cudaMalloc_t((void **) &m_device_dp, sizeof(DynamicParams), m_mem);
        cudaMemcpy(m_device_dp, &m_host_dp, sizeof(DynamicParams), cudaMemcpyHostToDevice);

        int device;
        cudaGetDevice(&device);

        cudaDeviceProp prop{};
        cudaGetDeviceProperties(&prop, device);

        m_threadNum = prop.maxThreadsPerBlock;
        m_blockNum = (m_host_cp.totalParticleNum + m_threadNum - 1) / m_threadNum;

        m_isInit &= cudaGetLastError_t("NeighborSearchUGB::initialize():: init m_device_dp failed!");
    }

    void NeighborSearchUGB::update(float3 *device_pos) {
        update_cuda(m_host_cp, m_host_dp, m_device_cp, m_device_dp, m_blockNum, m_threadNum, device_pos);
    }

    void NeighborSearchUGB::destroy() {
        if (m_isInit) {
            cudaFree(m_device_cp);

            cudaFree(m_host_dp.particleIndices);
            cudaFree(m_host_dp.cellIndices);
            cudaFree(m_host_dp.cellStart);
            cudaFree(m_host_dp.cellStart);
            cudaFree(m_host_dp.neighborNum);
            cudaFree(m_host_dp.neighbors);

            cudaFree(m_device_dp);

            m_device_cp = nullptr;
            m_device_dp = nullptr;
        }
    }

    uint32_t *NeighborSearchUGB::getPartIndexDevicePtr() const {
        return m_device_dp->particleIndices;
    }

    uint32_t *NeighborSearchUGB::getNeighborsDevicePtr() const {
        return m_device_dp->neighbors;
    }

    void NeighborSearchUGB::dumpInfo() const {
        auto cellNum = m_host_cp.cellNum;
        auto particleNum = m_host_cp.totalParticleNum;
        auto *c_cellStart = new uint32_t[cellNum];
        auto *c_cellEnd = new uint32_t[cellNum];
        cudaMemcpy(c_cellStart, m_host_dp.cellStart, cellNum * sizeof(uint32_t), cudaMemcpyDeviceToHost);
        cudaMemcpy(c_cellEnd, m_host_dp.cellEnd, cellNum * sizeof(uint32_t), cudaMemcpyDeviceToHost);
        uint32_t cnt = 0;
        for (int i = 0; i < cellNum; ++i) {
            if (c_cellStart[i] != UINT_MAX) {
                cnt++;
            }
        }
        delete[] c_cellStart;
        delete[] c_cellEnd;

        std::cout << std::setw(25) << "UGBNeighborSearcher::" << std::setw(35) << "Particle Num: " << std::setw(20)
                  << m_host_cp.totalParticleNum << "\n";
        std::cout << std::setw(25) << "UGBNeighborSearcher::" << std::setw(35) << "Cell Num: " << std::setw(20)
                  << cellNum
                  << "\n";
        std::cout << std::setw(25) << "UGBNeighborSearcher::" << std::setw(35) << "Grid Size: " << std::setw(20)
                  << std::to_string(m_host_cp.gridSize.x) + " * " + std::to_string(m_host_cp.gridSize.y)
                     + " * " + std::to_string(m_host_cp.gridSize.z) << "\n";
        std::cout << std::setw(25) << "UGBNeighborSearcher::" << std::setw(35) << "Allocated Mem: " << std::setw(20)
                  << std::to_string(m_mem) + " MB" << "\n";
        std::cout << std::setw(25) << "UGBNeighborSearcher::" << std::setw(35) << "Average PartNum per cell: "
                  << std::setw(20)
                  << (particleNum / cnt) << "\n";
        std::cout << std::setw(25) << "UGBNeighborSearcher::" << std::setw(35) << "Activate Cell num: " << std::setw(20)
                  << cnt << "\n\n";
    }

}