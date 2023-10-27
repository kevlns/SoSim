//@author        : Long Shen
//@date          : 2023/10/11
//@description   : 
//@version       : 1.0

#include <cuda_runtime.h>
#include <iostream>

#include "Public/Shared/NeighborSearchUGB/neighbor_search_ugb.hpp"
#include "Public/Shared/CudaUtils/cuda_tool.hpp"


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

        m_host_cp.cellOffsets.resize(27);
        m_host_cp.cellOffsets = {
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


        cudaMalloc_t((void **) &m_device_cp, sizeof(ConstParams), m_mem);
        cudaMemcpy(m_device_cp, &m_host_cp, sizeof(ConstParams), cudaMemcpyHostToDevice);
        m_isInit &= cudaGetLastError_t("NeighborSearchUGB::initialize():: init m_device_cp failed!");

        m_device_dp = new DynamicParams;
        size_t size1 = m_host_cp.totalParticleNum * sizeof(uint32_t);
        size_t size2 = m_host_cp.cellNum * sizeof(uint32_t);
        size_t size3 = m_host_cp.totalParticleNum * m_host_cp.maxNeighborNum * sizeof(uint32_t);

        cudaMalloc_t((void **) &m_device_dp->particleIndices, size1, m_mem);
        cudaMalloc_t((void **) &m_device_dp->cellIndices, size1, m_mem);
        cudaMalloc_t((void **) &m_device_dp->neighborNum, size1, m_mem);
        cudaMalloc_t((void **) &m_device_dp->cellStart, size2, m_mem);
        cudaMalloc_t((void **) &m_device_dp->cellEnd, size2, m_mem);
        cudaMalloc_t((void **) &m_device_dp->neighbors, size3, m_mem);
        m_isInit &= cudaGetLastError_t("NeighborSearchUGB::initialize():: init m_device_dp failed!");
    }

    void NeighborSearchUGB::update(float3 *device_pos) {

    }

    void NeighborSearchUGB::destroy() {
        if (m_isInit) {
            cudaFree(m_device_cp);
            cudaFree(m_device_dp->particleIndices);
            cudaFree(m_device_dp->cellIndices);
            cudaFree(m_device_dp->cellStart);
            cudaFree(m_device_dp->cellStart);
            cudaFree(m_device_dp->neighborNum);
            cudaFree(m_device_dp->neighbors);

            if (cudaGetLastError_t("NeighborSearchUGB::destroy():: cudaFree failed!")) {
                delete m_device_dp;

                m_device_cp = nullptr;
                m_device_dp = nullptr;

                std::cout << "NeighborSearchUGB destructed.\n";
            }

        }
    }

    uint32_t *NeighborSearchUGB::getPartIndexDevicePtr() const {
        return m_device_dp->particleIndices;
    }

    uint32_t *NeighborSearchUGB::getNeighborsDevicePtr() const {
        return m_device_dp->neighbors;
    }

}