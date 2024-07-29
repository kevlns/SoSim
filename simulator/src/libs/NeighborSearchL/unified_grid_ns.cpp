//
// Created by ADMIN on 2024/3/5.
//
#include <cuda_runtime.h>
#include <iomanip>
#include <string>

#include "libs/NeighborSearchL/unified_grid_ns.hpp"
#include "unified_grid_ns_kernel_api.cuh"

namespace SoSim {

    void NeighborSearchUGParams::malloc(const NeighborSearchUGConfig &config) {
        // reset
        freeMemory();

        // init neighbor search
        auto size1 = config.particle_num * sizeof(unsigned);
        auto size2 = config.cellNum * sizeof(unsigned);
        auto size3 = config.particle_num * config.maxNeighborNum * sizeof(unsigned);
        auto size4 = 27 * sizeof(Vec3i);
        cudaMalloc((void **) &particleIndices_cuData, size1);
        cudaMalloc((void **) &cellIndices_cuData, size1);
        cudaMalloc((void **) &neighborNum_cuData, size1);
        cudaMalloc((void **) &cellStart_cuData, size2);
        cudaMalloc((void **) &cellEnd_cuData, size2);
        cudaMalloc((void **) &neighbors_cuData, size3);
        cudaMalloc((void **) &cellOffsets_cuData, size4);

        std::vector<Vec3i> offsets = {
                {-1, -1, -1},
                {-1, -1, 0},
                {-1, -1, 1},
                {-1, 0,  -1},
                {-1, 0,  0},
                {-1, 0,  1},
                {-1, 1,  -1},
                {-1, 1,  0},
                {-1, 1,  1},
                {0,  -1, -1},
                {0,  -1, 0},
                {0,  -1, 1},
                {0,  0,  -1},
                {0,  0,  0},
                {0,  0,  1},
                {0,  1,  -1},
                {0,  1,  0},
                {0,  1,  1},
                {1,  -1, -1},
                {1,  -1, 0},
                {1,  -1, 1},
                {1,  0,  -1},
                {1,  0,  0},
                {1,  0,  1},
                {1,  1,  -1},
                {1,  1,  0},
                {1,  1,  1}};
        cudaMemcpy(cellOffsets_cuData, offsets.data(), size4, cudaMemcpyHostToDevice);
        cudaGetLastError();

        isInit = true;
    }

    void NeighborSearchUGParams::freeMemory() {
        if (isInit) {
            cudaFree(cellOffsets_cuData);
            cudaFree(particleIndices_cuData);
            cudaFree(cellIndices_cuData);
            cudaFree(cellStart_cuData);
            cudaFree(cellEnd_cuData);
            cudaFree(neighborNum_cuData);
            cudaFree(neighbors_cuData);

            cudaGetLastError();

            isInit = false;
        }
    }

    void NeighborSearchUG::setConfig(NeighborSearchUGConfig config) {
        h_config = config;
        h_config.gridSize = {
                static_cast<uint32_t>(std::ceil(h_config.sceneSize.x / h_config.cellLength)),
                static_cast<uint32_t>(std::ceil(h_config.sceneSize.y / h_config.cellLength)),
                static_cast<uint32_t>(std::ceil(h_config.sceneSize.z / h_config.cellLength))};
        h_config.cellNum = h_config.gridSize.x * h_config.gridSize.y * h_config.gridSize.z;
        if (d_config)
            cudaFree(d_config);

        cudaMalloc((void **) &d_config, sizeof(NeighborSearchUGConfig));
        cudaMemcpy(d_config, &h_config, sizeof(NeighborSearchUGConfig), cudaMemcpyHostToDevice);
        cudaGetLastError();
    }

    void NeighborSearchUG::malloc() {
        if (d_params)
            cudaFree(d_params);

        h_params.malloc(h_config);
        cudaMalloc((void **) &d_params, sizeof(NeighborSearchUGParams));
        cudaMemcpy(d_params, &h_params, sizeof(NeighborSearchUGParams), cudaMemcpyHostToDevice);
        cudaGetLastError();
    }

    void NeighborSearchUG::update(Vec3f *pos_cuData) {
        cu_update(h_config, h_params, d_config, d_params, pos_cuData);
    }

    void NeighborSearchUG::dump() const {
        auto cellNum = h_config.cellNum;
        auto particleNum = h_config.particle_num;
        auto *c_cellStart = new uint32_t[cellNum];
        auto *c_cellEnd = new uint32_t[cellNum];
        cudaMemcpy(c_cellStart, h_params.cellStart_cuData, cellNum * sizeof(uint32_t), cudaMemcpyDeviceToHost);
        cudaMemcpy(c_cellEnd, h_params.cellEnd_cuData, cellNum * sizeof(uint32_t), cudaMemcpyDeviceToHost);
        uint32_t cnt = 0;
        for (int i = 0; i < cellNum; ++i) {
            if (c_cellStart[i] != UINT_MAX) {
                cnt++;
            }
        }
        delete[] c_cellStart;
        delete[] c_cellEnd;

        std::cout << std::setw(25) << "UGBNeighborSearcher::" << std::setw(35) << "Particle Num: " << std::setw(20)
                  << h_config.particle_num << "\n";
        std::cout << std::setw(25) << "UGBNeighborSearcher::" << std::setw(35) << "Cell Num: " << std::setw(20)
                  << cellNum
                  << "\n";
        std::cout << std::setw(25) << "UGBNeighborSearcher::" << std::setw(35) << "Grid Size: " << std::setw(20)
                  << std::to_string(h_config.gridSize.x) + " * " + std::to_string(h_config.gridSize.y)
                     + " * " + std::to_string(h_config.gridSize.z) << "\n";
//        std::cout << std::setw(25) << "UGBNeighborSearcher::" << std::setw(35) << "Allocated Mem: " << std::setw(20)
//                  << std::to_string(m_mem) + " MB" << "\n";
        std::cout << std::setw(25) << "UGBNeighborSearcher::" << std::setw(35) << "Average PartNum per cell: "
                  << std::setw(20)
                  << (particleNum / cnt) << "\n";
        std::cout << std::setw(25) << "UGBNeighborSearcher::" << std::setw(35) << "Activate Cell num: " << std::setw(20)
                  << cnt << "\n\n";
    }

    void NeighborSearchUG::freeMemory() {
        if (d_config)
            cudaFree(d_config);

        h_params.freeMemory();
        if (d_params)
            cudaFree(d_params);

        cudaGetLastError();
    }
}