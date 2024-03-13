//
// Created by ADMIN on 2024/3/8.
//

#ifndef SOSIM_DFSPH_CUDA_API_CUH
#define SOSIM_DFSPH_CUDA_API_CUH

#include "solvers/DFSPH/dfsph_parameters.hpp"
#include "libs/NeighborSearchL/unified_grid_ns.hpp"

namespace SoSim {

    __device__ inline float cubic_value(const Vec3f &r, float h);

    __device__ inline Vec3f cubic_gradient(const Vec3f &r, float h);

    __global__ void
    init(DFSPHConstantParams *d_const,
         DFSPHDynamicParams *d_data,
         NeighborSearchUGConfig *d_nsConfig,
         NeighborSearchUGParams *d_nsParams);

    __global__ void
    computeRigidParticleVolume(DFSPHConstantParams *d_const,
                               DFSPHDynamicParams *d_data,
                               NeighborSearchUGConfig *d_nsConfig,
                               NeighborSearchUGParams *d_nsParams);

    __global__ void
    computeExtForce(DFSPHConstantParams *d_const,
                    DFSPHDynamicParams *d_data,
                    NeighborSearchUGConfig *d_nsConfig,
                    NeighborSearchUGParams *d_nsParams);

    __global__ void
    computeDensity(DFSPHConstantParams *d_const,
                   DFSPHDynamicParams *d_data,
                   NeighborSearchUGConfig *d_nsConfig,
                   NeighborSearchUGParams *d_nsParams);

    __global__ void
    predictDensity(DFSPHConstantParams *d_const,
                   DFSPHDynamicParams *d_data,
                   NeighborSearchUGConfig *d_nsConfig,
                   NeighborSearchUGParams *d_nsParams);

    __global__ void
    computeDivErr(DFSPHConstantParams *d_const,
                  DFSPHDynamicParams *d_data,
                  NeighborSearchUGConfig *d_nsConfig,
                  NeighborSearchUGParams *d_nsParams);

    __global__ void
    computeDFSPHAlpha(DFSPHConstantParams *d_const,
                      DFSPHDynamicParams *d_data,
                      NeighborSearchUGConfig *d_nsConfig,
                      NeighborSearchUGParams *d_nsParams);

    __global__ void
    adaptVelAdv_1(DFSPHConstantParams *d_const,
                  DFSPHDynamicParams *d_data,
                  NeighborSearchUGConfig *d_nsConfig,
                  NeighborSearchUGParams *d_nsParams);

    __global__ void
    advectPos(DFSPHConstantParams *d_const,
              DFSPHDynamicParams *d_data,
              NeighborSearchUGConfig *d_nsConfig,
              NeighborSearchUGParams *d_nsParams);

    __global__ void
    adaptVelAdv_2(DFSPHConstantParams *d_const,
                  DFSPHDynamicParams *d_data,
                  NeighborSearchUGConfig *d_nsConfig,
                  NeighborSearchUGParams *d_nsParams);

    __global__ void
    advectVel(DFSPHConstantParams *d_const,
              DFSPHDynamicParams *d_data,
              NeighborSearchUGConfig *d_nsConfig,
              NeighborSearchUGParams *d_nsParams);
}

#endif //SOSIM_DFSPH_CUDA_API_CUH
