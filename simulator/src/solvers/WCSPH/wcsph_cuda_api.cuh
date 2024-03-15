//
// Created by ADMIN on 2024/3/14.
//

#ifndef SOSIM_WCSPH_CUDA_API_CUH
#define SOSIM_WCSPH_CUDA_API_CUH

#include "solvers/WCSPH/wcsph_solver.hpp"
#include "libs/NeighborSearchL/unified_grid_ns.hpp"

namespace SoSim {
    /*
     *  Cuda API
     */
    __device__ inline float cubic_value(const Vec3f &r, float h);

    __device__ inline Vec3f cubic_gradient(const Vec3f &r, float h);

    __global__ void
    init_cuda(WCSPHConstantParams *d_const,
              WCSPHDynamicParams *d_data,
              NeighborSearchUGConfig *d_nsConfig,
              NeighborSearchUGParams *d_nsParams);

    __global__ void
    computeDensityAndPressure_cuda(WCSPHConstantParams *d_const,
                                   WCSPHDynamicParams *d_data,
                                   NeighborSearchUGConfig *d_nsConfig,
                                   NeighborSearchUGParams *d_nsParams);


    __global__ void
    computeGravityForce_cuda(WCSPHConstantParams *d_const,
                             WCSPHDynamicParams *d_data,
                             NeighborSearchUGConfig *d_nsConfig,
                             NeighborSearchUGParams *d_nsParams);

    __global__ void
    computePressureForce_cuda(WCSPHConstantParams *d_const,
                              WCSPHDynamicParams *d_data,
                              NeighborSearchUGConfig *d_nsConfig,
                              NeighborSearchUGParams *d_nsParams);

    __global__ void
    computeViscousForce_cuda(WCSPHConstantParams *d_const,
                             WCSPHDynamicParams *d_data,
                             NeighborSearchUGConfig *d_nsConfig,
                             NeighborSearchUGParams *d_nsParams);

    __global__ void
    advect_cuda(WCSPHConstantParams *d_const,
                WCSPHDynamicParams *d_data,
                NeighborSearchUGConfig *d_nsConfig,
                NeighborSearchUGParams *d_nsParams);

    __global__ void
    wcsphOverall_cuda();

}

namespace SoSim {
    /* ============================================================================================================
     * Host Invoke API
     */
    __host__ void
    computeDensityAndPressure();

    __host__ void
    computeGravityForce();

    __host__ void
    computePressureForce();

    __host__ void
    computeViscousForce();

    __host__ void
    advect();

    __host__ void
    wcsphOverall();

}

#endif //SOSIM_WCSPH_CUDA_API_CUH
