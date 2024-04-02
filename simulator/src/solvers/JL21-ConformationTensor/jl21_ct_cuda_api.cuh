//
// Created by ADMIN on 2024/3/22.
//

#ifndef SOSIM_JL21_CT_CUDA_API_CUH
#define SOSIM_JL21_CT_CUDA_API_CUH

#include "solvers/JL21-ConformationTensor/jl21_ct_parameters.hpp"
#include "libs/NeighborSearchL/unified_grid_ns.hpp"

namespace SoSim {

    __host__ void
    InitOtherData(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data);

    __host__ void
    EstimateDensity(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                    NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams);

    __host__ void
    ComputeVelGrad(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                   NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams);

    __host__ void
    ComputeCompliance(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                      NeighborSearchUGParams *d_nsParams);

    __host__ void
    ComputeConformationTensorForce(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const,
                                   JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
                                   NeighborSearchUGParams *d_nsParams);

    __host__ void
    ComputeFluidVisForce(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                         NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams);

    __host__ void
    ComputeCohesionForce(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                         NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams);

    __host__ void
    ComputeFluidPressureForce(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                              NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams);

    __host__ void
    ComputePhaseVel(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                    NeighborSearchUGParams *d_nsParams);

    __host__ void
    ComputeMixVel(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                  NeighborSearchUGParams *d_nsParams);

    __host__ void
    ComputeDragForce(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                     NeighborSearchUGParams *d_nsParams);

    __host__ void
    Advect(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
           NeighborSearchUGParams *d_nsParams);

    __host__ void
    Reset(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
          NeighborSearchUGParams *d_nsParams);

    __host__ void
    PhaseTransport(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                   NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams);

    __host__ void
    UpdateColor(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                NeighborSearchUGParams *d_nsParams);

}

#endif //SOSIM_JL21_CT_CUDA_API_CUH
