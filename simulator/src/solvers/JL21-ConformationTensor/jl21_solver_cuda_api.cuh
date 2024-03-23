//@author        : Long Shen
//@date          : 2024/1/14
//@description   :
//@version       : 1.0

#ifndef JL21_V1_JL21_SOLVER_CUDA_API_CUH
#define JL21_V1_JL21_SOLVER_CUDA_API_CUH

#include "solvers/JL21-ConformationTensor/jl21_ct_parameters.hpp"
#include "libs/NeighborSearchL/unified_grid_ns.hpp"

namespace SoSim {

    extern __global__ void
    checkParticleActivity(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                          NeighborSearchUGConfig *d_nsConst);

    extern __global__ void
    initOtherData(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data);

    extern __global__ void
    computeBoundaryVolume(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
                          NeighborSearchUGParams *d_nsParams);

    extern __global__ void
    estimateDensity(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
                    NeighborSearchUGParams *d_nsParams);

    extern __global__ void
    computeFluidVisForce(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
                         NeighborSearchUGParams *d_nsParams);

    extern __global__ void
    computeFluidPressureForce(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                              NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams);

    extern __global__ void
    computePhaseVel(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams);

    extern __global__ void
    computeMixVel(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams);

    extern __global__ void
    computeDragForce(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams);

    extern __global__ void
    advect(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams);

    extern __global__ void
    reset(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams);

    extern __global__ void
    computePhaseTransport(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
                          NeighborSearchUGParams *d_nsParams);

    extern __global__ void
    updateLambda(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams);

    extern __global__ void
    phaseTransportFinalize(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                           NeighborSearchUGParams *d_nsParams);

    extern __global__ void
    updateColor(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams);

    extern __global__ void
    computeVelGrad(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
                   NeighborSearchUGParams *d_nsParams);

    extern __global__ void
    computeStress(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams);

    extern __global__ void
    computeFluidStress(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
                       NeighborSearchUGParams *d_nsParams);

    extern __global__ void
    computeElastoplasticStress(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                               NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams);

    extern __global__ void
    computeStressForce(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
                       NeighborSearchUGParams *d_nsParams);

    extern __global__ void
    updateConformationTensor(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                             NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams);

    extern __global__ void
    computeElasticStress(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams);

    extern __global__ void
    computeCompliance(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams);

    extern __global__ void
    computeCohesionForce(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
                         NeighborSearchUGParams *d_nsParams);

    extern __global__ void
    matTransfer(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data);

    extern __global__ void
    updateTime(JL21CTConstantParams *d_const);

    extern __global__ void
    cylinderRotate(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data);

    extern __global__ void
    moveStirringRod(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data);

    extern __global__ void
    moveStirringRod_test(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data);

    extern __global__ void
    moveBowl(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data);

    extern __global__ void
    moveRigid(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data);

}
#endif//JL21_V1_JL21_SOLVER_CUDA_API_CUH
