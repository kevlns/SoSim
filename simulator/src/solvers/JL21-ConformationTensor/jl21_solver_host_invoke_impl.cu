//@author        : Long Shen
//@date          : 2024/1/14
//@description   :
//@version       : 1.0

#include "jl21_solver_host_invoke_api.cuh"

namespace SoSim {

    extern void
    CheckParticleActivity(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                          NeighborSearchUGConfig *d_nsConst) {
        checkParticleActivity<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data, d_nsConst);
    }

    extern void
    InitOtherData(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data) {
        initOtherData<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data);
    }

    extern void
    ComputeBoundaryVolume(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                          NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams) {
        computeBoundaryVolume<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data, d_nsConst, d_nsParams);
    }

    extern void
    EstimateDensity(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                    NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams) {
        estimateDensity<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data, d_nsConst, d_nsParams);
    }

    extern void
    ComputeFluidVisForce(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                         NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams) {
        computeFluidVisForce<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data, d_nsConst, d_nsParams);
    }

    extern void
    ComputeFluidPressureForce(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                              NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams) {
        computeFluidPressureForce<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data, d_nsConst, d_nsParams);
    }

    extern void
    ComputePhaseVel(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                    NeighborSearchUGParams *d_nsParams) {
        computePhaseVel<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data, d_nsParams);
    }

    extern void
    ComputeMixVel(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                  NeighborSearchUGParams *d_nsParams) {
        computeMixVel<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data, d_nsParams);
    }

    extern void
    ComputeDragForce(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                     NeighborSearchUGParams *d_nsParams) {
        computeDragForce<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data, d_nsParams);
    }

    extern void
    Advect(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
           NeighborSearchUGParams *d_nsParams) {
        advect<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data, d_nsParams);
    }

    extern void
    Reset(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams) {
        reset<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data, d_nsParams);
    }

    extern void
    PhaseTransport(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                   NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams) {
        for (int i = 0; i < 2; ++i) {
            computePhaseTransport<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data, d_nsConst, d_nsParams);

            updateLambda<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data, d_nsParams);
        }

        phaseTransportFinalize<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data, d_nsParams);
    }

    extern void
    UpdateColor(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                NeighborSearchUGParams *d_nsParams) {
        updateColor<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data, d_nsParams);
    }

    extern void
    ComputeStressForce(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                       NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams) {
        computeVelGrad<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data, d_nsConst, d_nsParams);

        computeStress<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data, d_nsParams);

        computeStressForce<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data, d_nsConst, d_nsParams);
    }

    extern void
    ComputeVelGrad(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                   NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams) {
        computeVelGrad<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data, d_nsConst, d_nsParams);
    }

    extern void
    ComputeCompliance(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                      NeighborSearchUGParams *d_nsParams) {
        computeCompliance<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data, d_nsParams);
    }

    extern void
    ComputeConformationTensorForce(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                                   NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams) {
        updateConformationTensor<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data, d_nsConst, d_nsParams);

        computeElasticStress<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data, d_nsParams);

        computeStressForce<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data, d_nsConst, d_nsParams);
    }

    extern void
    ComputeCohesionForce(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
                         NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams) {
        computeCohesionForce<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data, d_nsConst, d_nsParams);
    }

    extern void
    CylinderRotate(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data) {
        cylinderRotate<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data);
    }

    extern void
    MoveStirringRod(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data) {
        moveStirringRod<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data);
        //    moveStirringRod_test<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data);
    }

    extern void
    MoveBowl(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data) {
        moveBowl<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data);
    }

    extern void
    MatTransfer(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data) {
        matTransfer<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data);
    }

    extern void
    MoveRigid(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data) {
        moveRigid<<<h_const.block_num, h_const.thread_num>>>(d_const, d_data);
    }

    extern void
    UpdateTime(JL21CTConstantParams *d_const) {
        updateTime<<<1, 1>>>(d_const);
    }

}