//@author        : Long Shen
//@date          : 2024/1/14
//@description   :
//@version       : 1.0

#ifndef JL21_V1_JL21_SOLVER_HOST_INVOKE_API_CUH
#define JL21_V1_JL21_SOLVER_HOST_INVOKE_API_CUH

#include "jl21_solver_cuda_api.cuh"

namespace SoSim{

extern void
CheckParticleActivity(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst);

extern void
InitOtherData(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data);

extern void
ComputeBoundaryVolume(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams);

extern void
EstimateDensity(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams);

extern void
ComputeFluidVisForce(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams);

extern void
ComputeFluidPressureForce(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams);

extern void
ComputePhaseVel(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams);

extern void
ComputeMixVel(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams);

extern void
ComputeDragForce(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams);

extern void
Advect(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams);

extern void
Reset(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams);

extern void
PhaseTransport(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams);

extern void
UpdateColor(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams);

extern void
ComputeStressForce(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams);

extern void
ComputeVelGrad(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams);

extern void
ComputeCompliance(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams);

extern void
ComputeConformationTensorForce(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams);

extern void
ComputeCohesionForce(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams);

extern void
CylinderRotate(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data);

extern void
MoveStirringRod(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data);

extern void
MoveBowl(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data);

extern void
MatTransfer(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data);

extern void
MoveRigid(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data);

extern void
UpdateTime(JL21CTConstantParams *d_const);

}

#endif//JL21_V1_JL21_SOLVER_HOST_INVOKE_API_CUH
