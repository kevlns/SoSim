//
// Created by ADMIN on 2024/3/13.
//

#ifndef SOSIM_DFSPH_CUDA_INVOKE_API_CUH
#define SOSIM_DFSPH_CUDA_INVOKE_API_CUH

#include "solvers/DFSPH/dfsph_solver.hpp"

namespace SoSim {

    void Init(DFSPHSolverConfig *config,
              DFSPHConstantParams *d_const,
              DFSPHDynamicParams *d_data,
              NeighborSearchUGConfig *d_nsConfig,
              NeighborSearchUGParams *d_nsParams);

    void ComputeRigidParticleVolume(DFSPHSolverConfig *config,
                                    DFSPHConstantParams *d_const,
                                    DFSPHDynamicParams *d_data,
                                    NeighborSearchUGConfig *d_nsConfig,
                                    NeighborSearchUGParams *d_nsParams);

    void ComputeExtForce(DFSPHSolverConfig *config,
                         DFSPHConstantParams *d_const,
                         DFSPHDynamicParams *d_data,
                         NeighborSearchUGConfig *d_nsConfig,
                         NeighborSearchUGParams *d_nsParams);

    void ComputeDensityAndDFSPHAlpha(DFSPHSolverConfig *config,
                                     DFSPHConstantParams *d_const,
                                     DFSPHDynamicParams *d_data,
                                     NeighborSearchUGConfig *d_nsConfig,
                                     NeighborSearchUGParams *d_nsParams);

    void CorrectDivErr(DFSPHSolverConfig *config,
                       const DFSPHDynamicParams &h_data,
                       DFSPHConstantParams *d_const,
                       DFSPHDynamicParams *d_data,
                       NeighborSearchUGConfig *d_nsConfig,
                       NeighborSearchUGParams *d_nsParams);

    void CorrectDensityErr(DFSPHSolverConfig *config,
                           const DFSPHDynamicParams &h_data,
                           DFSPHConstantParams *d_const,
                           DFSPHDynamicParams *d_data,
                           NeighborSearchUGConfig *d_nsConfig,
                           NeighborSearchUGParams *d_nsParams);

    void AdvectPos(DFSPHSolverConfig *config,
                   DFSPHConstantParams *d_const,
                   DFSPHDynamicParams *d_data,
                   NeighborSearchUGConfig *d_nsConfig,
                   NeighborSearchUGParams *d_nsParams);

    void AdvectVel(DFSPHSolverConfig *config,
                   DFSPHConstantParams *d_const,
                   DFSPHDynamicParams *d_data,
                   NeighborSearchUGConfig *d_nsConfig,
                   NeighborSearchUGParams *d_nsParams);
}

#endif //SOSIM_DFSPH_CUDA_INVOKE_API_CUH
