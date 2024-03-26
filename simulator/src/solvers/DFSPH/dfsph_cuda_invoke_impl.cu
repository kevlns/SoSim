//
// Created by ADMIN on 2024/3/13.
//

#include "dfsph_cuda_invoke_api.cuh"
#include "dfsph_cuda_api.cuh"
#include "libs/AnalysisL/statistic_util.hpp"

namespace SoSim {

    void Init(DFSPHSolverConfig *config,
              DFSPHConstantParams *d_const,
              DFSPHDynamicParams *d_data,
              NeighborSearchUGConfig *d_nsConfig,
              NeighborSearchUGParams *d_nsParams) {
        init<<<config->kernel_blocks, config->kernel_threads>>>(d_const,
                                                                d_data,
                                                                d_nsConfig,
                                                                d_nsParams);
        cudaGetLastError();
    }

    void ComputeRigidParticleVolume(DFSPHSolverConfig *config,
                                    DFSPHConstantParams *d_const,
                                    DFSPHDynamicParams *d_data,
                                    NeighborSearchUGConfig *d_nsConfig,
                                    NeighborSearchUGParams *d_nsParams) {
        computeRigidParticleVolume<<<config->kernel_blocks, config->kernel_threads>>>(d_const,
                                                                                      d_data,
                                                                                      d_nsConfig,
                                                                                      d_nsParams);
        cudaGetLastError();
    }

    void ComputeExtForce(DFSPHSolverConfig *config,
                         DFSPHConstantParams *d_const,
                         DFSPHDynamicParams *d_data,
                         NeighborSearchUGConfig *d_nsConfig,
                         NeighborSearchUGParams *d_nsParams) {
        computeExtForce<<<config->kernel_blocks, config->kernel_threads>>>(d_const,
                                                                           d_data,
                                                                           d_nsConfig,
                                                                           d_nsParams);
        cudaGetLastError();
    }

    void ComputeDensityAndDFSPHAlpha(DFSPHSolverConfig *config,
                                     DFSPHConstantParams *d_const,
                                     DFSPHDynamicParams *d_data,
                                     NeighborSearchUGConfig *d_nsConfig,
                                     NeighborSearchUGParams *d_nsParams) {
        computeDensity<<<config->kernel_blocks, config->kernel_threads>>>(d_const,
                                                                          d_data,
                                                                          d_nsConfig,
                                                                          d_nsParams);
        computeDFSPHAlpha<<<config->kernel_blocks, config->kernel_threads>>>(d_const,
                                                                             d_data,
                                                                             d_nsConfig,
                                                                             d_nsParams);
        cudaGetLastError();
    }

    void CorrectDivErr(DFSPHSolverConfig *config,
                       const DFSPHDynamicParams &h_data,
                       DFSPHConstantParams *d_const,
                       DFSPHDynamicParams *d_data,
                       NeighborSearchUGConfig *d_nsConfig,
                       NeighborSearchUGParams *d_nsParams) {

        float div_err_0 = 1e-5;
        float div_err_avg;
        unsigned iter = 0;
        std::vector<float> div_err(config->particle_num_copy);

        do {
            div_err_avg = 0;
            computeDivErr<<<config->kernel_blocks, config->kernel_threads>>>(d_const,
                                                                             d_data,
                                                                             d_nsConfig,
                                                                             d_nsParams);

            adaptVelAdv_1<<<config->kernel_blocks, config->kernel_threads>>>(d_const,
                                                                             d_data,
                                                                             d_nsConfig,
                                                                             d_nsParams);

            cudaMemcpy(div_err.data(), h_data.div_err, config->particle_num_copy * sizeof(float),
                       cudaMemcpyDeviceToHost);
            for (auto &err: div_err)
                div_err_avg += err;
            div_err_avg /= config->particle_num_copy;

            iter++;
        } while (div_err_avg > div_err_0);

        cudaGetLastError();
    }

    void CorrectDensityErr(DFSPHSolverConfig *config,
                           const DFSPHDynamicParams &h_data,
                           DFSPHConstantParams *d_const,
                           DFSPHDynamicParams *d_data,
                           NeighborSearchUGConfig *d_nsConfig,
                           NeighborSearchUGParams *d_nsParams) {

        float density_err_0 = 1e-3;
        float density_avg;
        unsigned iter = 0;
        std::vector<float> density(config->particle_num_copy);

        do {
            predictDensity<<<config->kernel_blocks, config->kernel_threads>>>(d_const,
                                                                              d_data,
                                                                              d_nsConfig,
                                                                              d_nsParams);

            adaptVelAdv_2<<<config->kernel_blocks, config->kernel_threads>>>(d_const,
                                                                             d_data,
                                                                             d_nsConfig,
                                                                             d_nsParams);

            density_avg = 0;
            std::vector<float> density_sph(config->particle_num_copy);
            cudaMemcpy(density_sph.data(), h_data.density_sph, config->particle_num_copy * sizeof(float),
                       cudaMemcpyDeviceToHost);
            for (auto &dens: density_sph)
                density_avg += dens;
            density_avg /= config->particle_num_copy;

            iter++;
        } while ((density_avg - config->rest_density) > density_err_0);

        cudaGetLastError();
    }

    void AdvectPos(DFSPHSolverConfig *config,
                   DFSPHConstantParams *d_const,
                   DFSPHDynamicParams *d_data,
                   NeighborSearchUGConfig *d_nsConfig,
                   NeighborSearchUGParams *d_nsParams) {
        advectPos<<<config->kernel_blocks, config->kernel_threads>>>(d_const,
                                                                     d_data,
                                                                     d_nsConfig,
                                                                     d_nsParams);
        cudaGetLastError();
    }

    void AdvectVel(DFSPHSolverConfig *config,
                   DFSPHConstantParams *d_const,
                   DFSPHDynamicParams *d_data,
                   NeighborSearchUGConfig *d_nsConfig,
                   NeighborSearchUGParams *d_nsParams) {
        advectVel<<<config->kernel_blocks, config->kernel_threads>>>(d_const,
                                                                     d_data,
                                                                     d_nsConfig,
                                                                     d_nsParams);
        cudaGetLastError();
    }

}
