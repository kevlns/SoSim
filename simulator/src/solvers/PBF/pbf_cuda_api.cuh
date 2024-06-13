//
// Created by ADMIN on 2024/6/13.
//

#ifndef SOSIM_PBF_CUDA_API_CUH
#define SOSIM_PBF_CUDA_API_CUH

#include "solvers/PBF/pbf_parameters.hpp"
#include "libs/NeighborSearchL/unified_grid_ns.hpp"

namespace SoSim {

    __host__ void
    init_data(PBFConstantParams &h_const,
              PBFConstantParams *d_const,
              PBFDynamicParams *d_data,
              NeighborSearchUGConfig *d_nsConfig,
              NeighborSearchUGParams *d_nsParams);

    __host__ void
    compute_sph_density_and_error(PBFConstantParams &h_const,
                                  PBFConstantParams *d_const,
                                  PBFDynamicParams *d_data,
                                  NeighborSearchUGConfig *d_nsConfig,
                                  NeighborSearchUGParams *d_nsParams);

    __host__ void
    update_lamb(PBFConstantParams &h_const,
                PBFConstantParams *d_const,
                PBFDynamicParams *d_data,
                NeighborSearchUGConfig *d_nsConfig,
                NeighborSearchUGParams *d_nsParams);

    __host__ void
    compute_dx(PBFConstantParams &h_const,
               PBFConstantParams *d_const,
               PBFDynamicParams *d_data,
               NeighborSearchUGConfig *d_nsConfig,
               NeighborSearchUGParams *d_nsParams);

    __host__ void
    apply_ext_force(PBFConstantParams &h_const,
                    PBFConstantParams *d_const,
                    PBFDynamicParams *d_data,
                    NeighborSearchUGConfig *d_nsConfig,
                    NeighborSearchUGParams *d_nsParams);

    __host__ void
    apply_dx(PBFConstantParams &h_const,
               PBFConstantParams *d_const,
               PBFDynamicParams *d_data,
               NeighborSearchUGConfig *d_nsConfig,
               NeighborSearchUGParams *d_nsParams);

    __host__ void
    post_correct(PBFConstantParams &h_const,
                 PBFConstantParams *d_const,
                 PBFDynamicParams *d_data,
                 NeighborSearchUGConfig *d_nsConfig,
                 NeighborSearchUGParams *d_nsParams);
}

#endif //SOSIM_PBF_CUDA_API_CUH
