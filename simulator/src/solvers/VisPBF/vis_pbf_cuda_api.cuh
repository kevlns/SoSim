//
// Created by ADMIN on 2024/6/13.
//

#ifndef SOSIM_VIS_PBF_CUDA_API_CUH
#define SOSIM_VIS_PBF_CUDA_API_CUH

#include "solvers/VisPBF/vis_pbf_parameters.hpp"
#include "libs/NeighborSearchL/unified_grid_ns.hpp"

namespace SoSim {

    __host__ void
    init_data(VisPBFConstantParams &h_const,
              VisPBFConstantParams *d_const,
              VisPBFDynamicParams *d_data,
              NeighborSearchUGConfig *d_nsConfig,
              NeighborSearchUGParams *d_nsParams);

    __host__ void
    compute_sph_density_and_error(VisPBFConstantParams &h_const,
                                  VisPBFConstantParams *d_const,
                                  VisPBFDynamicParams *d_data,
                                  NeighborSearchUGConfig *d_nsConfig,
                                  NeighborSearchUGParams *d_nsParams);

    __host__ void
    update_lamb(VisPBFConstantParams &h_const,
                VisPBFConstantParams *d_const,
                VisPBFDynamicParams *d_data,
                NeighborSearchUGConfig *d_nsConfig,
                NeighborSearchUGParams *d_nsParams);

    __host__ void
    compute_dx(VisPBFConstantParams &h_const,
               VisPBFConstantParams *d_const,
               VisPBFDynamicParams *d_data,
               NeighborSearchUGConfig *d_nsConfig,
               NeighborSearchUGParams *d_nsParams);

    __host__ void
    apply_ext_force(VisPBFConstantParams &h_const,
                    VisPBFConstantParams *d_const,
                    VisPBFDynamicParams *d_data,
                    NeighborSearchUGConfig *d_nsConfig,
                    NeighborSearchUGParams *d_nsParams);

    __host__ void
    apply_dx(VisPBFConstantParams &h_const,
               VisPBFConstantParams *d_const,
               VisPBFDynamicParams *d_data,
               NeighborSearchUGConfig *d_nsConfig,
               NeighborSearchUGParams *d_nsParams);

    __host__ void
    post_correct(VisPBFConstantParams &h_const,
                 VisPBFConstantParams *d_const,
                 VisPBFDynamicParams *d_data,
                 NeighborSearchUGConfig *d_nsConfig,
                 NeighborSearchUGParams *d_nsParams);
}

#endif //SOSIM_VIS_PBF_CUDA_API_CUH
