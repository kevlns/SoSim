//
// Created by ADMIN on 2024/3/26.
//

#ifndef SOSIM_DFSPH_CUDA_API_CUH
#define SOSIM_DFSPH_CUDA_API_CUH

#include "solvers/DFSPH/dfsph_solver.hpp"
#include "solvers/DFSPH/dfsph_parameters.hpp"
#include "libs/NeighborSearchL/unified_grid_ns.hpp"

namespace SoSim {

    __host__ void
    init_data(DFSPHConstantParams &h_const,
              DFSPHConstantParams *d_const,
              DFSPHDynamicParams *d_data,
              NeighborSearchUGParams *d_nsParams);

    __host__ void
    prepare_dfsph(DFSPHConstantParams &h_const,
                DFSPHConstantParams *d_const,
                DFSPHDynamicParams *d_data,
                NeighborSearchUGConfig *d_nsConfig,
                NeighborSearchUGParams *d_nsParams);

    __host__ void
    sph_precompute(DFSPHConstantParams &h_const,
                   DFSPHConstantParams *d_const,
                   DFSPHDynamicParams *d_data,
                   NeighborSearchUGConfig *d_nsConfig,
                   NeighborSearchUGParams *d_nsParams);

    __host__ void
    vfsph_div(DFSPHConstantParams &h_const,
              DFSPHDynamicParams &h_data,
              Vec3ui &obj_part_index,
              DFSPHConstantParams *d_const,
              DFSPHDynamicParams *d_data,
              NeighborSearchUGConfig *d_nsConfig,
              NeighborSearchUGParams *d_nsParams,
              bool &crash);

    __host__ void
    apply_pressure_acc(DFSPHConstantParams &h_const,
                       DFSPHConstantParams *d_const,
                       DFSPHDynamicParams *d_data,
                       NeighborSearchUGParams *d_nsParams);

    __host__ void
    dfsph_gravity_vis_surface(DFSPHConstantParams &h_const,
                              DFSPHConstantParams *d_const,
                              DFSPHDynamicParams *d_data,
                              NeighborSearchUGConfig *d_nsConfig,
                              NeighborSearchUGParams *d_nsParams);

    __host__ void
    vfsph_incomp(DFSPHConstantParams &h_const,
                 DFSPHDynamicParams &h_data,
                 Vec3ui &obj_part_index,
                 DFSPHConstantParams *d_const,
                 DFSPHDynamicParams *d_data,
                 NeighborSearchUGConfig *d_nsConfig,
                 NeighborSearchUGParams *d_nsParams,
                 bool &crash);

    __host__ void
    update_pos(DFSPHConstantParams &h_const,
               DFSPHConstantParams *d_const,
               DFSPHDynamicParams *d_data,
               NeighborSearchUGParams *d_nsParams);

    __host__ void
    artificial_vis_bound(DFSPHConstantParams &h_const,
                         DFSPHConstantParams *d_const,
                         DFSPHDynamicParams *d_data,
                         NeighborSearchUGConfig *d_nsConfig,
                         NeighborSearchUGParams *d_nsParams);

}

#endif //SOSIM_DFSPH_CUDA_API_CUH
