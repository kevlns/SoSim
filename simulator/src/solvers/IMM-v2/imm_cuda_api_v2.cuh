//
// Created by ADMIN on 2024/3/26.
//

#ifndef SOSIM_IMM_CUDA_API_v2_CUH
#define SOSIM_IMM_CUDA_API_v2_CUH

#include "solvers/IMM-v2/imm_solver_v2.hpp"
#include "solvers/IMM-v2/imm_parameters_v2.hpp"
#include "libs/NeighborSearchL/unified_grid_ns.hpp"

namespace SoSim {

    __host__ void
    init_data(IMMConstantParams_v2 &h_const,
              IMMConstantParams_v2 *d_const,
              IMMDynamicParams_v2 *d_data,
              float *d_phase_density,
              Vec3f *d_phase_color,
              float *d_phase_vis);

    __host__ void
    prepare_imm(IMMConstantParams_v2 &h_const,
                IMMConstantParams_v2 *d_const,
                IMMDynamicParams_v2 *d_data,
                NeighborSearchUGConfig *d_nsConfig,
                NeighborSearchUGParams *d_nsParams);

    __host__ void
    sph_precompute(IMMConstantParams_v2 &h_const,
                   IMMConstantParams_v2 *d_const,
                   IMMDynamicParams_v2 *d_data,
                   NeighborSearchUGConfig *d_nsConfig,
                   NeighborSearchUGParams *d_nsParams);

    __host__ void
    vfsph_div(IMMConstantParams_v2 &h_const,
              IMMDynamicParams_v2 &h_data,
              Vec3ui &obj_part_index,
              IMMConstantParams_v2 *d_const,
              IMMDynamicParams_v2 *d_data,
              NeighborSearchUGConfig *d_nsConfig,
              NeighborSearchUGParams *d_nsParams,
              bool &crash);

    __host__ void
    apply_pressure_acc(IMMConstantParams_v2 &h_const,
                       IMMConstantParams_v2 *d_const,
                       IMMDynamicParams_v2 *d_data,
                       NeighborSearchUGParams *d_nsParams);

    __host__ void
    imm_gravity_vis_surface(IMMConstantParams_v2 &h_const,
                            IMMConstantParams_v2 *d_const,
                            IMMDynamicParams_v2 *d_data,
                            NeighborSearchUGConfig *d_nsConfig,
                            NeighborSearchUGParams *d_nsParams);

    __host__ void
    vfsph_incomp(IMMConstantParams_v2 &h_const,
                 IMMDynamicParams_v2 &h_data,
                 Vec3ui &obj_part_index,
                 IMMConstantParams_v2 *d_const,
                 IMMDynamicParams_v2 *d_data,
                 NeighborSearchUGConfig *d_nsConfig,
                 NeighborSearchUGParams *d_nsParams,
                 bool &crash);

    __host__ void
    update_pos(IMMConstantParams_v2 &h_const,
               IMMConstantParams_v2 *d_const,
               IMMDynamicParams_v2 *d_data,
               NeighborSearchUGParams *d_nsParams);

    __host__ void
    phase_transport_ism(IMMConstantParams_v2 &h_const,
                        IMMConstantParams_v2 *d_const,
                        IMMDynamicParams_v2 *d_data,
                        NeighborSearchUGConfig *d_nsConfig,
                        NeighborSearchUGParams *d_nsParams,
                        bool &crash);

    __host__ void
    update_mass_and_vel(IMMConstantParams_v2 &h_const,
                        IMMConstantParams_v2 *d_const,
                        IMMDynamicParams_v2 *d_data,
                        NeighborSearchUGParams *d_nsParams);

    __host__ void
    update_color(IMMConstantParams_v2 &h_const,
                 IMMConstantParams_v2 *d_const,
                 IMMDynamicParams_v2 *d_data,
                 NeighborSearchUGParams *d_nsParams);
}

#endif //SOSIM_IMM_CUDA_API_v2_CUH
