//
// Created by ADMIN on 2024/3/26.
//

#ifndef SOSIM_IMM_CUDA_API_CUH
#define SOSIM_IMM_CUDA_API_CUH

#include "solvers/IMM/imm_parameters.hpp"
#include "libs/NeighborSearchL/unified_grid_ns.hpp"

namespace SoSim {

    __host__ void
    init_data(IMMConstantParams &h_const,
              IMMConstantParams *d_const,
              IMMDynamicParams *d_data,
              NeighborSearchUGParams *d_nsParams);

    __host__ void
    prepare_ims(IMMConstantParams &h_const,
                IMMConstantParams *d_const,
                IMMDynamicParams *d_data,
                NeighborSearchUGConfig *d_nsConfig,
                NeighborSearchUGParams *d_nsParams);

    __host__ void
    sph_precompute(IMMConstantParams &h_const,
                   IMMConstantParams *d_const,
                   IMMDynamicParams *d_data,
                   NeighborSearchUGConfig *d_nsConfig,
                   NeighborSearchUGParams *d_nsParams);

    __host__ void
    vfsph_div(IMMConstantParams &h_const,
              IMMDynamicParams &h_data,
              Vec3ui &obj_part_index,
              IMMConstantParams *d_const,
              IMMDynamicParams *d_data,
              NeighborSearchUGConfig *d_nsConfig,
              NeighborSearchUGParams *d_nsParams,
              bool &crash);

    __host__ void
    apply_pressure_acc(IMMConstantParams &h_const,
                       IMMConstantParams *d_const,
                       IMMDynamicParams *d_data,
                       NeighborSearchUGParams *d_nsParams);

    __host__ void
    ism_gravity_vis_surface(IMMConstantParams &h_const,
                            IMMConstantParams *d_const,
                            IMMDynamicParams *d_data,
                            NeighborSearchUGConfig *d_nsConfig,
                            NeighborSearchUGParams *d_nsParams);

    __host__ void
    dfsph_gravity_vis_surface(IMMConstantParams &h_const,
                              IMMConstantParams *d_const,
                              IMMDynamicParams *d_data,
                              NeighborSearchUGConfig *d_nsConfig,
                              NeighborSearchUGParams *d_nsParams);

    __host__ void
    vfsph_incomp(IMMConstantParams &h_const,
                 IMMDynamicParams &h_data,
                 Vec3ui &obj_part_index,
                 IMMConstantParams *d_const,
                 IMMDynamicParams *d_data,
                 NeighborSearchUGConfig *d_nsConfig,
                 NeighborSearchUGParams *d_nsParams,
                 bool &crash);

    __host__ void
    update_pos(IMMConstantParams &h_const,
               IMMConstantParams *d_const,
               IMMDynamicParams *d_data,
               NeighborSearchUGParams *d_nsParams);

    __host__ void
    phase_transport_ism(IMMConstantParams &h_const,
                        IMMConstantParams *d_const,
                        IMMDynamicParams *d_data,
                        NeighborSearchUGConfig *d_nsConfig,
                        NeighborSearchUGParams *d_nsParams,
                        bool &crash);

    __host__ void
    update_mass_and_vel(IMMConstantParams &h_const,
                        IMMConstantParams *d_const,
                        IMMDynamicParams *d_data,
                        NeighborSearchUGParams *d_nsParams);

    __host__ void
    update_color(IMMConstantParams
                 &h_const,
                 IMMConstantParams *d_const,
                 IMMDynamicParams
                 *d_data,
                 NeighborSearchUGParams *d_nsParams);

    __host__ void
    ism_viscoelastic(IMMConstantParams &h_const,
                     IMMConstantParams *d_const,
                     IMMDynamicParams *d_data,
                     NeighborSearchUGConfig *d_nsConfig,
                     NeighborSearchUGParams *d_nsParams);

    __host__ void
    update_CT_parameters(IMMConstantParams &h_const,
                         IMMConstantParams *d_const,
                         IMMDynamicParams *d_data,
                         NeighborSearchUGConfig *d_nsConfig,
                         NeighborSearchUGParams *d_nsParams);

    __host__ void
    artificial_vis_bound(IMMConstantParams &h_const,
                         IMMConstantParams *d_const,
                         IMMDynamicParams *d_data,
                         NeighborSearchUGConfig *d_nsConfig,
                         NeighborSearchUGParams *d_nsParams);

}

namespace SoSim { // extra func
    __host__ void
    stirring(IMMConstantParams &h_const,
             IMMConstantParams *d_const,
             IMMDynamicParams *d_data,
             NeighborSearchUGParams *d_nsParams);

    __host__ void
    rotate_bowl(IMMConstantParams &h_const,
                IMMConstantParams *d_const,
                IMMDynamicParams *d_data,
                NeighborSearchUGParams *d_nsParams);

    __host__ void
    buckling(IMMConstantParams &h_const,
             IMMConstantParams *d_const,
             IMMDynamicParams *d_data,
             NeighborSearchUGParams *d_nsParams);
}

#endif //SOSIM_IMM_CUDA_API_CUH
