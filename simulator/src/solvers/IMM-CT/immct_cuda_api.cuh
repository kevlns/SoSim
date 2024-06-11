//
// Created by ADMIN on 2024/3/26.
//

#ifndef SOSIM_IMMCT_CUDA_API_CUH
#define SOSIM_IMMCT_CUDA_API_CUH

#include "solvers/IMM-CT/immct_parameters.hpp"
#include "libs/NeighborSearchL/unified_grid_ns.hpp"

namespace SoSim {

    __host__ void
    init_data(IMMCTConstantParams &h_const,
              IMMCTConstantParams *d_const,
              IMMCTDynamicParams *d_data,
              NeighborSearchUGParams *d_nsParams);

    __host__ void
    prepare_ims(IMMCTConstantParams &h_const,
                IMMCTConstantParams *d_const,
                IMMCTDynamicParams *d_data,
                NeighborSearchUGConfig *d_nsConfig,
                NeighborSearchUGParams *d_nsParams);

    __host__ void
    sph_precompute(IMMCTConstantParams &h_const,
                   IMMCTConstantParams *d_const,
                   IMMCTDynamicParams *d_data,
                   NeighborSearchUGConfig *d_nsConfig,
                   NeighborSearchUGParams *d_nsParams);

    __host__ void
    vfsph_div(IMMCTConstantParams &h_const,
              IMMCTDynamicParams &h_data,
              Vec3ui &obj_part_index,
              IMMCTConstantParams *d_const,
              IMMCTDynamicParams *d_data,
              NeighborSearchUGConfig *d_nsConfig,
              NeighborSearchUGParams *d_nsParams,
              bool &crash);

    __host__ void
    apply_pressure_acc(IMMCTConstantParams &h_const,
                       IMMCTConstantParams *d_const,
                       IMMCTDynamicParams *d_data,
                       NeighborSearchUGParams *d_nsParams);

    __host__ void
    ism_gravity_vis_surface(IMMCTConstantParams &h_const,
                            IMMCTConstantParams *d_const,
                            IMMCTDynamicParams *d_data,
                            NeighborSearchUGConfig *d_nsConfig,
                            NeighborSearchUGParams *d_nsParams);

    __host__ void
    dfsph_gravity_vis_surface(IMMCTConstantParams &h_const,
                              IMMCTConstantParams *d_const,
                              IMMCTDynamicParams *d_data,
                              NeighborSearchUGConfig *d_nsConfig,
                              NeighborSearchUGParams *d_nsParams);

    __host__ void
    vfsph_incomp(IMMCTConstantParams &h_const,
                 IMMCTDynamicParams &h_data,
                 Vec3ui &obj_part_index,
                 IMMCTConstantParams *d_const,
                 IMMCTDynamicParams *d_data,
                 NeighborSearchUGConfig *d_nsConfig,
                 NeighborSearchUGParams *d_nsParams,
                 bool &crash);

    __host__ void
    update_pos(IMMCTConstantParams &h_const,
               IMMCTConstantParams *d_const,
               IMMCTDynamicParams *d_data,
               NeighborSearchUGParams *d_nsParams);

    __host__ void
    phase_transport_ism(IMMCTConstantParams &h_const,
                        IMMCTConstantParams *d_const,
                        IMMCTDynamicParams *d_data,
                        NeighborSearchUGConfig *d_nsConfig,
                        NeighborSearchUGParams *d_nsParams,
                        bool &crash);

    __host__ void
    update_mass_and_vel(IMMCTConstantParams &h_const,
                        IMMCTConstantParams *d_const,
                        IMMCTDynamicParams *d_data,
                        NeighborSearchUGParams *d_nsParams);

    __host__ void
    update_color(IMMCTConstantParams
                 &h_const,
                 IMMCTConstantParams *d_const,
                 IMMCTDynamicParams
                 *d_data,
                 NeighborSearchUGParams *d_nsParams);

    __host__ void
    ism_viscoelastic(IMMCTConstantParams &h_const,
                     IMMCTConstantParams *d_const,
                     IMMCTDynamicParams *d_data,
                     NeighborSearchUGConfig *d_nsConfig,
                     NeighborSearchUGParams *d_nsParams);

    __host__ void
    update_CT_parameters(IMMCTConstantParams &h_const,
                         IMMCTConstantParams *d_const,
                         IMMCTDynamicParams *d_data,
                         NeighborSearchUGConfig *d_nsConfig,
                         NeighborSearchUGParams *d_nsParams);

    __host__ void
    artificial_vis_bound(IMMCTConstantParams &h_const,
                         IMMCTConstantParams *d_const,
                         IMMCTDynamicParams *d_data,
                         NeighborSearchUGConfig *d_nsConfig,
                         NeighborSearchUGParams *d_nsParams);

}

namespace SoSim { // extra func
    __host__ void
    stirring(IMMCTConstantParams &h_const,
             IMMCTConstantParams *d_const,
             IMMCTDynamicParams *d_data,
             NeighborSearchUGParams *d_nsParams);

    __host__ void
    rotate_bowl(IMMCTConstantParams &h_const,
                IMMCTConstantParams *d_const,
                IMMCTDynamicParams *d_data,
                NeighborSearchUGParams *d_nsParams);

    __host__ void
    buckling(IMMCTConstantParams &h_const,
             IMMCTConstantParams *d_const,
             IMMCTDynamicParams *d_data,
             NeighborSearchUGParams *d_nsParams);
}

#endif //SOSIM_IMMCT_CUDA_API_CUH
