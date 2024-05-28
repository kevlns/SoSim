//
// Created by ADMIN on 2024/3/26.
//

#ifndef SOSIM_ISM_CT_CUDA_API_CUH
#define SOSIM_ISM_CT_CUDA_API_CUH

#include "solvers/ISM-ConformationTensor/ism_ct_solver.hpp"
#include "libs/NeighborSearchL/unified_grid_ns.hpp"

namespace SoSim {

    __host__ void
    init_data(IMSCTConstantParams &h_const,
              IMSCTConstantParams *d_const,
              IMSCTDynamicParams *d_data,
              NeighborSearchUGParams *d_nsParams);

    __host__ void
    prepare_ims(IMSCTConstantParams &h_const,
                IMSCTConstantParams *d_const,
                IMSCTDynamicParams *d_data,
                NeighborSearchUGConfig *d_nsConfig,
                NeighborSearchUGParams *d_nsParams);

    __host__ void
    sph_precompute(IMSCTConstantParams &h_const,
                   IMSCTConstantParams *d_const,
                   IMSCTDynamicParams *d_data,
                   NeighborSearchUGConfig *d_nsConfig,
                   NeighborSearchUGParams *d_nsParams);

    __host__ void
    vfsph_div(IMSCTConstantParams &h_const,
              IMSCTDynamicParams &h_data,
              Vec3ui &obj_part_index,
              IMSCTConstantParams *d_const,
              IMSCTDynamicParams *d_data,
              NeighborSearchUGConfig *d_nsConfig,
              NeighborSearchUGParams *d_nsParams,
              bool &crash);

    __host__ void
    apply_pressure_acc(IMSCTConstantParams &h_const,
                       IMSCTConstantParams *d_const,
                       IMSCTDynamicParams *d_data,
                       NeighborSearchUGParams *d_nsParams);

    __host__ void
    ism_gravity_vis_surface(IMSCTConstantParams &h_const,
                            IMSCTConstantParams *d_const,
                            IMSCTDynamicParams *d_data,
                            NeighborSearchUGConfig *d_nsConfig,
                            NeighborSearchUGParams *d_nsParams);

    __host__ void
    dfsph_gravity_vis_surface(IMSCTConstantParams &h_const,
                              IMSCTConstantParams *d_const,
                              IMSCTDynamicParams *d_data,
                              NeighborSearchUGConfig *d_nsConfig,
                              NeighborSearchUGParams *d_nsParams);

    __host__ void
    vfsph_incomp(IMSCTConstantParams &h_const,
                 IMSCTDynamicParams &h_data,
                 Vec3ui &obj_part_index,
                 IMSCTConstantParams *d_const,
                 IMSCTDynamicParams *d_data,
                 NeighborSearchUGConfig *d_nsConfig,
                 NeighborSearchUGParams *d_nsParams,
                 bool &crash);

    __host__ void
    update_pos(IMSCTConstantParams &h_const,
               IMSCTConstantParams *d_const,
               IMSCTDynamicParams *d_data,
               NeighborSearchUGParams *d_nsParams);

    __host__ void
    phase_transport_ism(IMSCTConstantParams &h_const,
                        IMSCTConstantParams *d_const,
                        IMSCTDynamicParams *d_data,
                        NeighborSearchUGConfig *d_nsConfig,
                        NeighborSearchUGParams *d_nsParams,
                        bool &crash);

    __host__ void
    update_mass_and_vel(IMSCTConstantParams &h_const,
                        IMSCTConstantParams *d_const,
                        IMSCTDynamicParams *d_data,
                        NeighborSearchUGParams *d_nsParams);

    __host__ void
    update_color(IMSCTConstantParams
                 &h_const,
                 IMSCTConstantParams *d_const,
                 IMSCTDynamicParams
                 *d_data,
                 NeighborSearchUGParams *d_nsParams);

    __host__ void
    ism_viscoelastic(IMSCTConstantParams &h_const,
                     IMSCTConstantParams *d_const,
                     IMSCTDynamicParams *d_data,
                     NeighborSearchUGConfig *d_nsConfig,
                     NeighborSearchUGParams *d_nsParams);

    __host__ void
    update_CT_parameters(IMSCTConstantParams &h_const,
                         IMSCTConstantParams *d_const,
                         IMSCTDynamicParams *d_data,
                         NeighborSearchUGConfig *d_nsConfig,
                         NeighborSearchUGParams *d_nsParams);

    __host__ void
    artificial_vis_bound(IMSCTConstantParams &h_const,
                         IMSCTConstantParams *d_const,
                         IMSCTDynamicParams *d_data,
                         NeighborSearchUGConfig *d_nsConfig,
                         NeighborSearchUGParams *d_nsParams);

}

namespace SoSim { // extra func
    __host__ void
    stirring(IMSCTConstantParams &h_const,
             IMSCTConstantParams *d_const,
             IMSCTDynamicParams *d_data,
             NeighborSearchUGParams *d_nsParams);

    __host__ void
    rotate_bowl(IMSCTConstantParams &h_const,
                IMSCTConstantParams *d_const,
                IMSCTDynamicParams *d_data,
                NeighborSearchUGParams *d_nsParams);

    __host__ void
    buckling(IMSCTConstantParams &h_const,
             IMSCTConstantParams *d_const,
             IMSCTDynamicParams *d_data,
             NeighborSearchUGParams *d_nsParams);
}

#endif //SOSIM_ISM_CT_CUDA_API_CUH
