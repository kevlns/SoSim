//@author        : Long Shen
//@date          : 2023/11/12
//@description   :
//@version       : 1.0

#ifndef SOSIM_CUDA_API_CUH
#define SOSIM_CUDA_API_CUH

#include "Private/PhysicalSolvers/MixtureModelForNonNewtonFlow/data_pack.hpp"
#include "Public/Shared/NeighborSearchUGB/neighbor_search_ugb.hpp"

namespace SoSim::MMSPH { // device api

    __global__ void
    init_data_cuda(ConstParams *d_cp, DynamicParams *d_dp);

    __global__ void
    precompute_bPart_volume_cuda(ConstParams *d_cp, DynamicParams *d_dp, uint32_t *ns_particleIndices,
                                 uint32_t *ns_neighbors);

    __global__ void
    compute_ext_force_cuda(ConstParams *d_cp, DynamicParams *d_dp, uint32_t *ns_particleIndices,
                           uint32_t *ns_neighbors);

    __global__ void
    advect_pos_cuda(ConstParams *d_cp, DynamicParams *d_dp, uint32_t *ns_particleIndices,
                    uint32_t *ns_neighbors);

    __global__ void
    update_density_cuda(ConstParams *d_cp, DynamicParams *d_dp, uint32_t *ns_particleIndices, uint32_t *ns_neighbors);

    __global__ void
    compute_WC_pressure_force_cuda(ConstParams *d_cp, DynamicParams *d_dp, uint32_t *ns_particleIndices,
                                   uint32_t *ns_neighbors);

    __global__ void
    compute_viscous_force_cuda(ConstParams *d_cp, DynamicParams *d_dp, uint32_t *ns_particleIndices,
                               uint32_t *ns_neighbors);

}

namespace SoSim::MMSPH { // host api

    void
    init_data(ConstParams *d_cp, DynamicParams *d_dp, uint32_t block_num, uint32_t thread_num);

    void
    precompute_bPart_volume(ConstParams *d_cp, DynamicParams *d_dp, uint32_t *ns_particleIndices,
                            uint32_t *ns_neighbors, uint32_t block_num, uint32_t thread_num);

    void
    compute_ext_force(ConstParams *d_cp, DynamicParams *d_dp, uint32_t *ns_particleIndices,
                      uint32_t *ns_neighbors, uint32_t block_num, uint32_t thread_num);

    void
    advect_pos(ConstParams *d_cp, DynamicParams *d_dp, uint32_t *ns_particleIndices,
               uint32_t *ns_neighbors, uint32_t block_num, uint32_t thread_num);

    void
    update_density(ConstParams *d_cp, DynamicParams *d_dp, uint32_t *ns_particleIndices, uint32_t *ns_neighbors,
                   uint32_t block_num, uint32_t thread_num);

    void
    compute_WC_pressure_force(ConstParams *d_cp, DynamicParams *d_dp, uint32_t *ns_particleIndices,
                              uint32_t *ns_neighbors, uint32_t block_num, uint32_t thread_num);

    void
    compute_viscous_force(ConstParams *d_cp, DynamicParams *d_dp, uint32_t *ns_particleIndices,
                          uint32_t *ns_neighbors, uint32_t block_num, uint32_t thread_num);

}

#endif //SOSIM_CUDA_API_CUH
