//@author        : Long Shen
//@date          : 2023/11/12
//@description   :
//@version       : 1.0

#include "Private/PhysicalSolvers/MixtureModelForNonNewtonFlow/cuda_api.cuh"
#include "Public/Shared/SPHKernel/sph_kernel.cuh"
#include "Public/Shared/Math/helper_math.hpp"

namespace SoSim::MMSPH { // device impl

    __global__ void
    init_data_cuda(ConstParams *d_cp, DynamicParams *d_dp) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_cp->total_particle_num)
            return;

        d_dp->acc_mix[i] *= 0;
        d_dp->drift_vel_k1[i] *= 0;
        d_dp->drift_vel_k2[i] *= 0;
        d_dp->M_m[i] *= 0;
        d_dp->density_mix[i] *= 0;
        d_dp->density_sph[i] *= 0;
        d_dp->pressure_mix[i] *= 0;
        d_dp->bPart_volume[i] *= 0;
    }

    __global__ void
    precompute_bPart_volume_cuda(ConstParams *d_cp, DynamicParams *d_dp, uint32_t *ns_particleIndices,
                                 uint32_t *ns_neighbors) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_cp->total_particle_num)
            return;

        auto p_i = ns_particleIndices[i];

        if (d_dp->mat[p_i] != Material::BOUND)
            return;

        auto pos_i = d_dp->predictPos[p_i];
        float delta = 1e-6;
        for (uint32_t neib_0 = p_i * d_cp->max_neighbor_num, p_j = ns_neighbors[neib_0], t = 0;
             p_j != UINT_MAX && t < d_cp->max_neighbor_num; ++t, p_j = ns_neighbors[neib_0 + t]) {

            if (d_dp->mat[p_j] == BOUND) {
                auto pos_j = d_dp->predictPos[p_j];
                delta += cubic_value(pos_i - pos_j, d_cp->sph_h);
            }

        }

        d_dp->bPart_volume[p_i] = 1 / delta;

    }

    __global__ void
    compute_ext_force_cuda(ConstParams *d_cp, DynamicParams *d_dp, uint32_t *ns_particleIndices,
                           uint32_t *ns_neighbors) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_cp->total_particle_num)
            return;

        if (d_dp->mat[i] == BOUND)
            return;

        // gravity
//        d_dp->acc_mix[i] = d_cp->gravity;
        d_dp->acc_mix[i] += d_cp->gravity;

    }

    __global__ void
    advect_pos_cuda(ConstParams *d_cp, DynamicParams *d_dp, uint32_t *ns_particleIndices,
                    uint32_t *ns_neighbors) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_cp->total_particle_num)
            return;

        if (d_dp->mat[i] == BOUND)
            return;

        d_dp->vel_mix[i] += d_dp->acc_mix[i] * d_cp->dt;
        d_dp->predictPos[i] += d_dp->vel_mix[i] * d_cp->dt;
        d_dp->pos[i] = d_dp->predictPos[i];

    }

    __global__ void
    update_density_cuda(ConstParams *d_cp, DynamicParams *d_dp, uint32_t *ns_particleIndices, uint32_t *ns_neighbors) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_cp->total_particle_num)
            return;

        auto p_i = ns_particleIndices[i];

        if (d_dp->mat[p_i] == BOUND)
            return;

        // mixture density
        d_dp->density_mix[p_i] = dot(d_cp->rest_density, d_dp->alpha[p_i]);
        __syncthreads();

        // sph density
        d_dp->density_sph[p_i] *= 0;
        auto pos_i = d_dp->predictPos[p_i];
        for (uint32_t neib_0 = p_i * d_cp->max_neighbor_num, p_j = ns_neighbors[neib_0], t = 0;
             p_j != UINT_MAX && t < d_cp->max_neighbor_num; ++t, p_j = ns_neighbors[neib_0 + t]) {

            auto pos_j = d_dp->predictPos[p_j];

            if (length(pos_i - pos_j) < 1e-6)
                continue;

            if (d_dp->mat[p_j] == FLUID) {
                d_dp->density_sph[p_i] +=
                        d_dp->density_mix[p_j] * d_cp->rest_volume * cubic_value(pos_i - pos_j, d_cp->sph_h);
            } else if (d_dp->mat[p_j] == BOUND) {
                d_dp->density_sph[p_i] +=
                        d_dp->density_mix[p_i] * d_dp->bPart_volume[p_j] * cubic_value(pos_i - pos_j, d_cp->sph_h);
            }

        }

        d_dp->density_sph[p_i] = max(d_dp->density_sph[p_i], d_dp->density_mix[p_i]);

    }

    __global__ void
    compute_WC_pressure_force_cuda(ConstParams *d_cp, DynamicParams *d_dp, uint32_t *ns_particleIndices,
                                   uint32_t *ns_neighbors) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_cp->total_particle_num)
            return;

        auto p_i = ns_particleIndices[i];

        if (d_dp->mat[p_i] == BOUND)
            return;

        // reset acc
        d_dp->acc_mix[p_i] *= 0;

        // compute pressure
        float k = 10000;  // 220*220 = 48400
        d_dp->pressure_mix[p_i] =
                k * d_dp->density_mix[p_i] / 7 * (pow(d_dp->density_sph[p_i] / d_dp->density_mix[p_i], 7) - 1);
        if (d_dp->pressure_mix[p_i] < 0)
            d_dp->pressure_mix[p_i] = 0;
        __syncthreads();

        // compute pressure force
        auto pos_i = d_dp->predictPos[p_i];
        float dpi = d_dp->pressure_mix[p_i] / (pow(d_dp->density_sph[p_i], 2));
        float3 acc = {0, 0, 0};
        for (uint32_t neib_0 = p_i * d_cp->max_neighbor_num, p_j = ns_neighbors[neib_0], t = 0;
             p_j != UINT_MAX && t < d_cp->max_neighbor_num; ++t, p_j = ns_neighbors[neib_0 + t]) {

            auto pos_j = d_dp->predictPos[p_j];

            if (d_dp->mat[p_j] == FLUID) {
                float dpj = d_dp->pressure_mix[p_j] / (pow(d_dp->density_sph[p_j], 2));
                acc += -d_dp->density_mix[p_j] * d_cp->rest_volume * (dpi + dpj) *
                       cubic_grad(pos_i - pos_j, d_cp->sph_h); // TODO 计算质量的时候使用sph密度还是静止密度？
            } else if (d_dp->mat[p_j] == BOUND) {
                acc += -d_dp->density_mix[p_i] * d_dp->bPart_volume[p_j] * (dpi + dpi) *
                       cubic_grad(pos_i - pos_j, d_cp->sph_h);
            }
        }

        d_dp->acc_mix[p_i] += acc / d_dp->density_mix[p_i];

    }

    __global__ void
    compute_viscous_force_cuda(ConstParams *d_cp, DynamicParams *d_dp, uint32_t *ns_particleIndices,
                               uint32_t *ns_neighbors) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_cp->total_particle_num)
            return;

        auto p_i = ns_particleIndices[i];

        if (d_dp->mat[p_i] == BOUND)
            return;

        float niu = 0.001; // rest viscosity

        // compute viscous force
        auto pos_i = d_dp->predictPos[p_i];
        auto vel_i = d_dp->vel_mix[p_i];
        float3 acc = {0, 0, 0};
        for (uint32_t neib_0 = p_i * d_cp->max_neighbor_num, p_j = ns_neighbors[neib_0], t = 0;
             p_j != UINT_MAX && t < d_cp->max_neighbor_num; ++t, p_j = ns_neighbors[neib_0 + t]) {

            auto pos_j = d_dp->predictPos[p_j];
            auto vel_j = d_dp->vel_mix[p_j];

            if (d_dp->mat[p_j] == FLUID) {
                acc += 10 * d_cp->rest_volume * niu *
                       dot(vel_i - vel_j, pos_i - pos_j) /
                       (0.01 * d_cp->sph_h * d_cp->sph_h + pow(length(pos_j - pos_i), 2)) *
                       cubic_grad(pos_i - pos_j, d_cp->sph_h);
            }
        }

        d_dp->acc_mix[p_i] += acc / d_dp->density_mix[p_i];
    }

}

namespace SoSim::MMSPH { // host impl

    void
    init_data(ConstParams *d_cp, DynamicParams *d_dp, uint32_t block_num, uint32_t thread_num) {
        init_data_cuda<<<block_num, thread_num>>>(d_cp, d_dp);
    }

    void
    precompute_bPart_volume(ConstParams *d_cp, DynamicParams *d_dp, uint32_t *ns_particleIndices,
                            uint32_t *ns_neighbors, uint32_t block_num, uint32_t thread_num) {
        precompute_bPart_volume_cuda<<<block_num, thread_num>>>(d_cp, d_dp, ns_particleIndices, ns_neighbors);
    }

    void
    compute_ext_force(ConstParams *d_cp, DynamicParams *d_dp, uint32_t *ns_particleIndices,
                      uint32_t *ns_neighbors, uint32_t block_num, uint32_t thread_num) {
        compute_ext_force_cuda<<<block_num, thread_num>>>(d_cp, d_dp, ns_particleIndices, ns_neighbors);
    }

    void
    advect_pos(ConstParams *d_cp, DynamicParams *d_dp, uint32_t *ns_particleIndices,
               uint32_t *ns_neighbors, uint32_t block_num, uint32_t thread_num) {
        advect_pos_cuda<<<block_num, thread_num>>>(d_cp, d_dp, ns_particleIndices, ns_neighbors);
    }

    void
    update_density(ConstParams *d_cp, DynamicParams *d_dp, uint32_t *ns_particleIndices, uint32_t *ns_neighbors,
                   uint32_t block_num, uint32_t thread_num) {
        update_density_cuda<<<block_num, thread_num>>>(d_cp, d_dp, ns_particleIndices, ns_neighbors);
    }

    void
    compute_WC_pressure_force(ConstParams *d_cp, DynamicParams *d_dp, uint32_t *ns_particleIndices,
                              uint32_t *ns_neighbors, uint32_t block_num, uint32_t thread_num) {
        compute_WC_pressure_force_cuda<<<block_num, thread_num>>>(d_cp, d_dp, ns_particleIndices, ns_neighbors);
    }

    void
    compute_viscous_force(ConstParams *d_cp, DynamicParams *d_dp, uint32_t *ns_particleIndices,
                          uint32_t *ns_neighbors, uint32_t block_num, uint32_t thread_num) {
        compute_viscous_force_cuda<<<block_num, thread_num>>>(d_cp, d_dp, ns_particleIndices, ns_neighbors);
    }

}



