//@author        : Long Shen
//@date          : 2023/11/7
//@description   :
//@version       : 1.0
//@author        : Long Shen
//@date          : 2023/10/26
//@description   :
//@version       : 1.0

#include "Private/MSMSPH/cuda_api.cuh"
#include "Public/Shared/Math/helper_math.hpp"
#include "Public/Shared/SPHKernel/sph_kernel.cuh"

namespace SoSim::MSMSPH { // cuda kernels
    __global__ void
    init_data_cuda(ConstParams *d_const, DynamicParams *d_data) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        d_data->v_mk1[i] *= 0;
        d_data->v_mk2[i] *= 0;
        d_data->acc_m[i] *= 0;
        d_data->M_m[i] *= 0;
        d_data->delta_alpha[i] *= 0;

        switch (d_data->original_phase[i]) {
            case PHASE1:
                d_data->alpha_k[i] = {1, 0};
                d_data->density_m[i] = d_const->rest_density.x;
                break;
            case PHASE2:
                d_data->alpha_k[i] = {0, 1};
                d_data->density_m[i] = d_const->rest_density.y;
                break;
        }

        d_data->pressure_m[i] *= 0;
        d_data->pressure_k[i] *= 0;
        d_data->density_sph[i] *= 0;
    }

    __global__ void
    advect_gravity_cuda(ConstParams *d_const, DynamicParams *d_data) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        d_data->v_m[i] += d_const->gravity * d_const->dt;
        d_data->predictPos[i] += d_data->v_m[i] * d_const->dt;
    }

    __global__ void
    estimate_density_and_pressure_cuda(ConstParams *d_const, DynamicParams *d_data, uint32_t *partIndex,
                                       uint32_t *neighbors) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        auto p_i = partIndex[i];
        auto pos_i = d_data->predictPos[p_i];

        d_data->density_sph[p_i] = 0;
        for (uint32_t nb_0 = p_i * d_const->ns_maxNeighborNum, t = 0, p_j = neighbors[nb_0 + t];
             t < d_const->ns_maxNeighborNum && p_j != UINT_MAX; ++t, p_j = neighbors[nb_0 + t]) {

            auto pos_j = d_data->predictPos[p_j];

            d_data->density_sph[p_i] +=
                    d_data->density_m[p_j] * d_const->rest_volume * cubic_value(pos_i - pos_j, d_const->sph_h);

        }

        float k = 220 * 220 * d_data->density_m[p_i];
        d_data->pressure_m[p_i] =
                k / 7 *
                (pow(d_data->density_sph[p_i] / d_data->density_m[p_i], 7) - 1);
        if (d_data->pressure_m[p_i] < 1e-6)
            d_data->pressure_m[p_i] = 0;

        d_data->pressure_k[p_i] = d_data->pressure_m[p_i] * d_data->alpha_k[p_i];
    }

    __global__ void
    compute_pressure_force_cuda(ConstParams *d_const, DynamicParams *d_data, uint32_t *partIndex,
                                uint32_t *neighbors) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        auto p_i = partIndex[i];
        auto pos_i = d_data->predictPos[p_i];

        d_data->acc_m[p_i] *= 0;

        float3 acc = {0, 0, 0};
        float gamma = -(d_data->alpha_k[p_i].x * d_data->density_m[p_i] / d_const->rest_density.x +
                        d_data->alpha_k[p_i].y * d_data->density_m[p_i] / d_const->rest_density.y);
        for (uint32_t nb_0 = p_i * d_const->ns_maxNeighborNum, t = 0, p_j = neighbors[nb_0 + t];
             t < d_const->ns_maxNeighborNum && p_j != UINT_MAX; ++t, p_j = neighbors[nb_0 + t]) {

            auto pos_j = d_data->predictPos[p_j];

            acc += d_data->density_m[p_j] * d_const->rest_volume *
                   (d_data->pressure_m[p_j] / pow(d_data->density_m[p_j], 2) +
                    d_data->pressure_m[p_i] / pow(d_data->density_m[p_i], 2)) *
                   cubic_grad(pos_i - pos_j, d_const->sph_h);

        }

        d_data->acc_m[p_i] += gamma * acc;

    }


    __global__ void
    advect_pos_cuda(ConstParams *d_const, DynamicParams *d_data) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        if (d_data->alpha_k[i].y == 1)
            return;

        d_data->v_m[i] += d_data->acc_m[i] * d_const->dt;
        d_data->predictPos[i] += d_data->v_m[i] * d_const->dt;
        d_data->pos[i] = d_data->predictPos[i];
    }

}

namespace SoSim::MSMSPH { // host api

    __host__ void
    init_data(ConstParams *d_const, DynamicParams *d_data, uint32_t blockNum, uint32_t threadNum) {
        init_data_cuda<<< blockNum, threadNum >>>(d_const, d_data);
    }

    __host__ void
    advect_gravity(ConstParams *d_const, DynamicParams *d_data, uint32_t blockNum, uint32_t threadNum) {
        advect_gravity_cuda<<< blockNum, threadNum >>>(d_const, d_data);
    }

    __host__ void
    estimate_density_and_pressure(ConstParams *d_const, DynamicParams *d_data, uint32_t *partIndex, uint32_t *neighbors,
                                  uint32_t blockNum, uint32_t threadNum) {
        estimate_density_and_pressure_cuda<<< blockNum, threadNum >>>(d_const, d_data, partIndex, neighbors);
    }

    __host__ void
    compute_pressure_force(ConstParams *d_const, DynamicParams *d_data, uint32_t *partIndex, uint32_t *neighbors,
                           uint32_t blockNum, uint32_t threadNum) {

    }


    __host__ void
    advect_pos(ConstParams *d_const, DynamicParams *d_data, uint32_t blockNum, uint32_t threadNum) {
        advect_pos_cuda<<< blockNum, threadNum >>>(d_const, d_data);
    }

}