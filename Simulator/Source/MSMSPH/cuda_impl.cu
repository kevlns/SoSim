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

        switch (d_data->original_phase[i]) {
            case Phase::PHASE1:
                d_data->volF_k[i] = {1, 0};
                break;
            case Phase::PHASE2:
                d_data->volF_k[i] = {0, 1};
                break;
            case Phase::PHASE3:
                break;
        }

        d_data->predictPos[i] = d_data->pos[i];
        d_data->v_mk1[i] *= 0;
        d_data->v_mk2[i] *= 0;
        d_data->acc[i] *= 0;
        d_data->d_volF[i] *= 0;
        d_data->pressure_k[i] *= 0;
        d_data->pressure[i] *= 0;
        d_data->mass[i] *= 0;
        d_data->density[i] *= 0;
        d_data->density_ba[i] *= 0;

    }

    __global__ void
    compute_drift_vel_cuda(ConstParams *d_const, DynamicParams *d_data, uint32_t *partIndex, uint32_t *neighbors) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

        if (i >= d_const->total_particle_num)
            return;

        uint32_t p_i = partIndex[i];

        if (d_data->mat[p_i] == Material::BOUND)
            return;

        auto pos_i = d_data->predictPos[p_i];
        float tao = 1e-7;
        float sigma = 1e-4;

        // compute gradient of pressure and volume fraction
        float3 dp_k1 = {0, 0, 0};
        float3 dp_k2 = {0, 0, 0};
        float3 da_k1 = {0, 0, 0};
        float3 da_k2 = {0, 0, 0};

        for (uint32_t nb_0 = p_i * d_const->ns_maxNeighborNum, t = 0, p_j = neighbors[nb_0 + t];
             t < d_const->ns_maxNeighborNum && neighbors[p_j] != UINT_MAX; ++t) {

            auto pos_j = d_data->predictPos[p_j];

            if (d_data->mat[p_j] == Material::BOUND)
                continue;

            dp_k1 += d_data->density[p_j] * d_const->rest_volume /
                     d_data->density_ba[p_j] * (d_data->pressure_k[p_j].x - d_data->pressure_k[p_i].x) *
                     cubic_grad(pos_i - pos_j, d_const->sph_h);

            dp_k2 += d_data->density[p_j] * d_const->rest_volume /
                     d_data->density_ba[p_j] * (d_data->pressure_k[p_j].y - d_data->pressure_k[p_i].y) *
                     cubic_grad(pos_i - pos_j, d_const->sph_h);

            da_k1 += d_data->density[p_j] * d_const->rest_volume /
                     d_data->density_ba[p_j] * (d_data->volF_k[p_j].x - d_data->volF_k[p_i].x) *
                     cubic_grad(pos_i - pos_j, d_const->sph_h);

            da_k2 += d_data->density[p_j] * d_const->rest_volume /
                     d_data->density_ba[p_j] * (d_data->volF_k[p_j].y - d_data->volF_k[p_i].y) *
                     cubic_grad(pos_i - pos_j, d_const->sph_h);
        }

        // 3 sum-items in eq(9)
        float factor1_1 = 0;
        float factor1_2 = 0;
        float3 factor2_1 = {0, 0, 0};
        float3 factor2_2 = {0, 0, 0};
        float3 factor3_1 = {0, 0, 0};
        float3 factor3_2 = {0, 0, 0};

        for (uint32_t nb_0 = p_i * d_const->ns_maxNeighborNum, t = 0, p_j = neighbors[nb_0 + t];
             t < d_const->ns_maxNeighborNum && neighbors[p_j] != UINT_MAX; ++t) {

            auto pos_j = d_data->predictPos[p_j];

            float ck1 = d_data->volF_k[p_i].x * d_const->rest_density.x / d_data->density[p_i];
            float ck2 = d_data->volF_k[p_i].y * d_const->rest_density.y / d_data->density[p_i];

            factor1_1 += ck1 * d_const->rest_density.x;
            factor1_2 += ck2 * d_const->rest_density.y;

            factor2_1 += ck1 * dp_k1;
            factor2_2 += ck2 * dp_k2;

            factor3_1 += ck1 * da_k1 / d_data->volF_k[p_i].x;
            factor3_2 += ck2 * da_k2 / d_data->volF_k[p_i].y;
        }

        // compute acc, eq(10)
        float3 acc = {0, 0, 0};
        float3 dp_m = {0, 0, 0};
        float niu = 0.001;

        for (uint32_t nb_0 = p_i * d_const->ns_maxNeighborNum, t = 0, p_j = neighbors[nb_0 + t];
             t < d_const->ns_maxNeighborNum && neighbors[p_j] != UINT_MAX; ++t) {

            auto pos_j = d_data->predictPos[p_j];
            auto pos_ij = pos_i - pos_j;
            auto gradW = cubic_grad(pos_ij, d_const->sph_h);

            dp_m += d_data->density[p_j] * d_const->rest_volume * (d_data->pressure[p_i] + d_data->pressure[p_j]) / 2 *
                    gradW;

            acc += d_data->density[p_j] * d_const->rest_volume / d_data->density_ba[p_j] * niu *
                   (d_data->v_m[p_j] - d_data->v_m[p_i]) *
                   dot(-pos_ij, gradW) / (length(pos_ij) * length(pos_ij));
        }
        acc = dp_m / d_data->density[p_i] - acc;

        // compute v_k, eq(9)
        d_data->v_mk1[p_i] = tao * (d_const->rest_density.x - factor1_1) * acc - tao * (dp_k1 - factor2_1) -
                             sigma * (da_k1 / d_data->volF_k[p_i].x - factor3_1);
        d_data->v_mk2[p_i] = tao * (d_const->rest_density.y - factor1_2) * acc - tao * (dp_k2 - factor2_2) -
                             sigma * (da_k2 / d_data->volF_k[p_i].y - factor3_2);
    }

    __global__ void
    advect_volFrac_cuda(ConstParams *d_const, DynamicParams *d_data, uint32_t *partIndex, uint32_t *neighbors) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

        if (i >= d_const->total_particle_num)
            return;

        uint32_t p_i = partIndex[i];

        if (d_data->mat[p_i] == Material::BOUND)
            return;

        auto pos_i = d_data->predictPos[p_i];
        d_data->d_volF[p_i] *= 0;
        float a1 = 0;
        float a2 = 0;

        // all p_i neighbors
        float den_ = 0;
        for (uint32_t nb_0 = p_i * d_const->ns_maxNeighborNum, t = 0, p_j = neighbors[nb_0 + t];
             t < d_const->ns_maxNeighborNum && neighbors[p_j] != UINT_MAX; ++t) {

            auto pos_j = d_data->predictPos[p_j];

            float3 gradW = cubic_grad(pos_i - pos_j, d_const->sph_h);
            float f1 = d_data->density[p_j] / d_data->density_ba[p_j] * d_const->rest_volume;
            float f2 = dot(d_data->v_m[p_j] - d_data->v_m[p_i], gradW);


            // for phase1
            a1 += f1 * (d_data->volF_k[p_j].x + d_data->volF_k[p_i].x) / 2 * f2;
            a1 += f1 *
                  dot(d_data->volF_k[p_j].x * d_data->v_mk1[p_j] + d_data->volF_k[p_i].x * d_data->v_mk1[p_i], gradW);

            // for phase2
            a2 += f1 * (d_data->volF_k[p_j].y + d_data->volF_k[p_i].y) / 2 * f2;
            a2 += f1 *
                  dot(d_data->volF_k[p_j].y * d_data->v_mk2[p_j] + d_data->volF_k[p_i].y * d_data->v_mk2[p_i], gradW);
        }

        d_data->d_volF[p_i] = {-a1, -a2};

        __syncthreads();

        d_data->volF_k[p_i] += d_data->d_volF[p_i];
    }

    __global__ void
    update_density_and_pressure_cuda(ConstParams *d_const, DynamicParams *d_data, uint32_t *partIndex,
                                     uint32_t *neighbors) {

        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

        if (i >= d_const->total_particle_num)
            return;

        uint32_t p_i = partIndex[i];

        if (d_data->mat[p_i] == Material::BOUND)
            return;

        auto pos_i = d_data->predictPos[p_i];

        // all p_i neighbors
        d_data->density_ba[p_i] *= 0;
        for (uint32_t nb_0 = p_i * d_const->ns_maxNeighborNum, t = 0, p_j = neighbors[nb_0 + t];
             t < d_const->ns_maxNeighborNum && neighbors[p_j] != UINT_MAX; ++t) {

            auto pos_j = d_data->predictPos[p_j];
            d_data->density_ba[p_i] +=
                    d_data->density[p_j] * d_const->rest_volume * cubic_value(pos_i - pos_j, d_const->sph_h);

        }

        d_data->density[p_i] = dot(d_data->volF_k[p_i], d_const->rest_density);
        d_data->mass[p_i] = d_data->density[p_i] * d_const->rest_volume;
        d_data->density_ba[p_i] = max(d_data->density[p_i], d_data->density_ba[p_i]);
        float k = 220 * 220;
        d_data->pressure[p_i] =
                k * d_data->density[p_i] / 7 * (pow(d_data->density_ba[p_i] / d_data->density[p_i], 7) - 1);
        d_data->pressure_k[p_i] = d_data->pressure[p_i] * d_data->volF_k[p_i];
    }

    // not used yet; not finished yet
    __global__ void
    correct_pressure_cuda(ConstParams *d_const, DynamicParams *d_data, uint32_t *partIndex, uint32_t *neighbors) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

        if (i >= d_const->total_particle_num)
            return;

        uint32_t p_i = partIndex[i];

        if (d_data->mat[p_i] == Material::BOUND)
            return;

        // bound volFrac
        d_data->volF_k[p_i] += d_data->d_volF[p_i];
        if (d_data->volF_k[p_i].x < 1e-6)
            d_data->volF_k[p_i].x = 0;
        if (d_data->volF_k[p_i].y < 1e-6)
            d_data->volF_k[p_i].y = 0;
        d_data->volF_k[p_i].x /= dot({1, 1}, d_data->volF_k[p_i]); // re-scale
        d_data->volF_k[p_i].y /= dot({1, 1}, d_data->volF_k[p_i]); // re-scale

        __syncthreads();

        float k = 220 * 220;
        float d1 = 0;
        float d2 = 0;

        for (uint32_t nb_0 = p_i * d_const->ns_maxNeighborNum, t = 0, p_j = neighbors[nb_0 + t];
             t < d_const->ns_maxNeighborNum && neighbors[p_j] != UINT_MAX; ++t) {

            // for phase1
            d1 += -k * d_const->rest_density.x / 7 * (6 * pow(d_data->density_ba[p_i] / d_data->density[p_i], 7) + 1) *
                  d_data->d_volF[p_i].x;

            // for phase2
            d2 += -k * d_const->rest_density.y / 7 * (6 * pow(d_data->density_ba[p_i] / d_data->density[p_i], 7) + 1) *
                  d_data->d_volF[p_i].y;

        }

    }

    __global__ void
    compute_overall_acc_cuda(ConstParams *d_const, DynamicParams *d_data, uint32_t *partIndex, uint32_t *neighbors) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

        if (i >= d_const->total_particle_num)
            return;

        uint32_t p_i = partIndex[i];

        if (d_data->mat[p_i] == Material::BOUND)
            return;

        auto pos_i = d_data->predictPos[p_i];

        // compute acc, eq(10)
        float3 vis_acc = {0, 0, 0};
        float3 dp_m = {0, 0, 0};
        float niu = 0.001;

        for (uint32_t nb_0 = p_i * d_const->ns_maxNeighborNum, t = 0, p_j = neighbors[nb_0 + t];
             t < d_const->ns_maxNeighborNum && neighbors[p_j] != UINT_MAX; ++t) {

            auto pos_j = d_data->predictPos[p_j];
            auto pos_ij = pos_i - pos_j;
            auto gradW = cubic_grad(pos_ij, d_const->sph_h);

            if (d_data->mat[p_j] != Material::BOUND) {

                if (length(pos_ij) < 1e-6)
                    continue;

//                dp_m += d_data->density[p_j] * d_const->rest_volume * (d_data->pressure[p_i] + d_data->pressure[p_j]) /
//                        (2 * d_data->density_ba[p_j]) * gradW;

                vis_acc += d_data->density[p_j] * d_const->rest_volume / d_data->density_ba[p_j] * niu *
                           (d_data->v_m[p_j] - d_data->v_m[p_i]) *
                           dot(-pos_ij, gradW) / pow(length(pos_ij), 2);
            }
        }

        d_data->acc[p_i] = -dp_m / d_data->density[p_i] + vis_acc / d_data->density[p_i] + d_const->gravity;

    }

    __global__ void
    advect_particles_cuda(ConstParams *d_const, DynamicParams *d_data, uint32_t *partIndex) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

        if (i >= d_const->total_particle_num)
            return;

        uint32_t p_i = partIndex[i];

        if (d_data->mat[p_i] == Material::BOUND)
            return;

        d_data->v_m[p_i] += d_data->acc[p_i] * d_const->dt;
        d_data->predictPos[p_i] += d_data->v_m[p_i] * d_const->dt;
        d_data->pos[p_i] = d_data->predictPos[p_i];
    }

}

namespace SoSim::MSMSPH { // host invoke

    __host__ void
    init_data(ConstParams *d_const, DynamicParams *d_data, uint32_t blockNum, uint32_t threadNum) {
        init_data_cuda<<< blockNum, threadNum>>>(d_const, d_data);
    }

    __host__ void
    compute_drift_vel(ConstParams *d_const, DynamicParams *d_data, uint32_t *partIndex, uint32_t *neighbors,
                      uint32_t blockNum, uint32_t threadNum) {
        compute_drift_vel_cuda<<< blockNum, threadNum>>>(d_const, d_data, partIndex, neighbors);
    }

    __host__ void
    advect_volFrac(ConstParams *d_const, DynamicParams *d_data, uint32_t *partIndex, uint32_t *neighbors,
                   uint32_t blockNum, uint32_t threadNum) {
        advect_volFrac_cuda<<< blockNum, threadNum>>>(d_const, d_data, partIndex, neighbors);
    }

    __host__ void
    update_density_and_pressure(ConstParams *d_const, DynamicParams *d_data, uint32_t *partIndex, uint32_t *neighbors,
                                uint32_t blockNum, uint32_t threadNum) {
        update_density_and_pressure_cuda<<< blockNum, threadNum>>>(d_const, d_data, partIndex, neighbors);
    }

    __host__ void
    compute_overall_acc(ConstParams *d_const, DynamicParams *d_data, uint32_t *partIndex, uint32_t *neighbors,
                        uint32_t blockNum, uint32_t threadNum) {
        compute_overall_acc_cuda<<< blockNum, threadNum>>>(d_const, d_data, partIndex, neighbors);
    }

    __host__ void
    advect_particles(ConstParams *d_const, DynamicParams *d_data, uint32_t *partIndex, uint32_t blockNum,
                     uint32_t threadNum) {
        advect_particles_cuda<<< blockNum, threadNum>>>(d_const, d_data, partIndex);
    }

}