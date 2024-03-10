//
// Created by ADMIN on 2024/3/8.
//

#include "dfsph_cuda_api.cuh"

namespace SoSim {

    __device__ inline float
    cubic_value(const Vec3f &r, float h) {
        float r_norm = r.length();
        const float PI = 3.14159265;
        const float cubicSigma = 8.f / PI / static_cast<float>(std::pow(h, 3));

        float res = 0.0;
        float invH = 1 / h;
        float q = r_norm * invH;

        if (q <= 1) {
            if (q <= 0.5) {
                auto q2 = q * q;
                auto q3 = q2 * q;
                res = static_cast<float>(cubicSigma * (6.0 * q3 - 6.0 * q2 + 1));
            } else {
                res = static_cast<float>(cubicSigma * 2 * std::pow(1 - q, 3));
            }
        }

        return res;
    }

    __device__ inline Vec3f
    cubic_gradient(const Vec3f &r, float h) {
        const float PI = 3.14159265;
        const float cubicSigma = 8.f / PI / static_cast<float>(std::pow(h, 3));

        auto res = Vec3f();
        float invH = 1 / h;
        float q = r.length() * invH;

        if (q < 1e-6 || q > 1)
            return res;

        Vec3f grad_q = r / (r.length() * h);
        if (q <= 0.5)
            res = (6 * (3 * q * q - 2 * q)) * grad_q * cubicSigma;
        else {
            auto factor = 1 - q;
            res = -6 * factor * factor * grad_q * cubicSigma;
        }

        return res;
    }

    __global__ void
    computeDensity(DFSPHConstantParams *d_const,
                   DFSPHDynamicParams *d_data,
                   NeighborSearchUGConfig *d_nsConfig,
                   NeighborSearchUGParams *d_nsParams) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->particle_num)
            return;

        auto p_i = d_nsParams->particleIndices_cuData[i];

        if (d_data->mat[p_i] != FLUID)
            return;

        d_data->density[p_i] *= 0.0;
        auto pos_i = d_data->pos_adv[p_i];
        auto neib_ind = p_i * d_nsConfig->maxNeighborNum;

        for (unsigned int p_j = d_nsParams->neighbors_cuData[neib_ind], t = 0;
             p_j != UINT_MAX && t < d_nsConfig->maxNeighborNum;
             ++t, p_j = d_nsParams->neighbors_cuData[neib_ind + t]) {

            auto pos_j = d_data->pos_adv[p_j];
            d_data->density[p_i] += d_data->mass[p_j] * cubic_value(pos_i - pos_j, d_const->h);
        }

        d_data->density[p_i] = max(d_const->rest_density, d_data->density[p_i]);
        d_data->mass[p_i] = d_data->density[p_i] * d_const->rest_volume;
    }

    __global__ void
    computeDFSPHAlpha(DFSPHConstantParams *d_const,
                      DFSPHDynamicParams *d_data,
                      NeighborSearchUGConfig *d_nsConfig,
                      NeighborSearchUGParams *d_nsParams) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->particle_num)
            return;

        auto p_i = d_nsParams->particleIndices_cuData[i];

        if (d_data->mat[p_i] != FLUID)
            return;

        auto pos_i = d_data->pos_adv[p_i];
        auto neib_ind = p_i * d_nsConfig->maxNeighborNum;

        Vec3f alpha_1{0, 0, 0};
        float alpha_2 = 0;
        for (unsigned int p_j = d_nsParams->neighbors_cuData[neib_ind], t = 0;
             p_j != UINT_MAX && t < d_nsConfig->maxNeighborNum;
             ++t, p_j = d_nsParams->neighbors_cuData[neib_ind + t]) {

            auto pos_j = d_data->pos_adv[p_j];

            Vec3f da;
            da = d_data->mass[p_j] * cubic_gradient(pos_i - pos_j, d_const->h);
            if (d_data->mat[p_j] == FIXED_BOUND)
                da = d_data->density[p_i] * d_data->volume[p_j] * cubic_gradient(pos_i - pos_j, d_const->h);
            alpha_1 += da;
            alpha_2 += da.length() * da.length();
        }

        d_data->dfsph_alpha[p_i] = alpha_1.length() * alpha_1.length() + alpha_2 + 1e-6f;
    }

}
