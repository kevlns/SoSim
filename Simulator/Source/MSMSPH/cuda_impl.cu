//@author        : Long Shen
//@date          : 2023/10/26
//@description   :
//@version       : 1.0

#include "Private/MSMSPH/cuda_api.cuh"
#include "Public/Shared/Math/helper_math.hpp"
#include "Public/Shared/SPHKernel/sph_kernel.cuh"

namespace SoSim::MSMSPH {

    extern __device__ void
    update_density(ConstParams *d_const, DynamicParams *d_data, uint32_t *f_partIndex, uint32_t *f_neighbors,
                   uint32_t *b_partIndex, uint32_t *b_neighbors, uint32_t p_i) {

        d_data->density[p_i] = dot(d_const->rest_density, d_data->volF_k[p_i]);
        d_data->mass[p_i] = d_data->density[p_i] * d_const->rest_volume;
    }

    extern __device__ void
    update_pressure(ConstParams *d_const, DynamicParams *d_data, uint32_t *f_partIndex, uint32_t *f_neighbors,
                    uint32_t *b_partIndex, uint32_t *b_neighbors, uint32_t p_i) {

        uint32_t neighb_ind0 = p_i * d_const->ns_maxNeighborNum;

        // compute fluid particles' contribution
        const float k = 1000 / d_data->density[p_i] * 220 * 220 * 1000;
        const float eta = 7;
        float den_ = 0.0; // density_ba
        auto x_i = d_data->predictPos[p_i];
        for (unsigned int p_j = f_neighbors[neighb_ind0], t = 0;
             p_j != UINT_MAX && t < d_const->ns_maxNeighborNum;
             ++t, p_j = f_neighbors[neighb_ind0 + t]) {

            auto x_j = d_data->predictPos[p_j];
            auto x_ij = x_i - x_j;

            den_ += d_data->mass[p_j] * poly6_value(x_ij, d_const->sph_h);
        }

        // TODO compute rigid particles' contribution
        // ...

        d_data->pressure[p_i] = k * (den_ - d_data->density[p_i]); // eq(26)

        // correct pressure
        float factor = -k / eta * (eta - 1) * (pow(den_ / d_data->density[p_i], eta) + 1);
        d_data->pressure[p_i] += dot(factor * d_const->rest_density, d_data->d_volF[p_i]); // eq(30,31)

    }

    extern __global__ void
    __step__(ConstParams *d_const, DynamicParams *d_data, uint32_t *f_partIndex, uint32_t *f_neighbors,
             uint32_t *b_partIndex, uint32_t *b_neighbors) {

        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;

        if (i >= d_const->total_particle_num)
            return;

        uint32_t p_i = f_partIndex[i];

        update_density(d_const, d_data, f_partIndex, f_neighbors, b_partIndex, b_neighbors, p_i);
        __syncthreads();

        update_pressure(d_const, d_data, f_partIndex, f_neighbors, b_partIndex, b_neighbors, p_i);
        __syncthreads();

    }

}