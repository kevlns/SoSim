//@author        : Long Shen
//@date          : 2023/10/26
//@description   :
//@version       : 1.0

#ifndef SOSIM_MSMSPH_CUDA_API_CUH
#define SOSIM_MSMSPH_CUDA_API_CUH

#include "data_pack.hpp"

namespace SoSim::MSMSPH {

    extern __device__ void
    update_density(ConstParams *d_const, DynamicParams *d_data, uint32_t *f_partIndex, uint32_t *f_neighbors,
                   uint32_t *b_partIndex, uint32_t *b_neighbors, uint32_t p_i);

    extern __device__ void
    update_pressure(ConstParams *d_const, DynamicParams *d_data, uint32_t *f_partIndex, uint32_t *f_neighbors,
                    uint32_t *b_partIndex, uint32_t *b_neighbors, uint32_t p_i);

    extern __device__ void update_vel();

    extern __device__ void correct_pressure();

    extern __device__ void update_vm();

    extern __device__ void update_vk();

    extern __device__ void update_pos();

    extern __global__ void
    __step__(ConstParams *d_const, DynamicParams *d_data, uint32_t *f_partIndex, uint32_t *f_neighbors,
             uint32_t *b_partIndex, uint32_t *b_neighbors);

    __host__ void
    step_cuda(ConstParams *d_const, DynamicParams *d_data, uint32_t *f_partIndex, uint32_t *f_neighbors,
              uint32_t *b_partIndex, uint32_t *b_neighbors, uint32_t blockNum, uint32_t threadNum) {
        __step__<<< blockNum, threadNum >>>(d_const, d_data, f_partIndex, f_neighbors, b_partIndex, b_neighbors);
    }

}


#endif //SOSIM_MSMSPH_CUDA_API_CUH
