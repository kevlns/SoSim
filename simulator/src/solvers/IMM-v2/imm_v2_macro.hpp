//
// Created by ADMIN on 2024/3/15.
//

#ifndef SOSIM_IMM_MACRO_HPP
#define SOSIM_IMM_MACRO_HPP

#define CHECK_THREAD() \
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x; \
    if (i >= d_const->particle_num)                     \
        return;                                         \
    auto p_i = d_nsParams->particleIndices_cuData[i];

#define FOR_EACH_NEIGHBOR_Pj() \
       auto neib_ind = p_i * d_nsConfig->maxNeighborNum;                        \
       for (unsigned int p_j = d_nsParams->neighbors_cuData[neib_ind], t = 0;   \
            p_j != UINT_MAX && t < d_nsConfig->maxNeighborNum;                  \
            ++t, p_j = d_nsParams->neighbors_cuData[neib_ind + t])

#define CONST_VALUE(name) \
        d_const->name

#define CONST_VALUE_PHASE(name, k) \
        d_data->name[k]

#define DATA_VALUE(name, index) \
        d_data->name[index]

#define DATA_VALUE_PHASE(name, index, phase_num, k) \
        d_data->name[index*phase_num + k]

#define FOR_EACH_PHASE_k() \
       for (int k=0; k < d_const->phase_num; ++k)     \

#define CUBIC_KERNEL_VALUE() \
        cubic_value(pos_i - pos_j, d_const->sph_h)

#define CUBIC_KERNEL_GRAD() \
        cubic_gradient(pos_i - pos_j, d_const->sph_h)

#endif //SOSIM_IMM_MACRO_HPP
