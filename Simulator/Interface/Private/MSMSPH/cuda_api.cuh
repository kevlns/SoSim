//@author        : Long Shen
//@date          : 2023/10/26
//@description   :
//@version       : 1.0

#ifndef SOSIM_MSMSPH_CUDA_API_CUH
#define SOSIM_MSMSPH_CUDA_API_CUH

#include "data_pack.hpp"
#include "Public/Shared/NeighborSearchUGB/neighbor_search_ugb.hpp"

using NeighborSearcher = SoSim::NSUGB::NeighborSearchUGB;

namespace SoSim::MSMSPH {

    __host__ void
    init_data(ConstParams *d_const, DynamicParams *d_data, uint32_t blockNum, uint32_t threadNum);

    __host__ void
    advect_gravity(ConstParams *d_const, DynamicParams *d_data, uint32_t blockNum, uint32_t threadNum);

    __host__ void
    estimate_density_and_pressure(ConstParams *d_const, DynamicParams *d_data, uint32_t *partIndex, uint32_t *neighbors,
                                  uint32_t blockNum, uint32_t threadNum);

    __host__ void
    compute_pressure_force(ConstParams *d_const, DynamicParams *d_data, uint32_t *partIndex, uint32_t *neighbors,
            uint32_t blockNum, uint32_t threadNum);





    __host__ void
    advect_pos(ConstParams *d_const, DynamicParams *d_data, uint32_t blockNum, uint32_t threadNum);
}


#endif //SOSIM_MSMSPH_CUDA_API_CUH
