//
// Created by ADMIN on 2024/3/5.
//

#ifndef SOSIM_UNIFIED_GRID_NS_KERNEL_API_CUH
#define SOSIM_UNIFIED_GRID_NS_KERNEL_API_CUH

#include "libs/NeighborSearchL/unified_grid_ns.hpp"

namespace SoSim {

    extern __host__ void
    resetDevPtr(NeighborSearchUGConfig &h_config, NeighborSearchUGParams &h_params);

    extern __global__ void
    calcParticleHashValue(NeighborSearchUGConfig *d_config, NeighborSearchUGParams *d_params, Vec3f *pos);

    extern __host__ void
    sortByHashValue(NeighborSearchUGConfig &h_config, NeighborSearchUGParams &h_params);

    extern __global__ void
    findCellRange(NeighborSearchUGConfig *d_config, NeighborSearchUGParams *d_params);

    extern __global__ void
    findNeighbors(NeighborSearchUGConfig *d_config, NeighborSearchUGParams *d_params, Vec3f *pos);

    extern __host__ void
    cu_update(NeighborSearchUGConfig &h_config, NeighborSearchUGParams &h_params, NeighborSearchUGConfig *d_config,
              NeighborSearchUGParams *d_params, Vec3f *pos);

}

#endif //SOSIM_UNIFIED_GRID_NS_KERNEL_API_CUH
