//@author        : Long Shen
//@date          : 2023/11/8
//@description   :
//@version       : 1.0

#ifndef SOSIM_NEIGHBOR_SEARCH_UGB_CONFIG_HPP
#define SOSIM_NEIGHBOR_SEARCH_UGB_CONFIG_HPP

#include <vector_types.h>
#include <iostream>

#include "Public/Framework/framework_config.hpp"

namespace SoSim::NSUGB {

    struct NeighborSearchConfig {
        uint32_t block_num{0};
        uint32_t thread_num{0};
        float3 scene_lb{0, 0, 0};
        float3 scene_size{0, 0, 0};
        uint32_t max_neighbor_num{0};
        uint32_t total_particle_num{0};
        float particle_radius{0};
    };

}

#endif //SOSIM_NEIGHBOR_SEARCH_UGB_CONFIG_HPP
