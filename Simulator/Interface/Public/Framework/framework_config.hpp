//@author        : Long Shen
//@date          : 2023/11/8
//@description   :
//@version       : 1.0

#ifndef SOSIM_FRAMEWORK_CONFIG_HPP
#define SOSIM_FRAMEWORK_CONFIG_HPP

#include <vector_types.h>
#include <iostream>

namespace SoSim {

    struct CudaConfig {
        uint32_t max_block_num_per_processor{0};
        uint32_t max_thread_num_per_block{0};
    };

    struct SceneConfig {
        uint32_t block_num{0};
        uint32_t thread_num{0};
        float3 scene_lb{0, 0, 0};
        float3 scene_size{0, 0, 0};
    };

    struct SolverConfig {
        uint32_t block_num{0};
        uint32_t thread_num{0};
        float3 scene_lb{0, 0, 0};
        float3 scene_size{0, 0, 0};
        float3 gravity{0.f, -9.8f, 0.0};
        float dt{1 / 60.f};
        float cur_sim_time{0.f};
        float unified_particle_radius{0.f};
    };

    struct ObjectConfig {
        // TODO
    };

}

#endif //SOSIM_FRAMEWORK_CONFIG_HPP
