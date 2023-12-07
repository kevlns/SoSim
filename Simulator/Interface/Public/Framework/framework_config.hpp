//@author        : Long Shen
//@date          : 2023/11/8
//@description   :
//@version       : 1.0

#ifndef SOSIM_FRAMEWORK_CONFIG_HPP
#define SOSIM_FRAMEWORK_CONFIG_HPP

#include <vector_types.h>
#include <vector>
#include <iostream>
#include <set>

#include "Public/Framework/enum.hpp"
#include "Public/PhysicalSolvers/solver_ref.hpp"

namespace SoSim {

    struct CudaConfig {
        uint32_t max_block_num_per_processor{0};
        uint32_t max_thread_num_per_block{0};
        uint32_t max_multiprocessor{0};
    };

    struct TransformConfig {
        float3 transform{0, 0, 0};
        float3 rotate{0, 0, 0};
        float3 scale{1, 1, 1};
    };

    struct SceneConfig {
        uint32_t max_block_num_per_processor{0};
        uint32_t max_thread_num_per_block{0};
        uint32_t max_multiprocessor{0};
        float3 scene_lb{-15, -15, -15};
        float3 scene_size{30, 30, 30};

        // for gui
        bool preview{false};
    };

    struct SolverConfig {
        uint32_t max_block_num_per_processor{0};
        uint32_t max_thread_num_per_block{0};
        uint32_t max_multiprocessor{0};
        dim3 block_num{0};
        dim3 thread_num{0};
        float3 scene_lb{0, 0, 0};
        float3 scene_size{0, 0, 0};
        float3 gravity{0.f, -9.8f, 0.0};
        float dt{1 / 60.f};
        float cur_sim_time{0.f};
        float unified_particle_radius{0.f};
        SolverInstanceType instanceType{SOLVER_NONE};

        // for gui
        int solver_index{0};

        inline void computeCudaBlockNum(uint32_t particle_num) {
            block_num = dim3(std::min((particle_num + max_thread_num_per_block - 1) / max_thread_num_per_block,
                                      max_multiprocessor * max_block_num_per_processor));
        }
    };

    struct ObjectConfig {
        std::string model_file_path{"None"};
        std::string shape;
        MaterialType mat;
        PhaseType phase;
        std::set<ComponentType> components;
        float3 lb{0, 0, 0};
        float3 size{0, 0, 0};
        float3 color{1, 1, 1};
        float height{0};
        float3 top_center{0, 0, 0};
        float area_radius{0};
        int particle_layer_num{2};
        float particle_radius{0};
        TransformConfig t_config{};

        // for gui
        char model_file_path_cStr[256]{""};
        int shape_index{0};
        int mat_index{0};
        int phase_index{0};
        bool preview{false};
        std::vector<bool> component_selected;
    };

    struct ParticleObjectConfig {
        std::string model_file_path;
        std::string shape;
        float3 lb{0, 0, 0};
        float3 size{0, 0, 0};
        float3 color{1, 1, 1};
        float height{0};
        float3 top_center{0, 0, 0};
        float area_radius{0};
        float particle_radius{0};
        TransformConfig t_config{};
    };
}

#endif //SOSIM_FRAMEWORK_CONFIG_HPP
