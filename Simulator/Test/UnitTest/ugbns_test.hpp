//@author        : Long Shen
//@date          : 2023/10/27
//@description   :
//@version       : 1.0

#ifndef SOSIM_UGBNS_TEST_HPP
#define SOSIM_UGBNS_TEST_HPP

#include "Public/Shared/NeighborSearchUGB/neighbor_search_ugb.hpp"
#include "Public/Shared/ModelUtils/model_tool.hpp"
#include "Public/Shared/CudaUtils/cuda_tool.hpp"

using namespace SoSim;

void UGBNS_Test() {

    NSUGB::NeighborSearcher neighborSearchUgb;

    float part_radius = 0.05;

    auto cube = generate_cube({-4, -4, -4}, {8, 8, 8}, part_radius);

    float3 *d_pos;
    cudaMalloc((void **) &d_pos, cube.size() * sizeof(float3));
    cudaMemcpy(d_pos, cube.data(), cube.size() * sizeof(float3), cudaMemcpyHostToDevice);
    cudaGetLastError_t("d_pos malloc failed.");

    NSUGB::NeighborSearchConfig ns_config;
    ns_config.block_num = 16;
    ns_config.thread_num = 1024;
    ns_config.scene_lb = {-20, -20, -20};
    ns_config.scene_size = {40, 40, 40};
    ns_config.max_neighbor_num = 35;
    ns_config.total_particle_num = cube.size();
    ns_config.particle_radius = part_radius;
    neighborSearchUgb.setConfig(&ns_config);

    neighborSearchUgb.update(d_pos);

    neighborSearchUgb.dumpInfo();

    neighborSearchUgb.destroy();

    cudaFree(d_pos);
    d_pos = nullptr;

}

#endif //SOSIM_UGBNS_TEST_HPP
