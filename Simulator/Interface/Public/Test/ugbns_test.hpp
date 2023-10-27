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

    NSUGB::NeighborSearchUGB neighborSearchUgb;

    float part_radius = 0.05;

    auto cube = generate_cube({-4, -4, -4}, {8, 8, 8}, part_radius);

    neighborSearchUgb.initialize({-20, -20, -20}, {40, 40, 40}, cube.size(), 4 * part_radius);

    float3 *d_pos;
    cudaMalloc((void **) &d_pos, cube.size() * sizeof(float3));
    cudaMemcpy(d_pos, cube.data(), cube.size() * sizeof(float3), cudaMemcpyHostToDevice);
    cudaGetLastError_t("d_pos malloc failed.");

    neighborSearchUgb.update(d_pos);

    neighborSearchUgb.dumpInfo();

    neighborSearchUgb.destroy();

    cudaFree(d_pos);
    d_pos = nullptr;

}

#endif //SOSIM_UGBNS_TEST_HPP
