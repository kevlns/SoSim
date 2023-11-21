//@author        : Long Shen
//@date          : 2023/11/21
//@description   :
//@version       : 1.0

#ifndef SOSIM_MMSPH_TEST_HPP
#define SOSIM_MMSPH_TEST_HPP

#include "Public/Framework/simulator.hpp"
#include "Public/Framework/scene.hpp"
#include "Public/Framework/solver.hpp"

using namespace SoSim;

void test() {
    CudaConfig cudaConfig{};
    int device;
    cudaGetDevice(&device);
    cudaDeviceProp prop{};
    cudaGetDeviceProperties_v2(&prop, device);
    cudaConfig.max_block_num_per_processor = prop.maxBlocksPerMultiProcessor;
    cudaConfig.max_thread_num_per_block = prop.maxThreadsPerBlock;

    SceneConfig sceneConfig_1{};
    sceneConfig_1.block_num = cudaConfig.max_block_num_per_processor;
    sceneConfig_1.thread_num = cudaConfig.max_thread_num_per_block;
    sceneConfig_1.scene_lb = {-20, -20, -20};
    sceneConfig_1.scene_size = {40, 40, 40};
    auto scene1 = Simulator::createScene(sceneConfig_1);

    SolverConfig solverConfig_1{};
    solverConfig_1.dt = 0.0004;
    solverConfig_1.cur_sim_time = 0.0f;
    solverConfig_1.gravity = {0.f, -9.8f, 0.0};
    auto mmsphSolver = scene1->createSolver(MMSPH_MMSPHSolver, &solverConfig_1);

    mmsphSolver->addParticles(
            "E:\\Projects\\SoSim\\Simulator\\Source\\PhysicalSolvers\\MixtureModelForNonNewtonFlow\\overall_objs.json");

    Simulator::run();

    Simulator::terminate();
}

#endif //SOSIM_MMSPH_TEST_HPP
