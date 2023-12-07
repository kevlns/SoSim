//@author        : Long Shen
//@date          : 2023/11/9
//@description   :
//@version       : 1.0

#include <cuda_runtime.h>

#include "Public/Framework/simulator.hpp"

namespace SoSim {

    Simulator::Simulator() {
        m_cudaConfig = new CudaConfig;

        int dev;
        cudaGetDevice(&dev);

        cudaDeviceProp prop{};
        cudaGetDeviceProperties_v2(&prop, dev);

        m_cudaConfig->max_block_num_per_processor = prop.maxBlocksPerMultiProcessor;
        m_cudaConfig->max_thread_num_per_block = prop.maxThreadsPerBlock;
        m_cudaConfig->max_multiprocessor = prop.multiProcessorCount;
    }

    Simulator::~Simulator() {
        delete m_cudaConfig;
        m_cudaConfig = nullptr;
    }

    void Simulator::addSceneDefault() {
        auto scene = new Scene(m_cudaConfig);
        auto sceneConfig = scene->getConfig();
        m_scene_map[sceneConfig] = scene;
    }

    void Simulator::addToSceneRemoveList(SceneConfig *sceneConfig) {
        m_scenes_to_remove.insert(sceneConfig);
    }

    void Simulator::clearSceneRemoveList() {
        for (auto sceneConfig: m_scenes_to_remove) {
            auto scene = m_scene_map[sceneConfig];
            m_scene_map.erase(sceneConfig);
            removeScene(scene);
        }
        m_scenes_to_remove.clear();
    }

    void Simulator::removeScene(Scene *scene) {
        scene->destroy();
        delete scene;
        scene = nullptr;
    }

    void Simulator::refresh() {
        for (auto scene: m_scene_map)
            scene.second->refresh();
    }

    void Simulator::runTimeRange(double t) {
        for (auto scene: m_scene_map)
            scene.second->runTimeRange(t);
    }

    void Simulator::runSingleStep() {
        for (auto scene: m_scene_map)
            scene.second->runSingleStep();
    }

    void Simulator::terminate() {
        for (auto scene: m_scene_map) {
            scene.second->destroy();
            delete scene.second;
        }
        m_scene_map.clear();
        m_scenes_to_remove.clear();
    }

    std::unordered_map<SceneConfig *, Scene *> &Simulator::getSceneMap() {
        return m_scene_map;
    }
}