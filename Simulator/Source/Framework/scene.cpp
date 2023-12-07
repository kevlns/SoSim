//@author        : Long Shen
//@date          : 2023/10/11
//@description   : 
//@version       : 1.0

#include "Public/Framework/scene.hpp"
#include "Public/Framework/object.hpp"
#include "Public/PhysicalSolvers/SampleSolver/sample_solver.hpp"

namespace SoSim {
    Scene::Scene(CudaConfig *cudaConfig) {
        m_scene_config = new SceneConfig;
        m_scene_config->max_block_num_per_processor = cudaConfig->max_block_num_per_processor;
        m_scene_config->max_thread_num_per_block = cudaConfig->max_thread_num_per_block;
        m_scene_config->max_multiprocessor = cudaConfig->max_multiprocessor;

        std::cout << "Scene created.\n";
    }

    Scene::~Scene() {
        delete m_scene_config;
        m_scene_config = nullptr;
    }

    void Scene::addObjectDefault() {
        auto obj = new Object;
        auto objConfig = obj->getConfig();

        m_obj_map[objConfig] = obj;
    }

    void Scene::removeObject(Object *obj) {
        obj->destroy();
        delete obj;
        obj = nullptr;
    }

    void Scene::addToObjectRemoveList(ObjectConfig *objectConfig) {
        m_obj_to_remove.insert(objectConfig);
    }

    void Scene::clearObjectRemoveList() {
        if (m_obj_to_remove.empty())
            return;

        for (auto objConfig: m_obj_to_remove) {
            auto obj = m_obj_map[objConfig];
            m_obj_map.erase(objConfig);
            removeObject(obj);
        }
        m_obj_to_remove.clear();
    }

    void Scene::removeAllObjects() {
        destroy();
    }

    void Scene::refresh() {

        std::cout << "Scene refresh.\n";

        // TODO solver refresh
        // ...

        for (auto obj: m_obj_map)
            obj.second->refresh();
    }

    void Scene::destroy() {

        for (auto obj: m_obj_map) {
            obj.second->destroy();
            delete obj.second;
        }
        m_obj_map.clear();
        m_obj_to_remove.clear();

        std::cout << "Scene destroyed.\n";
    }

    void Scene::runTimeRange(double t) {
        // TODO
    }

    void Scene::runSingleStep() {
//        for (auto solver: m_solver_map)
//            solver.second->runSingleStep();
    }

    SceneConfig *Scene::getConfig() {
        return m_scene_config;
    }

    std::unordered_map<ObjectConfig *, Object *> &Scene::getObjectMap() {
        return m_obj_map;
    }

}
