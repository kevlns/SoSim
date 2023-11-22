//@author        : Long Shen
//@date          : 2023/11/9
//@description   :
//@version       : 1.0

#ifdef USE_SOSIM_GUI

#endif

#include "Public/Framework/simulator.hpp"

namespace SoSim {

    Scene *Simulator::createScene(const SceneConfig &sceneConfig) {
        auto *scene = new Scene(sceneConfig);
        m_scenes.insert(scene);
        return scene;
    }

    void Simulator::removeScene(Scene *scene) {
        if (m_scenes.count(scene) != 0) {
            scene->destroy();
            m_scenes.erase(scene);
        }
    }

    void Simulator::run() {
        // 串行执行
        for (auto scene: m_scenes)
            scene->run();
    }

    void Simulator::terminate() {
        for (auto scene: m_scenes)
            scene->destroy();
    }

    std::set<Scene *> Simulator::m_scenes;
}