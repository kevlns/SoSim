//@author        : Long Shen
//@date          : 2023/11/9
//@description   :
//@version       : 1.0

#ifndef SOSIM_SIMULATOR_HPP
#define SOSIM_SIMULATOR_HPP

#include <set>

#include "Public/Framework/scene.hpp"
#include "Public/Framework/framework_config.hpp"
#include "Public/GUI/gui.hpp"

#ifdef USE_SOSIM_GUI
const bool SOSIM_GUI = true;
#else
const bool SOSIM_GUI = false;
#endif

namespace SoSim {
    class Simulator {
    public:
        static Scene* createScene(const SceneConfig& sceneConfig);

        static void removeScene(Scene* scene);

        static void run(bool sosim_gui = SOSIM_GUI);

        static void terminate();

    private:
        static void runPure();

        static void runGUI();

    private:
        static std::set<Scene *> m_scenes;
    };
}

#endif //SOSIM_SIMULATOR_HPP
