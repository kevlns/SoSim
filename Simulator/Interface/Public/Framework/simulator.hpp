//@author        : Long Shen
//@date          : 2023/11/9
//@description   :
//@version       : 1.0

#ifndef SOSIM_SIMULATOR_HPP
#define SOSIM_SIMULATOR_HPP

#include <set>

#include "Public/Framework/scene.hpp"
#include "Public/Framework/framework_config.hpp"

#ifdef USE_SOSIM_GUI
const bool SOSIM_GUI = true;
#else
const bool SOSIM_GUI = false;
#endif

namespace SoSim {
    class Simulator {
    public:
        Simulator();

        ~Simulator();

        void addSceneDefault();

        void addToSceneRemoveList(SceneConfig *sceneConfig);

        void clearSceneRemoveList();

        void removeScene(Scene* scene);

        void refresh();

        void runTimeRange(double t);

        void runSingleStep();

        void terminate();

        std::unordered_map<SceneConfig *, Scene *>& getSceneMap();

    private:
        CudaConfig *m_cudaConfig{nullptr};

        std::set<SceneConfig *> m_scenes_to_remove;

        std::unordered_map<SceneConfig *, Scene *> m_scene_map;

    };
}

#endif //SOSIM_SIMULATOR_HPP
