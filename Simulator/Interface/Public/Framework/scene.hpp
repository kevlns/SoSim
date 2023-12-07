//@author        : Long Shen
//@date          : 2023/10/11
//@description   : 
//@version       : 1.0

#ifndef SOSIM_SCENE_HPP
#define SOSIM_SCENE_HPP

#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>

#include "Public/Framework/solver.hpp"
#include "Public/Framework/framework_config.hpp"
#include "Public/PhysicalSolvers/solver_ref.hpp"

namespace SoSim {

    /**
     * @brief TODO
     */
    class Scene {
    public:
        Scene(CudaConfig *cudaConfig);

        ~Scene();

        void addObjectDefault();

        void removeObject(Object *obj);

        void addToObjectRemoveList(ObjectConfig *objectConfig);

        void clearObjectRemoveList();

        void removeAllObjects();

        void refresh();

        void destroy();

        void runTimeRange(double t);

        void runSingleStep();

        SceneConfig *getConfig();

        std::unordered_map<ObjectConfig *, Object *> &getObjectMap();

    private:

        bool is_start{true};

        SceneConfig *m_scene_config{nullptr};

        std::set<ObjectConfig *> m_obj_to_remove;

        std::unordered_map<ObjectConfig *, Object *> m_obj_map;
    };
}

#endif //SOSIM_SCENE_HPP
