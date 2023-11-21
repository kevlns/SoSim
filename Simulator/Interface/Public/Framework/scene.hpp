//@author        : Long Shen
//@date          : 2023/10/11
//@description   : 
//@version       : 1.0

#ifndef SOSIM_SCENE_HPP
#define SOSIM_SCENE_HPP

#include <iostream>
#include <vector>
#include <set>

#include "Public/Framework/solver.hpp"
#include "Public/Framework/framework_config.hpp"
#include "Public/PhysicalSolvers/solver_header_list.hpp"

namespace SoSim {

    /**
     * @brief TODO
     */
    class Scene {

    public:
        explicit Scene(const SceneConfig &sceneConfig);

        ~Scene() = default;

        Solver *createSolver(SolverInstanceType solverInstance, SolverConfig *solverConfig);

        void addObject(Object* obj);

        void removeObject(Object *obj);

        void addSolver(Solver *solver);

        void removeSolver(Solver *solver);

        void destroy();

        void run();

    private:
        uint32_t m_id{};

        std::set<Solver *> m_solvers;

        std::set<Object *> m_objs;

        SceneConfig *m_config{nullptr};
    };

}

#endif //SOSIM_SCENE_HPP
