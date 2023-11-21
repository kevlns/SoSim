//@author        : Long Shen
//@date          : 2023/10/11
//@description   : 
//@version       : 1.0

#include "Public/Framework/scene.hpp"
#include "Public/PhysicalSolvers/solver_header_list.hpp"

namespace SoSim {

    Scene::Scene(const SoSim::SceneConfig &sceneConfig) {
        if (!m_config) {
            m_config = new SceneConfig;
            memcpy(m_config, &sceneConfig, sizeof(SceneConfig));
        }
    }

    Solver *Scene::createSolver(SolverInstanceType solverInstance, SolverConfig *solverConfig) {
        Solver *solver = nullptr;
        switch (solverInstance) {
            case SolverInstanceType::MMSPH_MMSPHSolver:
                solver = new MMSPH::MMSPHSolver;
                break;
            default:
                break;
        }

        if (solver) {
            solver->setSolverConfig(solverConfig, m_config);
            m_solvers.insert(solver);
        }

        return solver;
    }

    void Scene::addObject(Object *obj) {

    }

    void Scene::removeObject(Object *obj) {

    }

    void Scene::addSolver(Solver *solver) {
        if (m_solvers.count(solver) == 0) {
            m_solvers.insert(solver);
        }
    }

    void Scene::removeSolver(Solver *solver) {
        if (m_solvers.count(solver) != 0) {
            solver->destroy();
            m_solvers.erase(solver);
        }
    }

    void Scene::destroy() {
        if (!m_solvers.empty()) {
            for (auto solver: m_solvers)
                solver->destroy();
        }

        if (!m_objs.empty()) {
            for (auto obj: m_objs)
                obj->destroy();
        }

        delete m_config;
        m_config = nullptr;

        std::cout << "Scene destructed.\n";
    }

    void Scene::run() {
        // 串行执行
        for (auto solver: m_solvers)
            solver->run();
    }

}
