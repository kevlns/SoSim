//@author        : Long Shen
//@date          : 2023/11/27
//@description   :
//@version       : 1.0

#include "Public/PhysicalSolvers/SampleSolver/sample_solver.hpp"

namespace SoSim {

    SampleSolver::SampleSolver(SceneConfig *sceneConfig, SolverConfig *solverConfig) {
        m_solverConfig = new SolverConfig;
        memcpy(m_solverConfig, solverConfig, sizeof(SolverConfig));
        m_solverConfig->max_block_num_per_processor = sceneConfig->max_block_num_per_processor;
        m_solverConfig->max_thread_num_per_block = sceneConfig->max_thread_num_per_block;
        m_solverConfig->max_multiprocessor = sceneConfig->max_multiprocessor;
        m_solverConfig->scene_lb = sceneConfig->scene_lb;
        m_solverConfig->scene_size = sceneConfig->scene_size;
    }

    void SampleSolver::refresh() {
        // init block_num and thread_num
    }

    void SampleSolver::runTimeRange(double t) {

    }

    void SampleSolver::runSingleStep() {

    }

    void SampleSolver::destroy() {

    }

    void SampleSolver::attachObject(Object *obj) {

    }

    void SampleSolver::removeObject(Object *obj) {

    }

    SolverConfig *SampleSolver::getConfig() {
        return m_solverConfig;
    }

}
