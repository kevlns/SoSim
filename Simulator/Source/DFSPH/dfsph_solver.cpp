//@author        : Long Shen
//@date          : 2023/10/11
//@description   : 
//@version       : 1.0

#include "Public/DFSPH/dfsph_solver.hpp"

namespace SoSim::DFSPH {

    DFSPHSolver::DFSPHSolver() {
        m_neighborSearcher = new NeighborSearcher;
    }

    DFSPHSolver::~DFSPHSolver() {
        std::cout << "~DFSPHSolver()\n";
    }

    void DFSPHSolver::initialize() {
        m_neighborSearcher->initialize(m_host_cp.sceneLB, m_host_cp.sceneSize, m_host_cp.totalParticleNum,
                                       m_host_cp.sph_h);
    }

    bool DFSPHSolver::isInitialized() const {
        return m_isInit;
    }

    void DFSPHSolver::run() {

    }

    void DFSPHSolver::destroy() {
        if (!m_isInit)
            return;

        m_neighborSearcher->destroy();
    }

    void DFSPHSolver::setConfig(const Solver::SolverConfig &config) {

    }

    void DFSPHSolver::setParticleRadius(float particle_radius) {
        m_particleRadius = particle_radius;
        m_host_cp.sph_h = 4 * m_particleRadius;
    }

    void DFSPHSolver::setSceneInfo(float3 scene_lb, float3 scene_size) {
        m_host_cp.sceneLB = scene_lb;
        m_host_cp.sceneSize = scene_size;
    }

    void DFSPHSolver::step() {

    }

}