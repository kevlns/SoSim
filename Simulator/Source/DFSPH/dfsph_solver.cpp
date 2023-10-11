//@author        : Long Shen
//@date          : 2023/10/11
//@description   : 
//@version       : 1.0

#include "Public/DFSPH/dfsph_solver.hpp"

namespace SoSim::DFSPH {

    DFSPHSolver::DFSPHSolver() {
        std::cout << "DFSPHSolver()\n";
    }

    DFSPHSolver::~DFSPHSolver() {
        std::cout << "~DFSPHSolver()\n";
    }

    void DFSPHSolver::initialize() {

    }

    bool DFSPHSolver::isInitialized() const {
        return m_isInit;
    }

    void DFSPHSolver::run() {

    }

    void DFSPHSolver::destroy() {

    }

    void DFSPHSolver::setConfig(const Solver::SolverConfig &config) {

    }

    void DFSPHSolver::step() {

    }
}