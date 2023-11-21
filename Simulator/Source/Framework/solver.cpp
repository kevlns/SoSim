//@author        : Long Shen
//@date          : 2023/10/11
//@description   : 
//@version       : 1.0

#include "Public/Framework/solver.hpp"

#include <iostream>

namespace SoSim {

    void Solver::initialize() {}

    bool Solver::isInitialized() const { return false; }

    void Solver::run() {}

    void Solver::destroy() {}

    void Solver::setSolverConfig(SolverConfig *solverConfig, const SceneConfig *sceneConfig) {}

    void Solver::attachObject(const Object *obj) {}

    void Solver::addParticles(const std::string &obj_json_path) {}

    void Solver::step() {}

}