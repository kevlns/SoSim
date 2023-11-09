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

    void Solver::setSolverConfig(const SoSim::SolverConfig *config) {}

    void Solver::attachObject(Object *obj) {}

    void Solver::addParticles(const std::string &obj_json_path) {}

    void Solver::step() {}

}