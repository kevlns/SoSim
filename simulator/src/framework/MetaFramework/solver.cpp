//
// Created by ADMIN on 2024/3/5.
//

#include <iostream>

// solvers
#include "framework/MetaFramework/solver.hpp"


namespace SoSim {

    std::shared_ptr<SolverConfig> Solver::getConfig() {
        return nullptr;
    }

    void Solver::attachObject(std::shared_ptr<Object> object) {
        return;
    }

    void Solver::detachObject(std::shared_ptr<Object> object) {
        return;
    }

    void Solver::detachParticleEmitter(std::shared_ptr<ParticleEmitter> emitter) {
        return;
    }


    bool Solver::initialize() {
        return false;
    }

    void Solver::run(float total_time) {
        return;
    }

    void Solver::step() {
        return;
    }

    void Solver::attachParticleEmitter(std::shared_ptr<ParticleEmitter> emitter) {
        return;
    }
}
