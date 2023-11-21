//@author        : Long Shen
//@date          : 2023/10/11
//@description   : 
//@version       : 1.0

#ifndef SOSIM_SOLVER_HPP
#define SOSIM_SOLVER_HPP

#include <vector_types.h>
#include <iostream>

#include "Public/Framework/framework_config.hpp"
#include "Public/Framework/object.hpp"

namespace SoSim {

    class Solver {

    public:
        Solver() = default;

        virtual ~Solver() = default;

        virtual void initialize() = 0;

        virtual bool isInitialized() const = 0;

        virtual void run() = 0;

        virtual void destroy() = 0;

        virtual void setSolverConfig(SolverConfig *solverConfig, const SceneConfig *sceneConfig) = 0;

        virtual void attachObject(const Object *obj) = 0;

        virtual void addParticles(const std::string &obj_json_path) = 0;

    protected:
        virtual void step() = 0;

    };

}

#endif //SOSIM_SOLVER_HPP
