//@author        : Long Shen
//@date          : 2023/11/22
//@description   :
//@version       : 1.0

#ifndef SOSIM_SAMPLE_SOLVER_HPP
#define SOSIM_SAMPLE_SOLVER_HPP

#include <vector>
#include <iostream>
#include <set>

#include "Public/Framework/solver.hpp"
#include "Public/Framework/object.hpp"

namespace SoSim {

    namespace SAMPLE {

        class SampleSolver : public Solver {

        public:
            SampleSolver() = default;

            ~SampleSolver() override = default;

            void initialize() override;

            bool isInitialized() const override;

            void run() override;

            void destroy() override;

            void setSolverConfig(SolverConfig *solverConfig, const SceneConfig *sceneConfig) override;

            void attachObject(const Object *obj) override;

            void addParticles(const std::string &obj_json_path) override;

        protected:
            void step() override;

        private:
            std::set<Object *> m_objects;

        };
    }
}
#endif //SOSIM_SAMPLE_SOLVER_HPP
