//@author        : Long Shen
//@date          : 2023/11/27
//@description   :
//@version       : 1.0

#ifndef SOSIM_SAMPLE_SOLVER_HPP
#define SOSIM_SAMPLE_SOLVER_HPP

#include <set>

//#include "Public/Framework/object.hpp"
#include "Public/Framework/solver.hpp"
//#include "Public/Framework/framework_config.hpp"
//#include "Public/Framework/component.hpp"

namespace SoSim {

    class SampleSolver : public Solver {
    public:
        explicit SampleSolver(SceneConfig *sceneConfig, SolverConfig *solverConfig);

        SampleSolver() = delete;

        ~SampleSolver() override = default;

        void refresh() override;

        void runTimeRange(double t) override;

        void runSingleStep() override;

        void destroy() override;

        void attachObject(Object *obj) override;

        void removeObject(Object *obj) override;

        SolverConfig *getConfig() override;

    private:
        SolverConfig *m_solverConfig{nullptr};

        Object *m_fluid_obj{nullptr};
        Object *m_bound_obj{nullptr};
        Object *m_dyRigid_obj{nullptr};
        Object *m_elastic_obj{nullptr};

        std::set<Object *> m_attached_objects;

    };

}

#endif //SOSIM_SAMPLE_SOLVER_HPP
