//@author        : Long Shen
//@date          : 2023/10/11
//@description   : 
//@version       : 1.0

#ifndef SOSIM_SOLVER_HPP
#define SOSIM_SOLVER_HPP

namespace SoSim {

    class Solver {
    public:
        struct SolverConfig {
            double m_dt{0};
            double m_cur_time{0};
        };

    public:
        Solver() = default;

        virtual ~Solver() = default;

        virtual void initialize() = 0;

        virtual bool isInitialized() const = 0;

        virtual void run() = 0;

        virtual void destroy() = 0;

        virtual void setConfig(const SolverConfig &config) = 0;

    protected:
        virtual void step() = 0;
    };

}

#endif //SOSIM_SOLVER_HPP
