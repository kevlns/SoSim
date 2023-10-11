//@author        : Long Shen
//@date          : 2023/10/11
//@description   : 
//@version       : 1.0

#ifndef SOSIM_SCENE_HPP
#define SOSIM_SCENE_HPP

#include <iostream>

namespace SoSim {

    /**
     * @brief TODO
     */
    class Scene {

    public:
        Scene() = default;

        ~Scene() = default;

        void addSolver();

        void removeSolver();

        void destroy();

    private:
        uint32_t m_id;
    };

}

#endif //SOSIM_SCENE_HPP
