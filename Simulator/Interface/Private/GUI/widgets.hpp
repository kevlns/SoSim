//@author        : Long Shen
//@date          : 2023/11/26
//@description   :
//@version       : 1.0

#ifndef SOSIM_WIDGETS_HPP
#define SOSIM_WIDGETS_HPP

#include "Public/Framework/simulator.hpp"
#include "Public/Framework/scene.hpp"
#include "Public/Framework/framework_config.hpp"

namespace SoSim {

    extern void ShowMainWindow(Simulator *simulator);

    extern void ShowScenesConfig(Simulator *simulator);

    extern void ShowObjectDetail(Scene *scene);
}

#endif //SOSIM_WIDGETS_HPP
