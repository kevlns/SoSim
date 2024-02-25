//@author        : Long Shen
//@date          : 2023/12/8
//@description   :
//@version       : 1.0

#ifndef SOSIM_MACRO_HELPER_HPP
#define SOSIM_MACRO_HELPER_HPP

namespace SoSim {

#ifdef SOSIM_EXPORTS
#define SOSIM_API __declspec(dllexport)
#else
#define SOSIM_API __declspec(dllimport)
#endif

}

#endif //SOSIM_MACRO_HELPER_HPP
