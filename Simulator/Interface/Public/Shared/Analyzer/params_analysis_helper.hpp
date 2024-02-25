//@author        : Long Shen
//@date          : 2023/10/27
//@description   :
//@version       : 1.0

#ifndef SOSIM_PARAMS_ANALYSIS_HELPER_HPP
#define SOSIM_PARAMS_ANALYSIS_HELPER_HPP

#include "Private/BuildSystem/macro_helper.hpp"

namespace SoSim {

    template<typename ParamPType>
    SOSIM_API void dump_avg(unsigned int num, ParamPType d_params);

    template<typename ParamPType>
    SOSIM_API void dump_max(unsigned int num, ParamPType d_params);

    template<typename ParamPType>
    SOSIM_API void dump_min(unsigned int num, ParamPType d_params);

}

#endif //SOSIM_PARAMS_ANALYSIS_HELPER_HPP
