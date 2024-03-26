//
// Created by ADMIN on 2024/3/13.
//

#ifndef SOSIM_STATISTIC_UTIL_HPP
#define SOSIM_STATISTIC_UTIL_HPP

#include "core/math/matrix.hpp"

namespace SoSim {

    template<typename T>
    T dump_mean(T *d_ptr, unsigned raw_size, unsigned target_size, unsigned target_start = 0);

    template<typename T>
    T dump_max(T *d_ptr, unsigned raw_size, unsigned target_size, unsigned target_start = 0);

    template<typename T>
    T cal_mean(T *d_ptr, unsigned raw_size, unsigned target_size, unsigned target_start = 0);
}

#endif //SOSIM_STATISTIC_UTIL_HPP
