//
// Created by ADMIN on 2024/3/5.
//

#ifndef SOSIM_SIM_MATERIAL_HPP
#define SOSIM_SIM_MATERIAL_HPP

namespace SoSim {
    enum Material {
        // common
        COMMON_FLUID,
        FIXED_BOUND,
        DYNAMIC_RIGID,

        // JL21-CT
        JL21CT_NONNEWTON,
        ROTATE_RIGID,
        BOUND_BOWL,
        MOVE_RIGID,
        FLUID_PREPARE_1,

        // IMS-CT
        IMSCT_NONNEWTON,

    };

}

#endif //SOSIM_SIM_MATERIAL_HPP
