//
// Created by ADMIN on 2024/3/5.
//

#ifndef SOSIM_SIM_MATERIAL_HPP
#define SOSIM_SIM_MATERIAL_HPP

namespace SoSim {
    enum Material {
        // common
        COMMON_NEWTON,
        FIXED_BOUND,
        DYNAMIC_RIGID,
        Emitter_Particle,

        // JL21-CT
        JL21CT_NONNEWTON,
        ROTATE_RIGID,
        BOUND_BOWL,
        MOVE_RIGID,
        FLUID_PREPARE_1,

        // IMS-CT
        IMSCT_NONNEWTON,
        STIR_FAN,
        MOVING_BOWL,
        MOVING_TUBE,
        MOVING_COVER,

    };

}

#endif //SOSIM_SIM_MATERIAL_HPP
