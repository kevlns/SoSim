//@author        : Long Shen
//@date          : 2023/10/30
//@description   :
//@version       : 1.0

#ifndef SOSIM_ENUM_HPP
#define SOSIM_ENUM_HPP

namespace SoSim {
    enum MaterialType : uint8_t {
        MAT_NONE = 0,

        FLUID,
        RIGID,
        ELASTIC,
        BOUND,
        /**
         *  add your own Material blow
         */

        MAT_BOTTOM, // used for compute num of case
    };
    static inline const char *matItems[] = {"None",
                                            "fluid",
                                            "rigid",
                                            "elastic",
                                            "bound"};

    enum PhaseType : uint8_t {
        PHASE_NONE = 0,

        PHASE_1,
        PHASE_2,
        /**
         *  add your own PHASE blow
         */

        PHASE_BOTTOM, // used for compute num of case
    };
    static inline const char *phaseItems[] = {"None",
                                              "phase-1",
                                              "phase-2"};

    enum ComponentType : uint8_t {
        COMPONENT_NONE = 0,

        BASE_MOVE_COMPONENT,
        BASE_PARTICLE_COMPONENT,
        PARTICLE_GL_RRENDERABLE_COMPONENT,
        /**
         *  add your own component blow
         */

        COMPONENT_BOTTOM, // used for compute num of case
    };
    static inline const char *componentItems[] = {"None",
                                                  "Attach BaseMoveComponent",
                                                  "Attach BaseParticleComponent",
                                                  "Attach Particle_GLRenderableComponent"};
}

#endif //SOSIM_ENUM_HPP
