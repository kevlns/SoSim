//
// Created by ADMIN on 2024/3/30.
//

#ifndef SOSIM_RENDER_COMPONENT_HPP
#define SOSIM_RENDER_COMPONENT_HPP

#include "render_enum.hpp"
#include "framework/MetaFramework/object.hpp"

namespace SoSim {

    struct ObjectRenderComponent {
        bool visible{false};

        RenderPipeline render_system{NonePipeline};
        RenderMaterial render_material{NoneMat};

        // for BasicParticleRenderPipeline
        float particle_radius{0.01};
        Vec3f *device_particles{nullptr};
        Vec3f *device_particle_colors{nullptr};

        // for graphics api
        // gl ...
        // vulkan ...
    };

}

#endif //SOSIM_RENDER_COMPONENT_HPP
