//
// Created by ADMIN on 2024/3/30.
//

#ifndef SOSIM_RENDER_ENUM_HPP
#define SOSIM_RENDER_ENUM_HPP

namespace SoSim {

    enum RenderBackend {
        OpenGL,

        // others ...
        Vulkan,
        Metal,
    };

    enum RenderMaterial {
        NoneMat,

        // others ...
        Still,
        Wood
    };

    enum RenderPipeline {
        NonePipeline,

        BasicParticleRenderPipeline,
        BasicMeshRenderPipeline,

        CustomRenderPipeline
    };
}

#endif //SOSIM_RENDER_ENUM_HPP
