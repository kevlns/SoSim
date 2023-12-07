//@author        : Long Shen
//@date          : 2023/11/27
//@description   :
//@version       : 1.0

#ifndef SOSIM_RENDER_HPP
#define SOSIM_RENDER_HPP

#include <vector_types.h>
#include <unordered_map>

#include "Public/GUI/shader.hpp"
#include "Public/GUI/camera.hpp"
#include "Private/GUI/opengl_resource.hpp"
#include "Public/Framework/component.hpp"

namespace SoSim {

    class Renderer {
    public:
        Renderer();

        void destroy();

        void drawParticles_tmp(Component *component, float3 color);

    private:
        Camera *m_camera{nullptr};

        std::unordered_map<std::string, Shader *> m_shaders;

        std::unordered_map<Component *, OglResource *> m_oglResource;

    private:
        glm::mat4 m_model;
        glm::mat4 m_projection;

    };

}

#endif //SOSIM_RENDER_HPP
