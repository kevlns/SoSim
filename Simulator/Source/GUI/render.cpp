//@author        : Long Shen
//@date          : 2023/11/27
//@description   :
//@version       : 1.0

#include "Public/GUI/render.hpp"
#include "glad/glad.h"

namespace SoSim {


    Renderer::Renderer() {

        const unsigned int SCR_WIDTH = 1080;

        const unsigned int SCR_HEIGHT = 720;

        m_camera = new Camera(glm::vec3(0.0f, -12.f, 20.0f));

        m_model = glm::mat4(1.0f);
        m_projection = glm::perspective(glm::radians(m_camera->Zoom), (float) SCR_WIDTH / (float) SCR_HEIGHT, 0.1f,
                                        100.0f);

        // prepare shaders
        Shader *drawSphereShader = new Shader("Shaders/draw_sphere.vert", "Shaders/draw_sphere.frag");
        m_shaders["drawSphereShader"] = drawSphereShader;
    }

    void Renderer::drawParticles_tmp(Component *component, float3 color) {
        if (m_oglResource.count(component) == 0) {
            OglResource *oglResource = new OglResource;
            
            glGenVertexArrays(1, &oglResource->VAO);
            glBindVertexArray(oglResource->VAO);

            unsigned int idx[] = { 0, 1, 2, 2, 3, 0 };
            glGenBuffers(1, &oglResource->EBO);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, oglResource->EBO);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(idx), idx, GL_STATIC_DRAW);

            float vertices[] = {
                    -0.5f, -0.5f, 0.0f,
                    0.5f, -0.5f, 0.0f,
                    0.5f, 0.5f, 0.0f,
                    -0.5f, 0.5f, 0.0f
            };
            glGenBuffers(1, &oglResource->VBO);
            glBindBuffer(GL_ARRAY_BUFFER, oglResource->VBO);
            glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
            glEnableVertexAttribArray(0);
            glBindBuffer(GL_ARRAY_BUFFER, 0);

            float uvs[] = { 0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 1.0f };
            glGenBuffers(1, &oglResource->VUV);
            glBindBuffer(GL_ARRAY_BUFFER, oglResource->VUV);
            glBufferData(GL_ARRAY_BUFFER, sizeof(uvs), uvs, GL_STATIC_DRAW);

            glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
            glEnableVertexAttribArray(1);
            glBindBuffer(GL_ARRAY_BUFFER, 0);

            m_oglResource[component] = oglResource;
        }




//        auto oglResource = m_oglResource[component];
//        auto shader = m_shaders["drawSphereShader"];
//        shader->use();
//        shader->setMat4("model", m_model);
//        shader->setMat4("projection", m_projection);
//
//        glm::mat4 view = m_camera->GetViewMatrix();
//        shader->setMat4("view", view);
//
//        glBindBuffer(GL_ARRAY_BUFFER, oglResource->oglResource->VBO);
//        glDrawArrays(GL_POINTS, 0, );
//        glBindBuffer(GL_ARRAY_BUFFER, 0);

    }

    void Renderer::destroy() {
        for (auto shader: m_shaders) {
            shader.second->destroy();
            delete shader.second;
        }

        m_shaders.clear();
    }

}