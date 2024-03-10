//
// Created by ADMIN on 2024/3/5.
//

#include <utility>

#include "framework/object.hpp"
#include "libs/ModelL/model_loader.hpp"

namespace SoSim {

    Object::Object(unsigned int id) {
        m_id = id;
        m_name = "obj_" + std::to_string(m_id);
    }

    void Object::setConfig(ParticleObjectConfig *config) {
        m_particleObjectConfig = config;
    }

    void Object::createParticleObject() {
        if (m_particleObjectConfig.value()->shape.has_value()) {
            if (m_particleObjectConfig.value()->shape == "cube")
                m_particles = ModelLoader::createParticleCube(m_particleObjectConfig.value()->particle_radius,
                                                              m_particleObjectConfig.value()->lb,
                                                              m_particleObjectConfig.value()->size);
            if (m_particleObjectConfig.value()->shape == "box")
                m_particles = ModelLoader::createParticleBox(m_particleObjectConfig.value()->particle_radius,
                                                             m_particleObjectConfig.value()->lb,
                                                             m_particleObjectConfig.value()->size,
                                                             m_particleObjectConfig.value()->layer);
            if (m_particleObjectConfig.value()->shape == "plane")
                m_particles = ModelLoader::createParticlePlane(m_particleObjectConfig.value()->particle_radius,
                                                               m_particleObjectConfig.value()->lb,
                                                               m_particleObjectConfig.value()->size,
                                                               m_particleObjectConfig.value()->layer);
        }

        if (m_particleObjectConfig.value()->model_file.has_value()) {
            m_particles = ModelLoader::loadParticle3DModel(m_particleObjectConfig.value()->model_file.value());
        }

        // todo transform

    }

    std::vector<Vec3f> &Object::getParticles() {
        return m_particles;
    }

    ParticleObjectConfig *Object::getParticleObjectConfig() {
        return m_particleObjectConfig.value();
    }

    unsigned Object::getParticleNum() const {
        return m_particles.size();
    }

    void Object::update() {
        if (!m_particleObjectConfig.value()->shape.has_value() &&
            !m_particleObjectConfig.value()->model_file.has_value())
            m_particleObjectConfig.value()->shape = "cube";

        createParticleObject();
    }

    void Object::destroy() {
        // todo if m_particles_renderBuffer not empty, delete it

        std::cout << "Object: " << m_name << " destroyed.\n";
    }

    void Object::rename(std::string new_name) {
        m_name = std::move(new_name);
    }

    unsigned Object::getID() const {
        return m_id;
    }

    std::string Object::getName() const {
        return m_name;
    }

}
