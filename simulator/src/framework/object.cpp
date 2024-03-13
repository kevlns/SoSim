//
// Created by ADMIN on 2024/3/5.
//

#include <utility>
#include <memory>

#include "framework/object.hpp"
#include "libs/ModelL/model_loader.hpp"

namespace SoSim {

    Object::Object(unsigned int id) {
        m_id = id;
        m_particleObjectConfig = std::make_shared<ParticleObjectConfig>();
        m_particleObjectConfig->name = "obj_" + std::to_string(m_id);

        std::cout << "Create an object.\n";
    }

    Object::~Object() {
        destroy();
    }

    void Object::createParticleObject() {
        if (m_particleObjectConfig->shape.has_value()) {
            if (m_particleObjectConfig->shape == "cube")
                m_particles = ModelLoader::createParticleCube(m_particleObjectConfig->particle_radius,
                                                              m_particleObjectConfig->lb,
                                                              m_particleObjectConfig->size);
            if (m_particleObjectConfig->shape == "box")
                m_particles = ModelLoader::createParticleBox(m_particleObjectConfig->particle_radius,
                                                             m_particleObjectConfig->lb,
                                                             m_particleObjectConfig->size,
                                                             m_particleObjectConfig->layer);
            if (m_particleObjectConfig->shape == "plane")
                m_particles = ModelLoader::createParticlePlane(m_particleObjectConfig->particle_radius,
                                                               m_particleObjectConfig->lb,
                                                               m_particleObjectConfig->size,
                                                               m_particleObjectConfig->layer);
        }

        if (m_particleObjectConfig->model_file.has_value()) {
            m_particles = ModelLoader::loadParticle3DModel(m_particleObjectConfig->model_file.value());
        }

        // todo transform

    }

    std::vector<Vec3f> &Object::getParticles() {
        return m_particles;
    }

    std::shared_ptr<ParticleObjectConfig> Object::getParticleObjectConfig() {
        return m_particleObjectConfig;
    }

    unsigned Object::getParticleNum() const {
        return m_particles.size();
    }

    void Object::update() {
        if (!m_particleObjectConfig->shape.has_value() &&
            !m_particleObjectConfig->model_file.has_value())
            m_particleObjectConfig->shape = "cube";

        createParticleObject();
    }

    void Object::destroy() {
        // todo if m_particles_renderBuffer not empty, delete it

        std::cout << "Object: " << m_particleObjectConfig->name << " destroyed.\n";
    }

    void Object::setName(std::string new_name) {
        m_particleObjectConfig->name = std::move(new_name);
    }

    unsigned Object::getID() const {
        return m_id;
    }

    std::string Object::getName() const {
        return m_particleObjectConfig->name;
    }

}
