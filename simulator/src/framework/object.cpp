//
// Created by ADMIN on 2024/3/5.
//

#include <utility>

#include "framework/object.hpp"
#include "libs/ModelL/model_loader.hpp"

namespace SoSim {


    Object::~Object() {

    }

    void Object::createParticleObject(ParticleObjectConfig config) {
        m_particleObjectConfig = std::move(config);
        if (m_particleObjectConfig.shape != "None") {
            if (m_particleObjectConfig.shape == "cube")
                m_particles = ModelLoader::createParticleCube(m_particleObjectConfig.particle_radius,
                                                              m_particleObjectConfig.lb,
                                                              m_particleObjectConfig.size);
            if (m_particleObjectConfig.shape == "box")
                m_particles = ModelLoader::createParticleBox(m_particleObjectConfig.particle_radius,
                                                             m_particleObjectConfig.lb,
                                                             m_particleObjectConfig.size,
                                                             m_particleObjectConfig.layer);
            if (m_particleObjectConfig.shape == "plane")
                m_particles = ModelLoader::createParticlePlane(m_particleObjectConfig.particle_radius,
                                                               m_particleObjectConfig.lb,
                                                               m_particleObjectConfig.size,
                                                               m_particleObjectConfig.layer);
        } else {
            m_particles = ModelLoader::loadParticle3DModel(m_particleObjectConfig.model_file);
        }

        // todo transform

    }

    std::vector<Vec3f> &Object::getParticles() {
        return m_particles;
    }

    ParticleObjectConfig &Object::getParticleObjectConfig() {
        return m_particleObjectConfig;
    }
}
