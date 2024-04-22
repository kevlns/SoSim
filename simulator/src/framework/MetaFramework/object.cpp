//
// Created by ADMIN on 2024/3/5.
//

#include <utility>
#include <memory>

#include "framework/MetaFramework/object.hpp"
#include "libs/ModelL/model_loader.hpp"

namespace SoSim {

    Object::Object(unsigned int id) {
        m_id = id;
        m_name = "object_" + std::to_string(id);
        std::cout << "Create an object.\n";
    }

    Object::~Object() {
        destroy();
    }

    std::vector<Vec3f> &Object::getParticles() {
        return m_host_particles;
    }

    std::shared_ptr<ParticleObjectConfig> Object::getParticleObjectConfig() {
        if (!m_particle_obj_config)
            m_particle_obj_config = std::make_shared<ParticleObjectConfig>();

        return m_particle_obj_config;
    }

    std::shared_ptr<MeshObjectConfig> Object::getMeshObjectConfig() {
        if (!m_mesh_obj_config)
            m_mesh_obj_config = std::make_shared<MeshObjectConfig>();

        return m_mesh_obj_config;
    }

    unsigned Object::getParticleNum() const {
        return m_host_particles.size();
    }

    void Object::update() {
        m_host_particles = ModelHelper::create3DParticleModel(m_particle_obj_config);

        if (m_device_cuda_jit_particles)
            cudaFree(m_device_cuda_jit_particles);
        cudaMalloc((void **) &m_device_cuda_jit_particles, getParticleNum() * sizeof(Vec3f));
        cudaMemcpy(m_device_cuda_jit_particles, m_host_particles.data(), getParticleNum() * sizeof(Vec3f),
                   cudaMemcpyHostToDevice);
    }

    void Object::destroy() {
        // todo if m_particles_renderBuffer not empty, delete it

        std::cout << "Object: " << m_name << " destroyed.\n";
    }

    void Object::setName(std::string new_name) {
        m_name = std::move(new_name);
    }

    unsigned Object::getID() const {
        return m_id;
    }

    std::string Object::getName() const {
        return m_name;
    }

    void Object::exportAsPly(const std::string& dir, const std::string& name) {
        ModelHelper::export3DModelAsPly(m_host_particles,
                                        dir,
                                        name);
    }

}
