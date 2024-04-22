//
// Created by ADMIN on 2024/3/5.
//

#ifndef SOSIM_OBJECT_HPP
#define SOSIM_OBJECT_HPP

#include <string>
#include <vector>
#include <optional>

#include "sim_material.hpp"
#include "core/math/matrix.hpp"
#include "libs/ModelL/model_helper.hpp"

namespace SoSim {
    class Object {
    public:
        Object(unsigned id);

        ~Object();

        std::vector<Vec3f> &getParticles();

        std::shared_ptr<ParticleObjectConfig> getParticleObjectConfig();

        std::shared_ptr<MeshObjectConfig> getMeshObjectConfig();

        unsigned getParticleNum() const;

        void update();

        void setName(std::string new_name);

        std::string getName() const;

        unsigned getID() const;

        void exportAsPly(const std::string &dir, const std::string &name);

    private:
        void destroy();

    private:
        unsigned m_id;
        std::string m_name;

        std::shared_ptr<ParticleObjectConfig> m_particle_obj_config;
        std::shared_ptr<MeshObjectConfig> m_mesh_obj_config;
//        std::shared_ptr<ObjectRenderComponent> m_render_info;
        std::vector<Vec3f> m_host_particles;

    public:
        Vec3f *m_device_cuda_jit_particles;
    };

}

#endif //SOSIM_OBJECT_HPP
