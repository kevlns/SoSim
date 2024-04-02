//
// Created by ADMIN on 2024/3/5.
//

#ifndef SOSIM_OBJECT_HPP
#define SOSIM_OBJECT_HPP

#include <string>
#include <vector>
#include <optional>

#include "core/math/matrix.hpp"
#include "sim_material.hpp"

namespace SoSim {

    struct ParticleObjectConfig {
        // common
        float particle_radius;
        std::optional<Material> particle_mat;
        Vec3f vel_start{0, 0, 0};

        /* default building */
        std::optional<std::string> shape;
        // cube/box/plane
        Vec3f lb{-0.5, -0.5, -0.5};
        Vec3f size{1, 1, 1};
        // box/plane
        float layer{2};
        // cylinder/sphere
        Vec3f center{0, 0, 0};
        // cylinder
        float bottom_area_radius{1};
        float height{1};
        // sphere
        float volume_radius{1};

        /* load file */
        std::optional<std::string> model_file;

        /* transform */
        Vec3f transfer{0, 0, 0};
        Vec3f scale{1, 1, 1};
        Vec3f rotate{0, 0, 0};

        /* solver related */
        // multi-phase fluid solver
        std::vector<float> phases;          // size is phase num, value is phase volume fraction
    };

    struct MeshObjectConfig {
        // TODO
    };

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

    private:
        void createParticleObject();

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
