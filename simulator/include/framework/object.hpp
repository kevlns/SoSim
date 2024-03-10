//
// Created by ADMIN on 2024/3/5.
//

#ifndef SOSIM_OBJECT_HPP
#define SOSIM_OBJECT_HPP

#include <string>
#include <vector>
#include <optional>

#include "core/data_type.hpp"
#include "core/material.hpp"

namespace SoSim {

    struct ParticleObjectConfig {
        // common
        float particle_radius;
        std::optional<Material> particle_mat;
        Vec3f vel_start{0,0,0};

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
    };

    class Object {
    public:
        Object(unsigned id);

        void setConfig(ParticleObjectConfig* config);

        std::vector<Vec3f> &getParticles();

        ParticleObjectConfig* getParticleObjectConfig();

        unsigned getParticleNum() const;

        void update();

        void destroy();

        void rename(std::string new_name);

        std::string getName() const;

        unsigned getID() const;

    private:
        void createParticleObject();

    private:
        std::string m_name;
        unsigned m_id;
        std::optional<ParticleObjectConfig*> m_particleObjectConfig;
        std::vector<Vec3f> m_particles;
        Vec3f *m_particles_renderBuffer{nullptr};
    };

}

#endif //SOSIM_OBJECT_HPP
