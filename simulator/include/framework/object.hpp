//
// Created by ADMIN on 2024/3/5.
//

#ifndef SOSIM_OBJECT_HPP
#define SOSIM_OBJECT_HPP

#include <string>
#include <vector>

#include "core/data_type.hpp"
#include "core/material.hpp"

namespace SoSim {

    struct ParticleObjectConfig {
        // common
        float particle_radius{0.05};
        Material particle_mat;

        /* default building */
        std::string shape{"None"};
        // cube/box/plane
        Vec3f lb{-0.5, -0.5, -0.5};
        Vec3f size{1, 1, 1};
        // box/plane
        float layer{2};
        // cylinder/sphere
        Vec3f center;
        // cylinder
        float bottom_area_radius;
        float height;
        // sphere
        float volume_radius;

        /* load file */
        std::string model_file{"None"};

        /* transform */
        Vec3f transfer{0,0,0};
        Vec3f scale{1,1,1};
        Vec3f rotate{0,0,0};
    };

    class Object {
    public:
        ~Object();

        void createParticleObject(ParticleObjectConfig config);

        std::vector<Vec3f> &getParticles();

        ParticleObjectConfig &getParticleObjectConfig();

    private:
        ParticleObjectConfig m_particleObjectConfig;
        std::vector<Vec3f> m_particles;
    };

}

#endif //SOSIM_OBJECT_HPP
