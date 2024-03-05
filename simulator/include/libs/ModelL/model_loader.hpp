//
// Created by ADMIN on 2024/3/5.
//

#ifndef SOSIM_MODEL_LOADER_HPP
#define SOSIM_MODEL_LOADER_HPP

#include <string>
#include <vector>

#include "core/data_type.hpp"

namespace SoSim {

    class ModelLoader {
    public:
        static std::vector<Vec3f> createParticleCube(float particle_radius, Vec3f lb, Vec3f size);

        static std::vector<Vec3f> createParticleBox(float particle_radius, Vec3f lb, Vec3f size, float layer);

        static std::vector<Vec3f> createParticlePlane(float particle_radius, Vec3f lb, Vec3f size, float layer);

        static std::vector<Vec3f> loadParticle3DModel(std::string model_file);
    };

}

#endif //SOSIM_MODEL_LOADER_HPP
