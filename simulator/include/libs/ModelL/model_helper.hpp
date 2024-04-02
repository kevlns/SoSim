//
// Created by ADMIN on 2024/3/23.
//

/**
 * description: novel model util, haven't integrate into framework
 */

#ifndef SOSIM_MODEL_HELPER_HPP
#define SOSIM_MODEL_HELPER_HPP

#include <vector_types.h>
#include <iostream>
#include <vector>
#include <string>
#include <optional>

#include "core/math/matrix.hpp"

namespace SoSim {

    enum ObjectShape : uint8_t {
        Cube,
        Box,
        Plane,
        Cylinder,
        Sphere
    };

    struct ParticleModelConfig {
        std::optional<float> particle_radius;

        /* default building */
        std::optional<ObjectShape> shape;
        // cube/box/plane
        Vec3f lb{-1, -1, -1};
        Vec3f size{2, 2, 2};
        // box/plane
        float layer{2};
        // cylinder/sphere
        Vec3f center{0, 0, 0};
        // cylinder
        float bottom_area_radius{1};
        float height{1};
        // sphere
        float volume_radius{1};

        /* load 3D model */
        std::optional<std::string> ply_file;
    };

    class ModelHelper {
    public:
        static std::vector<Vec3f> create3DParticleModel(ParticleModelConfig &config);

        static std::vector<float> transformVec3fSetToFloatSet(std::vector<Vec3f> &vec3f_set);

        static void
        export3DModelAsPly(const std::vector<Vec3f> &particles, const std::string &dir, const std::string &file_name);

        static void
        export3DModelAsPly(const std::vector<Vec3f> &particles, const std::vector<Vec3f> &colors,
                           const std::string &dir, const std::string &file_name);

    private:
        static std::vector<Vec3f> create3DParticleCube(float particle_radius, Vec3f lb, Vec3f size);

        static std::vector<Vec3f> create3DParticleBox(float particle_radius, Vec3f lb, Vec3f size, float layer);

        static std::vector<Vec3f> create3DParticlePlane(float particle_radius, Vec3f lb, Vec3f size, float layer);

        static std::vector<Vec3f>
        create3DParticleCylinder(float particleRadius, Vec3f center, float height, float bottom_area_radius);

        static std::vector<Vec3f> create3DParticleSphere(float particle_radius, Vec3f center, float volume_radius);

        static std::vector<Vec3f> loadPly3DModel(std::string ply_file);
    };

}  // namespace SoSim

#endif  // SOSIM_MODEL_HELPER_HPP
