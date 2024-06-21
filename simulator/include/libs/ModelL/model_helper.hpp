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
#include "framework/MetaFramework/sim_material.hpp"

namespace SoSim {

    enum ObjectShape : uint8_t {
        Cube,
        Box,
        Plane,
        Cylinder,
        Sphere,
        Surface_Sample,
        Volume_Sample
    };

    struct ParticleObjectConfig {
        // common
        std::optional<float> particle_radius;
        std::optional<Material> particle_mat;
        Vec3f vel_start{0, 0, 0};

        /* default building */
        std::optional<ObjectShape> shape;
        // cube/box/plane
        Vec3f lb{-0.5, -0.5, -0.5};
        Vec3f size{1, 1, 1};
        // box/plane
        float layer{2};  // 1 or 2
        // cylinder/sphere
        Vec3f center{0, 0, 0};
        // cylinder
        float bottom_area_radius{1};
        float height{1};
        // sphere
        float volume_radius{1};

        /* load file */
        std::optional<std::string> model_file;
        // for surface sample
        float ratio{1.9};

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

    class ModelHelper {
    public:
        static std::vector<Vec3f> create3DParticleModel(std::shared_ptr<ParticleObjectConfig> config);

        static  std::vector<Vec3f> create3DParticleModel(const ParticleObjectConfig& config);

        static std::vector<float> transformVec3fSetToFloatSet(std::vector<Vec3f> &vec3f_set);

        static Vec3f loadEmitterAgentNormal(const std::string& agent_file);

        static void
        export3DModelAsPly(const std::vector<Vec3f> &particles, const std::string &dir, const std::string &file_name);

        static void
        export3DModelAsPly(const std::vector<Vec3f> &particles, const std::vector<Vec3f> &colors,
                           const std::string &dir, const std::string &file_name);

        static void
        export3DModelAsPly(const std::vector<Vec3f> &particles, const std::vector<Vec2f> &phases,
                           const std::string &dir, const std::string &file_name);

    private:
        static std::vector<Vec3f> create3DParticleCube(float particle_radius, Vec3f lb, Vec3f size);

        static std::vector<Vec3f> create3DParticleBox(float particle_radius, Vec3f lb, Vec3f size, float layer);

        static std::vector<Vec3f> create3DParticlePlane(float particle_radius, Vec3f lb, Vec3f size, float layer);

        static std::vector<Vec3f>
        create3DParticleCylinder(float particleRadius, Vec3f center, float height, float bottom_area_radius);

        static std::vector<Vec3f> create3DParticleSphere(float particle_radius, Vec3f center, float volume_radius);

        static std::vector<Vec3f> loadPly3DModel(std::string ply_file);

        static std::vector<Vec3f> sample3DSurfaceParticle(float particle_radius, std::string file_path, float ratio);

    };

}  // namespace SoSim

#endif  // SOSIM_MODEL_HELPER_HPP
