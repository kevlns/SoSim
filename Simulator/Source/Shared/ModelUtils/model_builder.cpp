//@author        : Long Shen
//@date          : 2023/11/26
//@description   :
//@version       : 1.0

#include "Public/Shared/ModelUtils/model_builder.hpp"
#include "Public/Shared/Math/helper_math.hpp"

namespace SoSim {

    extern void handleTransformer(std::vector<float3> &src_pos, const ObjectConfig *objectConfig) {
        // TODO
    }

    extern std::vector<float3> genParticleCube(const ObjectConfig *object_config) {
        std::vector<float3> pos;

        const auto particleRadius = object_config->particle_radius;
        const auto lb = object_config->lb;
        const auto size = object_config->size;
        const auto diameter = 2 * particleRadius;

        for (float z = particleRadius + lb.z; z < lb.z + size.z;) {
            for (float y = particleRadius + lb.y; y < lb.y + size.y;) {
                for (float x = particleRadius + lb.x; x < lb.x + size.x;) {
                    float3 _pos = make_float3(x, y, z);
                    pos.push_back(_pos);
                    x += diameter;
                }
                y += diameter;
            }
            z += diameter;
        }

        return pos;
    }

    extern std::vector<float3> genParticleBox(const ObjectConfig *object_config) {
        std::vector<float3> pos;

        const auto particleRadius = object_config->particle_radius;
        const auto lb = object_config->lb;
        const auto size = object_config->size;
        const auto diameter = 2 * particleRadius;

        int numParticles[] = {
                static_cast<int>(size.x / (2.0 * particleRadius)),
                static_cast<int>(size.y / (2.0 * particleRadius)),
                static_cast<int>(size.z / (2.0 * particleRadius))
        };

        for (int i = 0; i < numParticles[0]; ++i) {
            for (int j = 0; j < numParticles[1]; ++j) {
                for (int k = 0; k < numParticles[2]; ++k) {
                    // If this particle is in the first two or last two layers in any dimension...
                    if (i < 2 || i >= numParticles[0] - 2 || j < 2 || j >= numParticles[1] - 2 || k < 2 ||
                        k >= numParticles[2] - 2) {
                        float3 p;
                        p.x = lb.x + particleRadius + 2.0 * particleRadius * i;
                        p.y = lb.y + particleRadius + 2.0 * particleRadius * j;
                        p.z = lb.z + particleRadius + 2.0 * particleRadius * k;
                        pos.push_back(p);
                    }
                }
            }
        }

        return pos;
    }

    extern std::vector<float3> genParticlePlaneX(const ObjectConfig *object_config) {
        std::vector<float3> pos;

        const auto particleRadius = object_config->particle_radius;
        const auto lb = object_config->lb;
        const auto size = object_config->size;
        const auto diameter = 2 * particleRadius;

        for (float z = particleRadius + lb.z; z < lb.z + size.z; z += diameter) {
            for (float y = particleRadius + lb.y, cnt = 0; cnt < 2; y += diameter, cnt += 1) {
                for (float x = particleRadius + lb.x; x < lb.x + size.x; x += diameter) {
                    float3 _pos = make_float3(x, y, z);

                    pos.push_back(_pos);
                }
            }
        }

        return pos;
    }

    extern std::vector<float3> genParticlePlaneZ(const ObjectConfig *object_config) {
        std::vector<float3> pos;

        const auto particleRadius = object_config->particle_radius;
        const auto lb = object_config->lb;
        const auto size = object_config->size;
        const auto diameter = 2 * particleRadius;

        for (float z = particleRadius + lb.z, cnt = 0; cnt < 2; z += diameter, cnt += 1) {
            for (float y = particleRadius + lb.y; y < lb.y + size.y; y += diameter) {
                for (float x = particleRadius + lb.x; x < lb.x + size.x; x += diameter) {
                    float3 _pos = make_float3(x, y, z);

                    pos.push_back(_pos);
                }
            }
        }

        return pos;
    }

    extern std::vector<float3> genParticleCylinder(const ObjectConfig *object_config) {
        std::vector<float3> pos;

        const auto particleRadius = object_config->particle_radius;
        const auto diameter = 2 * particleRadius;
        const auto topCenter = object_config->top_center;
        const auto height = object_config->height;
        const auto areaRadius = object_config->area_radius;

        float y0 = topCenter.y;

        for (float y = y0 - particleRadius; y >= y0 - height; y -= diameter) {
            float x0 = topCenter.x - areaRadius;

            for (float x = x0 + particleRadius; x <= topCenter.x + areaRadius; x += diameter) {
                float m_cos = fabs(topCenter.x - x) / areaRadius;
                float length = areaRadius * sqrt(1 - m_cos * m_cos);
                float z0 = topCenter.z - length;
                for (float z = z0 + particleRadius; z <= topCenter.z + length; z += diameter) {
                    float3 _pos = {x, y, z};
                    pos.push_back(_pos);
                }
            }
        }

        return pos;
    }

    extern std::vector<float3> genFromObjectConfig(const ObjectConfig *objectConfig) {
        if (objectConfig->model_file_path != "None") {
            // TODO load model
            // ...
            return {};
        } else {
            if (objectConfig->shape == "cube")
                return genParticleCube(objectConfig);
            else if (objectConfig->shape == "box")
                return genParticleBox(objectConfig);
            else if (objectConfig->shape == "cylinder")
                return genParticleCylinder(objectConfig);
            else if (objectConfig->shape == "plane-x")
                return genParticlePlaneX(objectConfig);
            else if (objectConfig->shape == "plane-z")
                return genParticlePlaneZ(objectConfig);

            return {};
        }
    }
}
