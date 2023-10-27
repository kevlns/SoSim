//@author        : Long Shen
//@date          : 2023/10/27
//@description   :
//@version       : 1.0

#include "Public/Shared/ModelUtils/model_tool.hpp"
#include "Public/Shared/Math/helper_math.hpp"

namespace SoSim {

    extern void genObjFromJson() {/* TODO */};

    std::vector<float3>
    generate_cube(float3 cubeLB, float3 cubeSize, float particleRadius) {
        std::vector<float3> pos;
        auto diameter = 2 * particleRadius;

        for (float z = particleRadius + cubeLB.z; z < cubeLB.z + cubeSize.z; z += diameter) {
            for (float y = particleRadius + cubeLB.y; y < cubeLB.y + cubeSize.y; y += diameter) {
                for (float x = particleRadius + cubeLB.x; x < cubeLB.x + cubeSize.x; x += diameter) {
                    float3 _pos = make_float3(x, y, z);

                    pos.push_back(_pos);
                }
            }
        }

        return pos;
    }

    extern std::vector<float3>
    generate_box(float3 boxLB, float3 boxSize, float particleRadius) {
        std::vector<float3> pos;

        int numParticles[] = {
                static_cast<int>(boxSize.x / (2.0 * particleRadius)),
                static_cast<int>(boxSize.y / (2.0 * particleRadius)),
                static_cast<int>(boxSize.z / (2.0 * particleRadius))
        };

        for (int i = 0; i < numParticles[0]; ++i) {
            for (int j = 0; j < numParticles[1]; ++j) {
                for (int k = 0; k < numParticles[2]; ++k) {
                    // If this particle is in the first two or last two layers in any dimension...
                    if (i < 2 || i >= numParticles[0] - 2 || j < 2 || j >= numParticles[1] - 2 || k < 2 ||
                        k >= numParticles[2] - 2) {
                        float3 p;
                        p.x = boxLB.x + particleRadius + 2.0 * particleRadius * i;
                        p.y = boxLB.y + particleRadius + 2.0 * particleRadius * j;
                        p.z = boxLB.z + particleRadius + 2.0 * particleRadius * k;
                        pos.push_back(p);
                    }
                }
            }
        }

        return pos;
    }

    extern std::vector<float3>
    generate_plane_X(float3 planeLB, float3 planeSize, float particleRadius) {
        std::vector<float3> pos;
        auto diameter = 2 * particleRadius;

        for (float z = particleRadius + planeLB.z; z < planeLB.z + planeSize.z; z += diameter) {
            for (float y = particleRadius + planeLB.y, cnt = 0; cnt < 2; y += diameter, cnt += 1) {
                for (float x = particleRadius + planeLB.x; x < planeLB.x + planeSize.x; x += diameter) {
                    float3 _pos = make_float3(x, y, z);

                    pos.push_back(_pos);
                }
            }
        }

        return pos;
    }

    extern std::vector<float3>
    generate_plane_Z(float3 planeLB, float3 planeSize, float particleRadius) {
        std::vector<float3> pos;
        auto diameter = 2 * particleRadius;

        for (float z = particleRadius + planeLB.z, cnt = 0; cnt < 2; z += diameter, cnt += 1) {
            for (float y = particleRadius + planeLB.y; y < planeLB.y + planeSize.y; y += diameter) {
                for (float x = particleRadius + planeLB.x; x < planeLB.x + planeSize.x; x += diameter) {
                    float3 _pos = make_float3(x, y, z);

                    pos.push_back(_pos);
                }
            }
        }

        return pos;
    }

    extern std::vector<float3>
    generate_cylinder(float3 topCenter, float height, float areaRadius, float particleRadius) {

        std::vector<float3> pos;
        float diameter = 2 * particleRadius;
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

    extern std::pair<std::vector<float3>, float>
    import_ply_model(const std::string &file_name) {/*TODO*/ return {}; }
//    extern std::pair<std::vector<float3>, float>
//    import_ply_model(const std::string &file_name) {
//
//// 打开 PLY 文件
//        std::ifstream ss(file_name, std::ios::binary);
//        if (!ss) {
//            std::cerr << "Error:: opening file: " << file_name << std::endl;
//            return {{}, 0};
//        }
//
//        tinyply::PlyFile file;
//        file.parse_header(ss);
//
//        std::shared_ptr<tinyply::PlyData> vertices;
//        try {
//            vertices = file.request_properties_from_element("vertex", {"x", "y", "z"});
//        }
//        catch (const std::exception &e) {
//            std::cerr << "Error:: tinyply exception: " << e.what() << std::endl;
//            return {{}, 0};
//        }
//
//        file.read(ss);
//
//        std::vector<float3> pos;
//        for (int i = 0; i < vertices->count; ++i) {
//            float3 point;
//            point.x = reinterpret_cast<float *>(vertices->buffer.get())[3 * i + 0];
//            point.y = reinterpret_cast<float *>(vertices->buffer.get())[3 * i + 1];
//            point.z = reinterpret_cast<float *>(vertices->buffer.get())[3 * i + 2];
//            pos.push_back(point);
//        }
//
//        auto radius_min = 100.f;
//        auto radius_max = 0.f;
//        for (int i = 1; i < pos.size(); ++i) {
//            radius_max = std::max(radius_max, length(pos[i] - pos[i - 1]) / 2);
//            radius_min = std::min(radius_min, length(pos[i] - pos[i - 1]) / 2);
//        }
//
//        float radius = 0;
//        return {pos, radius};
//    }

}
