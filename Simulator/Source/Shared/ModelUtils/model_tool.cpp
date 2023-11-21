//@author        : Long Shen
//@date          : 2023/10/27
//@description   :
//@version       : 1.0

#include "Public/Shared/ModelUtils/model_tool.hpp"
#include "Public/Shared/Math/helper_math.hpp"

namespace SoSim {

    extern uint32_t loadObjFromFile(const std::string &file_path, std::vector<float3> &pos) { /*TODO*/return 0; }

    extern void genObjFromJson(const std::string &json_path, std::vector<float3> &pos, std::vector<float3> &vel,
                               std::vector<float> &den, std::vector<Material> &mat, std::vector<Phase> &phase,
                               float &radius) {

        using json = nlohmann::json;

        std::ifstream f(json_path);
        json config = json::parse(f);

        radius = config["unified_particle_radius"];

        for (auto obj: config["objs"]) {

            uint32_t newPartNum = 0;
            std::vector<float3> pos_;
            std::vector<float3> vel_;
            std::vector<float> den_;
            std::vector<Material> mat_;
            std::vector<Phase> phase_;

            if (!obj["source_file"].get<std::string>().empty()) { // incomplete
                newPartNum = loadObjFromFile(obj["source_file"], pos_);

            } else {

                if (obj["default"]["shape"] == "box") {
                    float3 lb = {obj["default"]["lb"][0], obj["default"]["lb"][1], obj["default"]["lb"][2]};
                    float3 size = {obj["default"]["size"][0], obj["default"]["size"][1], obj["default"]["size"][2]};
                    pos_ = generate_box(lb, size, radius);
                    newPartNum = pos_.size();
                } else if (obj["default"]["shape"] == "cube") {
                    float3 lb = {obj["default"]["lb"][0], obj["default"]["lb"][1], obj["default"]["lb"][2]};
                    float3 size = {obj["default"]["size"][0], obj["default"]["size"][1], obj["default"]["size"][2]};
                    pos_ = generate_cube(lb, size, radius);
                    newPartNum = pos_.size();
                } else if (obj["default"]["shape"] == "plane-x") {
                    float3 lb = {obj["default"]["lb"][0], obj["default"]["lb"][1], obj["default"]["lb"][2]};
                    float3 size = {obj["default"]["size"][0], obj["default"]["size"][1], obj["default"]["size"][2]};
                    pos_ = generate_plane_X(lb, size, radius);
                    newPartNum = pos_.size();
                } else if (obj["default"]["shape"] == "cylinder") {
                    // TODO
                    newPartNum = 0;

                }

            }

            float3 v_ = {obj["velStart"][0], obj["velStart"][1], obj["velStart"][2]};
            float d_ = obj["density"];
            vel_ = std::vector<float3>(pos_.size(), v_);
            den_ = std::vector<float>(pos_.size(), d_);

            pos.insert(pos.end(), pos_.begin(), pos_.end());
            vel.insert(vel.end(), vel_.begin(), vel_.end());
            den.insert(den.end(), den_.begin(), den_.end());

            if (obj["mat"].get<std::string>() == "fluid")
                mat_ = std::vector<Material>(newPartNum, Material::FLUID);
            else if (obj["mat"].get<std::string>() == "rigid")
                mat_ = std::vector<Material>(newPartNum, Material::RIGID);
            else if (obj["mat"].get<std::string>() == "elastic")
                mat_ = std::vector<Material>(newPartNum, Material::ELASTIC);
            else if (obj["mat"].get<std::string>() == "bound")
                mat_ = std::vector<Material>(newPartNum, Material::BOUND);
            mat.insert(mat.end(), mat_.begin(), mat_.end());

            if (obj["phase"] == 1)
                phase_ = std::vector<Phase>(newPartNum, Phase::PHASE1);
            else if (obj["phase"] == 2)
                phase_ = std::vector<Phase>(newPartNum, Phase::PHASE2);
            phase.insert(phase.end(), phase_.begin(), phase_.end());
        }

        f.close();
    }

    extern std::vector<float3>
    gen_pos(const nlohmann::json &config, float r) {
        std::vector<float3> pos;
        if (config["shape"] == "cube") {
            float3 lb = make_float3(config["lb"].get<std::vector<float>>());
            float3 size = make_float3(config["size"].get<std::vector<float>>());
            pos = generate_cube(lb, size, r);
        } else if (config["shape"] == "box") {
            float3 lb = make_float3(config["lb"].get<std::vector<float>>());
            float3 size = make_float3(config["size"].get<std::vector<float>>());
            pos = generate_box(lb, size, r);
        } else if (config["shape"] == "plane-x") {
            float3 lb = make_float3(config["lb"].get<std::vector<float>>());
            float3 size = make_float3(config["size"].get<std::vector<float>>());
            pos = generate_plane_X(lb, size, r);
        } else if (config["shape"] == "plane-z") {
            float3 lb = make_float3(config["lb"].get<std::vector<float>>());
            float3 size = make_float3(config["size"].get<std::vector<float>>());
            pos = generate_plane_Z(lb, size, r);
        } else if (config["shape"] == "cylinder") {
            float3 top_center = make_float3(config["top_center"].get<std::vector<float>>());
            float height = config["height"];
            float area_radius = config["area_radius"];
            pos = generate_cylinder(top_center, height, area_radius, r);
        }

        return pos;
    }

    extern std::vector<float3>
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

    extern void write_ply(const std::string &filename, const std::vector<float3> &points) {
        std::ofstream ofs(filename);

        ofs << "ply\n";
        ofs << "format ascii 1.0\n";
        ofs << "element vertex " << points.size() << "\n";
        ofs << "property float x\n";
        ofs << "property float y\n";
        ofs << "property float z\n";
        ofs << "end_header\n";

        for (const auto &point: points) {
            ofs << point.x << " " << point.y << " " << point.z << "\n";
        }

        ofs.close();
    }

    extern void write_ply(const std::string &filename, const float3 *d_pos, uint32_t begin, uint32_t end) {
        auto num = end - begin;
        std::vector<float3> pos(num);
        cudaMemcpy(pos.data(), d_pos, num * sizeof(float3), cudaMemcpyDeviceToHost);

        write_ply(filename, pos);
    }

}
