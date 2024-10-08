//
// Created by ADMIN on 2024/3/23.
//

#include <fstream>
#include <filesystem>

#include "libs/ModelL/model_helper.hpp"

#include "assimp/Importer.hpp"
#include "assimp/postprocess.h"
#include "assimp/scene.h"

namespace SoSim {

    std::vector<Vec3f> ModelHelper::create3DParticleModel(std::shared_ptr<ParticleObjectConfig> config) {
        static std::vector<Vec3f> particles;

        if (!config->particle_radius.has_value()) {
            config->particle_radius = 0.1f;
            std::cout << "Warning: Not specified particle radius, use default 0.1.\n";
        }

        if (!config->shape.has_value() && !config->model_file.has_value()) {
            config->shape = ObjectShape::Cube;
            std::cout << "Warning: Not specified object shape and ply model file, use default Cube shape.\n";
        }

        if (config->shape.has_value()) {

            switch (config->shape.value()) {
                case Cube:
                    particles = create3DParticleCube(config->particle_radius.value(),
                                                     config->lb,
                                                     config->size);
                    break;
                case Box:
                    particles = create3DParticleBox(config->particle_radius.value(),
                                                    config->lb,
                                                    config->size,
                                                    config->layer);
                    break;
                case Plane:
                    particles = create3DParticlePlane(config->particle_radius.value(),
                                                      config->lb,
                                                      config->size,
                                                      config->layer);
                    break;
                case Cylinder:
                    particles = create3DParticleCylinder(config->particle_radius.value(),
                                                         config->center,
                                                         config->height,
                                                         config->bottom_area_radius);
                    break;
                case Sphere:
                    particles = create3DParticleSphere(config->particle_radius.value(),
                                                       config->center,
                                                       config->volume_radius);
                    break;
                default:
                    std::cout << "Ops! No matching shape.\n";
                    break;
            }
        }

        if (config->model_file.has_value()) {
            std::filesystem::path fp(config->model_file.value());
            auto postfix = fp.extension().string();
            if (postfix == ".ply")
                particles = ModelHelper::loadPly3DModel(config->model_file.value());
        }

        for(auto &particle: particles){
            particle += config->transfer;
        }

        return particles;
    }

    void ModelHelper::export3DModelAsPly(const std::vector<Vec3f> &particles, const std::string &dir,
                                         const std::string &file_name) {
        auto dir_ = dir;
#if defined(WIN32)
        size_t pos = 0;
        while ((pos = dir_.find('/', pos)) != std::string::npos) {
            dir_.replace(pos, 1, "\\");
            pos += 1;
        }

        // TODO other platform file path transformation
        // ...
#endif

        if (!std::filesystem::exists(dir_))
            std::filesystem::create_directories(dir_);

        std::ofstream ofs(dir_ + "\\" + file_name + ".ply");

        ofs << "ply\n";
        ofs << "format ascii 1.0\n";
        ofs << "element vertex " << particles.size() << "\n";
        ofs << "property float x\n";
        ofs << "property float y\n";
        ofs << "property float z\n";
        ofs << "end_header\n";

        for (const auto &particle: particles) {
            ofs << particle.x << " " << particle.y << " " << particle.z << "\n";
        }

        ofs.close();
    }

    void ModelHelper::export3DModelAsPly(const std::vector<Vec3f> &particles, const std::vector<Vec3f> &colors,
                                         const std::string &dir, const std::string &file_name) {
        auto dir_ = dir;
#if defined(WIN32)
        size_t pos = 0;
        while ((pos = dir_.find('/', pos)) != std::string::npos) {
            dir_.replace(pos, 1, "\\");
            pos += 1;
        }

        // TODO other platform file path transformation
        // ...
#endif

        if (!std::filesystem::exists(dir_))
            std::filesystem::create_directories(dir_);

        std::ofstream ofs(dir_ + "\\" + file_name + ".ply");

        ofs << "ply\n";
        ofs << "format ascii 1.0\n";
        ofs << "element vertex " << particles.size() << "\n";
        ofs << "property float x\n";
        ofs << "property float y\n";
        ofs << "property float z\n";
        ofs << "property uchar red" << std::endl;
        ofs << "property uchar green" << std::endl;
        ofs << "property uchar blue" << std::endl;
        ofs << "end_header\n";

        for (int i = 0; i < particles.size(); ++i) {
            ofs << particles[i].x << " " << particles[i].y << " " << particles[i].z << " ";
            ofs << colors[i].x << " " << colors[i].y << " " << colors[i].z << "\n";
        }

        ofs.close();
    }

    void ModelHelper::export3DModelAsPly(
            const std::vector<Vec3f> &particles,
            const std::vector<Vec2f> &phases,
            const std::string &dir,
            const std::string &file_name) {
        auto dir_ = dir;
#if defined(WIN32)
        size_t pos = 0;
        while ((pos = dir_.find('/', pos)) != std::string::npos) {
            dir_.replace(pos, 1, "\\");
            pos += 1;
        }

        // TODO other platform file path transformation
        // ...
#endif

        if (!std::filesystem::exists(dir_))
            std::filesystem::create_directories(dir_);

        std::ofstream ofs(dir_ + "\\" + file_name + ".ply");

        ofs << "ply\n";
        ofs << "format ascii 1.0\n";
        ofs << "comment  ascii 1.0\n";
        ofs << "element vertex " << particles.size() << "\n";
        ofs << "property float x\n";
        ofs << "property float y\n";
        ofs << "property float z\n";
        ofs << "property float f1" << std::endl;
        ofs << "property float f2" << std::endl;
        ofs << "end_header\n";

        for (int i = 0; i < particles.size(); ++i) {
            ofs << particles[i].x << " " << particles[i].y << " " << particles[i].z << " ";
            ofs << phases[i].x << " " << phases[i].y << "\n";
        }

        ofs.close();
    }

    void ModelHelper::export3DModelAsPly(
            const std::vector<Vec3f> &particles,
            const std::vector<Vec2f> &phases,
            const std::vector<Vec3f> &velocity,
            const std::vector<Vec3f> &velocity_phase,
            const std::vector<Vec3f> &velocity_drift_phase,
            const std::vector<float> &phase_rest_density,
            const std::vector<float> &phase_vis,
            const float &rest_bound_density,
            const float &dt,
            const float &Cf,
            const float &Cd,
            const std::string &dir,
            const std::string &file_name
    ){
        auto dir_ = dir;
        int phase_num = phase_rest_density.size();
#if defined(WIN32)
        size_t pos = 0;
        while ((pos = dir_.find('/', pos)) != std::string::npos) {
            dir_.replace(pos, 1, "\\");
            pos += 1;
        }

        // TODO other platform file path transformation
        // ...
#endif

        if (!std::filesystem::exists(dir_))
            std::filesystem::create_directories(dir_);

        std::ofstream ofs(dir_ + "\\" + file_name + ".ply");

        ofs << "ply\n";
        ofs << "format ascii 1.0\n";
        ofs << "comment phase density:";
        for(auto phase_density: phase_rest_density){
            ofs << " " << phase_density ;
        }
        ofs << "\n";
        ofs << "comment phase viscosity:";
        for(auto phase_viscosity: phase_vis){
            ofs << " " << phase_viscosity ;
        }
        ofs << "\n";
        ofs << "comment rest bound density: " << rest_bound_density << "\n";
        ofs << "comment dt: " << dt << "\n";
        ofs << "comment cf cd: " << Cf << " " << Cd << "\n";
        ofs << "element vertex " << particles.size() << "\n";
        ofs << "property float x\n";
        ofs << "property float y\n";
        ofs << "property float z\n";
        ofs << "property float f1" << std::endl;
        ofs << "property float f2" << std::endl;
        ofs << "property float vx" << std::endl;
        ofs << "property float vy" << std::endl;
        ofs << "property float vz" << std::endl;
        for(int phase_index = 1; phase_index <= phase_num; ++phase_index){
            ofs << "property float vx_phase_" << phase_index << std::endl;
            ofs << "property float vy_phase_" << phase_index << std::endl;
            ofs << "property float vz_phase_" << phase_index << std::endl;
        }
        for(int phase_index = 1; phase_index <= phase_num; ++phase_index){
            ofs << "property float vx_drift_phase_" << phase_index << std::endl;
            ofs << "property float vy_drift_phase_" << phase_index << std::endl;
            ofs << "property float vz_drift_phase_" << phase_index << std::endl;
        }
        ofs << "end_header\n";

        for (int i = 0; i < particles.size(); ++i) {
            ofs << particles[i].x << " " << particles[i].y << " " << particles[i].z << " ";
            ofs << phases[i].x << " " << phases[i].y << " ";
            ofs << velocity[i].x << " " << velocity[i].y << " " << velocity[i].z;
            for(int phase_index = 0; phase_index < phase_num; ++phase_index){
                ofs << " " << velocity_phase[(i * phase_num) + phase_index].x;
                ofs << " " << velocity_phase[(i * phase_num) + phase_index].y;
                ofs << " " << velocity_phase[(i * phase_num) + phase_index].z;
            }
            for(int phase_index = 0; phase_index < phase_num; ++phase_index){
                ofs << " " << velocity_drift_phase[(i * phase_num) + phase_index].x;
                ofs << " " << velocity_drift_phase[(i * phase_num) + phase_index].y;
                ofs << " " << velocity_drift_phase[(i * phase_num) + phase_index].z;
            }
            ofs << "\n";
        }

        ofs.close();
    }

    std::vector<Vec3f> ModelHelper::create3DParticleCube(float particle_radius, Vec3f lb, Vec3f size) {
        std::vector<Vec3f> particles;
        auto diameter = 2 * particle_radius;

        float z = particle_radius + lb.z;
        while (z < lb.z + size.z) {
            float y = particle_radius + lb.y;
            while (y < lb.y + size.y) {
                float x = particle_radius + lb.x;
                while (x < lb.x + size.x) {
                    Vec3f _particles = {x, y, z};
                    particles.push_back(_particles);

                    x += diameter;
                }
                y += diameter;
            }
            z += diameter;
        }

        return particles;
    }

    std::vector<Vec3f> ModelHelper::create3DParticleBox(float particle_radius, Vec3f lb, Vec3f size, float layer) {
        std::vector<Vec3f> particles;

        int numParticles[] = {
                static_cast<int>((size.x + particle_radius) / (2.0 * particle_radius)),
                static_cast<int>((size.y + particle_radius) / (2.0 * particle_radius)),
                static_cast<int>((size.z + particle_radius) / (2.0 * particle_radius))
        };

        for (int i = 0; i <= numParticles[0]; ++i) {
            for (int j = 0; j <= numParticles[1]; ++j) {
                for (int k = 0; k <= numParticles[2]; ++k) {
                    // If this particle is in the first two or last two layers in any dimension...
                    if (i < layer || i >= numParticles[0] - layer || j < layer || j >= numParticles[1] - layer ||
                        k < layer || k >= numParticles[2] - layer) {
                        Vec3f p;
                        p.x = static_cast<float>(lb.x + particle_radius + 2.0 * particle_radius * i);
                        p.y = static_cast<float>(lb.y + particle_radius + 2.0 * particle_radius * j);
                        p.z = static_cast<float>(lb.z + particle_radius + 2.0 * particle_radius * k);
                        particles.push_back(p);
                    }
                }
            }
        }

        return particles;
    }

    std::vector<Vec3f> ModelHelper::create3DParticlePlane(float particle_radius, Vec3f lb, Vec3f size, float layer) {
        std::vector<Vec3f> particles;
        auto diameter = 2 * particle_radius;

        for (float z = particle_radius + lb.z; z < lb.z + size.z; z += diameter) {
            for (float y = particle_radius + lb.y, cnt = 0; cnt < layer; y += diameter, cnt += 1) {
                for (float x = particle_radius + lb.x; x < lb.x + size.x; x += diameter) {
                    Vec3f _particles = {x, y, z};

                    particles.push_back(_particles);
                }
            }
        }

        return particles;
    }

    std::vector<Vec3f> ModelHelper::create3DParticleCylinder(float particleRadius,
                                                             Vec3f center,
                                                             float height,
                                                             float bottom_area_radius) {
        std::vector<Vec3f> particles;
        float diameter = 2 * particleRadius;
        Vec3f topCenter = {center.x, center.y + height / 2, center.z};
        float y0 = topCenter.y;

        for (float y = y0 - particleRadius; y >= y0 - height; y -= diameter) {
            float x0 = topCenter.x - bottom_area_radius;

            for (float x = x0 + particleRadius; x <= topCenter.x + bottom_area_radius; x += diameter) {

                float m_cos = fabs(topCenter.x - x) / bottom_area_radius;
                float length = bottom_area_radius * sqrt(1 - m_cos * m_cos);
                float z0 = topCenter.z - length;
                for (float z = z0 + particleRadius; z <= topCenter.z + length; z += diameter) {
                    Vec3f particle = {x, y, z};
                    particles.push_back(particle);
                }
            }
        }

        return particles;
    }

    std::vector<Vec3f> ModelHelper::create3DParticleSphere(float particle_radius, Vec3f center, float volume_radius) {
        std::vector<Vec3f> particles;
        float gap = particle_radius * 2.0f;

        int num_particles_per_side = std::ceil(volume_radius / gap);
        for (int i = -num_particles_per_side; i <= num_particles_per_side; ++i) {
            for (int j = -num_particles_per_side; j <= num_particles_per_side; ++j) {
                for (int k = -num_particles_per_side; k <= num_particles_per_side; ++k) {
                    Vec3f particle = {float(i) * gap + center.x, float(j) * gap + center.y, float(k) * gap + center.z};

                    if ((particle.x - center.x) * (particle.x - center.x) +
                        (particle.y - center.y) * (particle.y - center.y) +
                        (particle.z - center.z) * (particle.z - center.z) <= volume_radius * volume_radius) {
                        particles.push_back(particle);
                    }
                }
            }
        }

        return particles;
    }

    std::vector<float> ModelHelper::transformVec3fSetToFloatSet(std::vector<Vec3f> &vec3f_set) {
        std::vector<float> particles;
        particles.resize(3 * vec3f_set.size());
        memcpy(particles.data(), vec3f_set.data(), vec3f_set.size() * sizeof(Vec3f));
        return particles;
    }

    std::vector<Vec3f> ModelHelper::loadPly3DModel(std::string ply_file) {
        std::vector<Vec3f> particles;

        Assimp::Importer importer;
        const aiScene *scene = importer.ReadFile(ply_file, aiProcess_Triangulate);

        if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
            std::cerr << "ERROR::ASSIMP::" << importer.GetErrorString() << std::endl;
            return {};
        }

        for (unsigned int m = 0; m < scene->mNumMeshes; m++) {
            const aiMesh *mesh = scene->mMeshes[m];
            for (unsigned int v = 0; v < mesh->mNumVertices; v++) {
                const aiVector3D &vertex = mesh->mVertices[v];
                Vec3f _particles = {vertex.x, vertex.y, vertex.z};
                particles.emplace_back(_particles);
            }
        }

        return particles;
    }

    std::vector<Vec3f> ModelHelper::create3DParticleModel(const ParticleObjectConfig &config) {
        auto *c_ptr = new ParticleObjectConfig;
        memcpy_s(c_ptr, sizeof(ParticleObjectConfig), &config, sizeof(ParticleObjectConfig));
        std::shared_ptr<ParticleObjectConfig> config_ptr(c_ptr);
        return create3DParticleModel(config_ptr);
    }

    Vec3f ModelHelper::loadEmitterAgentNormal(const std::string& agent_file) {
        std::vector<Vec3f> normals;

        Assimp::Importer importer;
        const aiScene *scene = importer.ReadFile(agent_file, aiProcess_Triangulate);

        if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
            std::cerr << "ERROR::ASSIMP::" << importer.GetErrorString() << std::endl;
            return {};
        }

        for (unsigned int m = 0; m < scene->mNumMeshes; m++) {
            const aiMesh *mesh = scene->mMeshes[m];
            for (unsigned int v = 0; v < mesh->mNumVertices; v++) {
                const aiVector3D &normal = mesh->mNormals[v];
                Vec3f _normal = {normal.x, normal.y, normal.z};
                normals.emplace_back(_normal);
            }
        }

        return normals[0];
    }

}  // namespace SoSim
