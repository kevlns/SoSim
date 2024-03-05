//
// Created by ADMIN on 2024/3/5.
//

#include "libs/ModelL/model_loader.hpp"
#include "assimp/Importer.hpp"
#include "assimp/postprocess.h"
#include "assimp/scene.h"

namespace SoSim {

    std::vector<Vec3f> ModelLoader::createParticleCube(float particle_radius, Vec3f lb, Vec3f size) {
        std::vector<Vec3f> pos;
        auto diameter = 2 * particle_radius;

        for (float z = particle_radius + lb.z; z < lb.z + size.z; z += diameter) {
            for (float y = particle_radius + lb.y; y < lb.y + size.y; y += diameter) {
                for (float x = particle_radius + lb.x; x < lb.x + size.x; x += diameter) {
                    Vec3f _pos = {x, y, z};

                    pos.push_back(_pos);
                }
            }
        }

        return pos;
    }

    std::vector<Vec3f> ModelLoader::createParticleBox(float particle_radius, Vec3f lb, Vec3f size, float layer) {
        std::vector<Vec3f> pos;

        int numParticles[] = {
                static_cast<int>(size.x / (2.0 * particle_radius)),
                static_cast<int>(size.y / (2.0 * particle_radius)),
                static_cast<int>(size.z / (2.0 * particle_radius))};

        for (int i = 0; i < numParticles[0]; ++i) {
            for (int j = 0; j < numParticles[1]; ++j) {
                for (int k = 0; k < numParticles[2]; ++k) {
                    // If this particle is in the first two or last two layers in any dimension...
                    if (i < layer || i >= numParticles[0] - layer || j < layer || j >= numParticles[1] - layer ||
                        k < layer ||
                        k >= numParticles[2] - layer) {
                        Vec3f p;
                        p.x = lb.x + particle_radius + 2.0 * particle_radius * i;
                        p.y = lb.y + particle_radius + 2.0 * particle_radius * j;
                        p.z = lb.z + particle_radius + 2.0 * particle_radius * k;
                        pos.push_back(p);
                    }
                }
            }
        }

        return pos;
    }

    std::vector<Vec3f> ModelLoader::createParticlePlane(float particle_radius, Vec3f lb, Vec3f size, float layer) {
        std::vector<Vec3f> pos;
        auto diameter = 2 * particle_radius;

        for (float z = particle_radius + lb.z; z < lb.z + size.z; z += diameter) {
            for (float y = particle_radius + lb.y, cnt = 0; cnt < layer; y += diameter, cnt += 1) {
                for (float x = particle_radius + lb.x; x < lb.x + size.x; x += diameter) {
                    Vec3f _pos = {x, y, z};

                    pos.push_back(_pos);
                }
            }
        }

        return pos;
    }

    std::vector<Vec3f> ModelLoader::loadParticle3DModel(std::string model_file) {
        std::vector<Vec3f> pos;

        Assimp::Importer importer;
        const aiScene *scene = importer.ReadFile(model_file, aiProcess_Triangulate);

        if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
            std::cerr << "ERROR::ASSIMP::" << importer.GetErrorString() << std::endl;
            return {};
        }

        for (unsigned int m = 0; m < scene->mNumMeshes; m++) {
            const aiMesh *mesh = scene->mMeshes[m];
            for (unsigned int v = 0; v < mesh->mNumVertices; v++) {
                const aiVector3D &vertex = mesh->mVertices[v];
                Vec3f _pos = {vertex.x, vertex.y, vertex.z};
                pos.emplace_back(_pos);
            }
        }

        return pos;
    }
}
