//
// Created by ADMIN on 2024/3/7.
//

#ifndef SOSIM_DFSPH_PARAMETERS_HPP
#define SOSIM_DFSPH_PARAMETERS_HPP

#include <cuda_runtime.h>

#include "core/data_type.hpp"
#include "core/material.hpp"

namespace SoSim {
    struct DFSPHConstantParams {
        float dt;
        float h;
        float rest_vis;
        float rest_density;
        float rest_volume;
        unsigned particle_num;
        Vec3f gravity;
    };

    struct DFSPHDynamicParams {
        float *density;
        float *mass;
        float *pressure;
        float *volume;
        float *dfsph_alpha;
        float *density_err;
        Vec3f *vel;
        Vec3f *pos;
        Vec3f *acc;
        Vec3f *pos_adv;
        Vec3f *vel_adv;
        Material *mat;

    private:
        bool is_init{false};

    public:
        inline void malloc(unsigned particle_num) {
            if (!is_init) {
                cudaMalloc((void **) &density, particle_num * sizeof(float));
                cudaMalloc((void **) &mass, particle_num * sizeof(float));
                cudaMalloc((void **) &pressure, particle_num * sizeof(float));
                cudaMalloc((void **) &volume, particle_num * sizeof(float));
                cudaMalloc((void **) &dfsph_alpha, particle_num * sizeof(float));
                cudaMalloc((void **) &density_err, particle_num * sizeof(float));
                cudaMalloc((void **) &vel, particle_num * sizeof(Vec3f));
                cudaMalloc((void **) &pos, particle_num * sizeof(Vec3f));
                cudaMalloc((void **) &acc, particle_num * sizeof(Vec3f));
                cudaMalloc((void **) &pos_adv, particle_num * sizeof(Vec3f));
                cudaMalloc((void **) &vel_adv, particle_num * sizeof(Vec3f));
                cudaMalloc((void **) &mat, particle_num * sizeof(Material));

                if (cudaGetLastError() == cudaSuccess)
                    is_init = true;
            }
        }

        inline void freeMemory() {
            if (is_init) {
                cudaFree(density);
                cudaFree(mass);
                cudaFree(pressure);
                cudaFree(volume);
                cudaFree(dfsph_alpha);
                cudaFree(density_err);
                cudaFree(vel);
                cudaFree(pos);
                cudaFree(acc);
                cudaFree(pos_adv);
                cudaFree(vel_adv);
                cudaFree(mat);
            }
        }
    };
}


#endif //SOSIM_DFSPH_PARAMETERS_HPP
