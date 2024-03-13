//
// Created by ADMIN on 2024/3/14.
//

#ifndef SOSIM_WCSPH_PARAMETERS_HPP
#define SOSIM_WCSPH_PARAMETERS_HPP

#include "core/data_type.hpp"
#include "core/material.hpp"

namespace SoSim {
    struct WCSPHConstantParams {
        float dt;
        float h;
        float rest_vis;
        float rest_density;
        float rest_volume;
        float rest_rigid_density;
        unsigned particle_num;
        Vec3f gravity;
    };

    struct WCSPHDynamicParams {
        float *density;
        float *pressure;
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
                cudaMalloc((void **) &pressure, particle_num * sizeof(float));
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
                cudaFree(pressure);
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

#endif //SOSIM_WCSPH_PARAMETERS_HPP
