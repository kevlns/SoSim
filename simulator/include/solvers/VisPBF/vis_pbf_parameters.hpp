//
// Created by ADMIN on 2024/6/13.
//

#ifndef SOSIM_VIS_PBF_PARAMETERS_HPP
#define SOSIM_VIS_PBF_PARAMETERS_HPP

#include <cuda_runtime.h>

#include "core/math/matrix.hpp"
#include "framework/MetaFramework/sim_material.hpp"

namespace SoSim {
    struct VisPBFConstantParams {
        float dt{0};
        float inv_dt{0};
        float inv_dt2{0};
        float cur_sim_time{0};
        Vec3f gravity{0, -9.8f, 0};
        unsigned block_num{1};
        unsigned thread_num{1};
        int particle_num{0};
        float particle_radius{0};

        // sph
        float rest_density{1000};
        float rest_volume{1};
        float rest_rigid_density{0};
        float rest_bound_density{0};
        float sph_h{0};
    };

    struct VisPBFDynamicParams {
        // common
        Material *mat{nullptr};
        Vec3f *pos{nullptr};
        Vec3f *dx{nullptr};
        Vec3f *vel{nullptr};
        Vec3f *color{nullptr};
        Vec3f *error_grad{nullptr};
        float *mass{nullptr};
        float *volume{nullptr};
        float *density_sph{nullptr};
        float *lamb{nullptr};
        float *error{nullptr};
        bool *is_alive{nullptr};

    private:
        bool isInit{false};

    public:
        void malloc(int particle_num) {
            freeMemory();

            // common
            cudaMalloc((void **) &mat, particle_num * sizeof(Material));
            cudaMalloc((void **) &pos, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &dx, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &vel, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &mass, particle_num * sizeof(float));
            cudaMalloc((void **) &volume, particle_num * sizeof(float));
            cudaMalloc((void **) &color, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &error_grad, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &density_sph, particle_num * sizeof(float));
            cudaMalloc((void **) &lamb, particle_num * sizeof(float));
            cudaMalloc((void **) &error, particle_num * sizeof(float));
            cudaMalloc((void **) &is_alive, particle_num * sizeof(bool));

            if (cudaGetLastError() == cudaSuccess)
                isInit = true;
        }

        void freeMemory() {
            if (isInit) {
                cudaFree(mat);
                cudaFree(pos);
                cudaFree(dx);
                cudaFree(vel);
                cudaFree(mass);
                cudaFree(volume);
                cudaFree(color);
                cudaFree(error_grad);
                cudaFree(density_sph);
                cudaFree(lamb);
                cudaFree(error);
                cudaFree(is_alive);

                if (cudaGetLastError() == cudaSuccess)
                    isInit = false;
            }
        }
    };
}

#endif //SOSIM_VIS_PBF_PARAMETERS_HPP
