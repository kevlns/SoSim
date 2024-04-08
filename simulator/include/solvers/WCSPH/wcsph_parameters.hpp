//
// Created by ADMIN on 2024/3/14.
//

#ifndef SOSIM_WCSPH_PARAMETERS_HPP
#define SOSIM_WCSPH_PARAMETERS_HPP

#include "core/math/matrix.hpp"
#include "framework/MetaFramework/sim_material.hpp"

namespace SoSim {
    struct WCSPHConstantParams {
        float dt;
        float h;
        float rest_vis;
        float rest_density;
        float rest_volume;
        float rest_rigid_density;
        float stiff;
        unsigned particle_num;
        Vec3f gravity;
    };

    struct WCSPHDynamicParams {
        float *m_device_density;
        float *m_device_pressure;
        float *m_device_volume;
        Vec3f *m_device_vel;
        Vec3f *m_device_pos;
        Vec3f *m_device_pos_adv;
        Vec3f *m_device_vel_adv;
        Material *m_device_mat;

    private:
        bool is_init{false};

    public:
        inline void malloc(unsigned particle_num) {
            if (!is_init) {
                cudaMalloc((void **) &m_device_density, particle_num * sizeof(float));
                cudaMalloc((void **) &m_device_pressure, particle_num * sizeof(float));
                cudaMalloc((void **) &m_device_volume, particle_num * sizeof(float));
                cudaMalloc((void **) &m_device_vel, particle_num * sizeof(Vec3f));
                cudaMalloc((void **) &m_device_pos, particle_num * sizeof(Vec3f));
                cudaMalloc((void **) &m_device_pos_adv, particle_num * sizeof(Vec3f));
                cudaMalloc((void **) &m_device_vel_adv, particle_num * sizeof(Vec3f));
                cudaMalloc((void **) &m_device_mat, particle_num * sizeof(Material));

                if (cudaGetLastError() == cudaSuccess)
                    is_init = true;
            }
        }

        inline void freeMemory() {
            if (is_init) {
                cudaFree(m_device_density);
                cudaFree(m_device_pressure);
                cudaFree(m_device_volume);
                cudaFree(m_device_vel);
                cudaFree(m_device_pos);
                cudaFree(m_device_pos_adv);
                cudaFree(m_device_vel_adv);
                cudaFree(m_device_mat);
            }
        }
    };
}

#endif //SOSIM_WCSPH_PARAMETERS_HPP
