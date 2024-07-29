//
// Created by ADMIN on 2024/3/7.
//

#ifndef SOSIM_IMM_PARAMETERS_v2_HPP
#define SOSIM_IMM_PARAMETERS_v2_HPP

#include <cuda_runtime.h>

#include "core/math/matrix.hpp"
#include "framework/MetaFramework/sim_material.hpp"

namespace SoSim {
    struct IMMConstantParams_v2 {
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
        float rest_volume{1};
        float rest_rigid_density{0};
        float rest_bound_density{0};
        float sph_h{0};
        float rest_viscosity{0.001};

        float div_free_threshold{1e-4};
        float incompressible_threshold{1e-4};

        // imm
        float Cf{0.1};
        float Cd{0.5};
        int phase_num{1};
        float *phase_rest_density;
        Vec3f *phase_color;
        float *phase_vis;
    };

    struct IMMDynamicParams_v2 {
        // common
        Material *mat{nullptr};
        Vec3f *pos{nullptr};
        Vec3f *pos_adv{nullptr};
        Vec3f *vel{nullptr};
        Vec3f *vel_adv{nullptr};
        Vec3f *acc{nullptr};
        Vec3f *surface_normal{nullptr};
        float *mass{nullptr};
        float *volume{nullptr};
        Vec3f *color{nullptr};
        float *density_sph{nullptr};
        bool *is_alive{nullptr};

        float *rest_density{nullptr};
        float *compression_ratio{nullptr};
        float *delta_compression_ratio{nullptr};
        float *df_alpha{nullptr};
        Vec3f *df_alpha_1{nullptr};
        float *df_alpha_2{nullptr};
        float *kappa_div{nullptr};
        float *kappa_incomp{nullptr};

        // imm
        float *vol_frac{nullptr};
        float *vol_frac_in{nullptr};
        float *vol_frac_out{nullptr};
        Vec3f *vel_phase{nullptr};
        Vec3f *acc_phase{nullptr};
        Vec3f *vel_drift_phase{nullptr};

        float *flag_negative_vol_frac{nullptr};

    private:
        bool isInit{false};

    public:
        void malloc(int particle_num, int phase_num = 1) {
            freeMemory();

            // common
            cudaMalloc((void **) &mat, particle_num * sizeof(Material));
            cudaMalloc((void **) &pos, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &pos_adv, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &vel, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &vel_adv, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &acc, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &surface_normal, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &mass, particle_num * sizeof(float));
            cudaMalloc((void **) &volume, particle_num * sizeof(float));
            cudaMalloc((void **) &color, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &density_sph, particle_num * sizeof(float));
            cudaMalloc((void **) &is_alive, particle_num * sizeof(bool));

            cudaMalloc((void **) &rest_density, particle_num * sizeof(float));

            cudaMalloc((void **) &compression_ratio, particle_num * sizeof(float));
            cudaMalloc((void **) &delta_compression_ratio, particle_num * sizeof(float));
            cudaMalloc((void **) &df_alpha, particle_num * sizeof(float));
            cudaMalloc((void **) &df_alpha_1, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &df_alpha_2, particle_num * sizeof(float));
            cudaMalloc((void **) &kappa_div, particle_num * sizeof(float));
            cudaMalloc((void **) &kappa_incomp, particle_num * sizeof(float));

            cudaMalloc((void **) &vol_frac, particle_num * phase_num * sizeof(float));
            cudaMalloc((void **) &vol_frac_in, particle_num * phase_num * sizeof(float));
            cudaMalloc((void **) &vol_frac_out, particle_num * phase_num * sizeof(float));
            cudaMalloc((void **) &vel_phase, particle_num * phase_num * sizeof(Vec3f));
            cudaMalloc((void **) &acc_phase, particle_num * phase_num * sizeof(Vec3f));
            cudaMalloc((void **) &vel_drift_phase, particle_num * phase_num * sizeof(Vec3f));

            cudaMalloc((void **) &flag_negative_vol_frac, particle_num * sizeof(float));

            if (cudaGetLastError() == cudaSuccess)
                isInit = true;
        }

        void freeMemory() {
            if (isInit) {
                cudaFree(mat);
                cudaFree(pos);
                cudaFree(pos_adv);
                cudaFree(vel);
                cudaFree(vel_adv);
                cudaFree(acc);
                cudaFree(surface_normal);
                cudaFree(mass);
                cudaFree(volume);
                cudaFree(color);
                cudaFree(density_sph);
                cudaFree(is_alive);

                cudaFree(rest_density);

                cudaFree(compression_ratio);
                cudaFree(delta_compression_ratio);
                cudaFree(df_alpha);
                cudaFree(df_alpha_1);
                cudaFree(df_alpha_2);
                cudaFree(kappa_div);
                cudaFree(kappa_incomp);

                cudaFree(vol_frac);
                cudaFree(vol_frac_in);
                cudaFree(vol_frac_out);
                cudaFree(vel_phase);
                cudaFree(acc_phase);
                cudaFree(vel_drift_phase);
                cudaFree(flag_negative_vol_frac);


                if (cudaGetLastError() == cudaSuccess)
                    isInit = false;
            }
        }
    };
}


#endif //SOSIM_IMM_PARAMETERS_v2_HPP
