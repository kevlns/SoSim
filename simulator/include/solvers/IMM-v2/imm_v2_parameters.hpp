//
// Created by ADMIN on 2024/3/26.
//

#ifndef SOSIM_IMM_v2_PARAMETERS_HPP
#define SOSIM_IMM_v2_PARAMETERS_HPP

#include <cuda_runtime.h>
#include <vector>

#include "core/math/matrix.hpp"
#include "framework/MetaFramework/sim_material.hpp"

namespace SoSim {
    struct IMMConstantParams_v2 {
        // common
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
        Set<float> phase_rest_densities;
        Set<Vec3f> phase_colors;
        Set<float> phase_vis;
        int phase_num{1};
        float rest_volume{1};
        float rest_rigid_density{0};
        float rest_bound_density{0};
        float sph_h{0};
        float rest_viscosity{0};

        // mixture model
        float Cf{0};

        // ISM
        float Cd{0};
        float div_free_threshold{1e-4};
        float incompressible_threshold{1e-4};
    };

    struct IMMDynamicParams_v2 {
    public:
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

        // imm
        float *vol_frac{nullptr};
        float *rest_density{nullptr};
        Vec3f *vel_phase{nullptr};
        Vec3f *acc_phase{nullptr};
        Vec3f *vel_drift_phase{nullptr};
        float *vis{nullptr};

        float *flag_negative_vol_frac{nullptr};
        float *compression_ratio{nullptr};
        float *delta_compression_ratio{nullptr};
        float *df_alpha{nullptr};
        Vec3f *df_alpha_1{nullptr};
        float *df_alpha_2{nullptr};
        float *kappa_div{nullptr};
        float *kappa_incomp{nullptr};

        float *vol_frac_in{nullptr};
        float *vol_frac_out{nullptr};

        float *Cd{nullptr};
        float *Cf{nullptr}; // TODO

        // Variable-length array
        float* const_phase_rest_densities{nullptr};
        Vec3f* const_phase_colors{nullptr};
        float* const_phase_vis{nullptr};

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

            // imm
            cudaMalloc((void **) &vol_frac, particle_num * sizeof(float) * phase_num);
            cudaMalloc((void **) &vol_frac_in, particle_num * sizeof(float) * phase_num);
            cudaMalloc((void **) &vol_frac_out, particle_num * sizeof(float) * phase_num);
            cudaMalloc((void **) &rest_density, particle_num * sizeof(float));
            cudaMalloc((void **) &vel_phase, particle_num * sizeof(Vec3f) * phase_num);
            cudaMalloc((void **) &acc_phase, particle_num * sizeof(Vec3f) * phase_num);
            cudaMalloc((void **) &vel_drift_phase, particle_num * sizeof(float) * phase_num);
            cudaMalloc((void **) &flag_negative_vol_frac, particle_num * sizeof(float));

            // ims
            cudaMalloc((void **) &compression_ratio, particle_num * sizeof(float));
            cudaMalloc((void **) &delta_compression_ratio, particle_num * sizeof(float));
            cudaMalloc((void **) &df_alpha, particle_num * sizeof(float));
            cudaMalloc((void **) &df_alpha_1, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &df_alpha_2, particle_num * sizeof(float));
            cudaMalloc((void **) &kappa_div, particle_num * sizeof(float));
            cudaMalloc((void **) &kappa_incomp, particle_num * sizeof(float));

            cudaMalloc((void **) &Cd, particle_num * sizeof(float));
//            cudaMalloc((void **) &vis, particle_num * sizeof(float));

            // Variable-length array
            cudaMalloc((void **) &const_phase_rest_densities, sizeof(float) * phase_num);
            cudaMalloc((void **) &const_phase_colors, sizeof(Vec3f) * phase_num);
            cudaMalloc((void **) &const_phase_vis, sizeof(float) * phase_num);

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

                cudaFree(vol_frac);
                cudaFree(rest_density);
                cudaFree(vel_phase);
                cudaFree(acc_phase);
                cudaFree(vel_drift_phase);
                cudaFree(flag_negative_vol_frac);

                cudaFree(compression_ratio);
                cudaFree(delta_compression_ratio);
                cudaFree(df_alpha);
                cudaFree(df_alpha_1);
                cudaFree(df_alpha_2);
                cudaFree(kappa_div);
                cudaFree(kappa_incomp);
                cudaFree(vol_frac_in);
                cudaFree(vol_frac_out);

                cudaFree(Cd);
//                cudaFree(vis);

                cudaFree(const_phase_rest_densities);
                cudaFree(const_phase_colors);
                cudaFree(const_phase_vis);

                if (cudaGetLastError() == cudaSuccess)
                    isInit = false;
            }
        }
    };
}

#endif //SOSIM_IMM_v2_PARAMETERS_HPP
