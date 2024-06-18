//
// Created by ADMIN on 2024/3/26.
//

#ifndef SOSIM_IMMCT_PARAMETERS_HPP
#define SOSIM_IMMCT_PARAMETERS_HPP

#include <cuda_runtime.h>

#include "core/math/matrix.hpp"
#include "framework/MetaFramework/sim_material.hpp"

namespace SoSim {
    struct IMMCTConstantParams {
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
        Vec3f phase1_color;
        Vec3f phase2_color;
        Vec2f rest_density;
        float rest_volume{1};
        float rest_rigid_density{0};
        float rest_bound_density{0};
        float sph_h{0};
        float rest_viscosity{0};

        // mixture model
        float Cf{0};

        // JL21
        float kd{0};

        // ISM
        float Cd{0};
        float div_free_threshold{1e-4};
        float incompressible_threshold{1e-4};

        // ISM-CT
        float Cd0{0.2};
        float ct_thinning_exp0{0};
        float solution_vis_base{0};
        float solution_vis_max{0};
        float ct_relaxation_time{0.1};
        float polymer_vol_frac0;
        float vis_bound_damp_factor{0.1};

        // method compare
        float phase1_vis;
        float phase2_vis;
    };

    struct IMMCTDynamicParams {
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

        // mixture model
        Vec2f *vol_frac{nullptr};
        float *rest_density{nullptr};
        Vec3f *vel_phase_1{nullptr};
        Vec3f *vel_phase_2{nullptr};
        Vec3f *acc_phase_1{nullptr};
        Vec3f *acc_phase_2{nullptr};
        Vec3f *vel_drift_phase_1{nullptr};
        Vec3f *vel_drift_phase_2{nullptr};
        float *flag_negative_vol_frac{nullptr};

        // ism
        float *compression_ratio{nullptr};
        float *delta_compression_ratio{nullptr};
        float *df_alpha{nullptr};
        Vec3f *df_alpha_1{nullptr};
        float *df_alpha_2{nullptr};
        float *kappa_div{nullptr};
        float *kappa_incomp{nullptr};
        Vec2f *vol_frac_in{nullptr};
        Vec2f *vol_frac_out{nullptr};

        // ism-ct
        float *Cd{nullptr};
        float *ct_thinning_exp{nullptr};
        Mat33f *CT{nullptr};
        Mat33f *vel_grad{nullptr};
        Mat33f *viscoelastic_stress{nullptr};
        float *shear_rate{nullptr};
        float *solution_vis{nullptr};
        float *ct_vis_increase_exp{nullptr};
        float *blocking_factor{nullptr};

        // method compare
        float *vis{nullptr};

    private:
        bool isInit{false};

    public:
        void malloc(int particle_num) {
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

            // mixture model
            cudaMalloc((void **) &vol_frac, particle_num * sizeof(Vec2f));
            cudaMalloc((void **) &rest_density, particle_num * sizeof(float));
            cudaMalloc((void **) &vel_phase_1, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &vel_phase_2, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &acc_phase_1, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &acc_phase_2, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &vel_drift_phase_1, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &vel_drift_phase_2, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &flag_negative_vol_frac, particle_num * sizeof(float));

            // ims
            cudaMalloc((void **) &compression_ratio, particle_num * sizeof(float));
            cudaMalloc((void **) &delta_compression_ratio, particle_num * sizeof(float));
            cudaMalloc((void **) &df_alpha, particle_num * sizeof(float));
            cudaMalloc((void **) &df_alpha_1, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &df_alpha_2, particle_num * sizeof(float));
            cudaMalloc((void **) &kappa_div, particle_num * sizeof(float));
            cudaMalloc((void **) &kappa_incomp, particle_num * sizeof(float));
            cudaMalloc((void **) &vol_frac_in, particle_num * sizeof(Vec2f));
            cudaMalloc((void **) &vol_frac_out, particle_num * sizeof(Vec2f));

            // ims-ct
            cudaMalloc((void **) &Cd, particle_num * sizeof(float));
            cudaMalloc((void **) &ct_thinning_exp, particle_num * sizeof(float));
            cudaMalloc((void **) &CT, particle_num * sizeof(Mat33f));
            cudaMalloc((void **) &vel_grad, particle_num * sizeof(Mat33f));
            cudaMalloc((void **) &viscoelastic_stress, particle_num * sizeof(Mat33f));
            cudaMalloc((void **) &shear_rate, particle_num * sizeof(float));
            cudaMalloc((void **) &solution_vis, particle_num * sizeof(float));
            cudaMalloc((void **) &ct_vis_increase_exp, particle_num * sizeof(float));
            cudaMalloc((void **) &blocking_factor, particle_num * sizeof(float));

            cudaMalloc((void **) &vis, particle_num * sizeof(float));

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
                cudaFree(vel_phase_1);
                cudaFree(vel_phase_2);
                cudaFree(acc_phase_1);
                cudaFree(acc_phase_2);
                cudaFree(vel_drift_phase_1);
                cudaFree(vel_drift_phase_2);
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
                cudaFree(ct_thinning_exp);
                cudaFree(CT);
                cudaFree(vel_grad);
                cudaFree(viscoelastic_stress);
                cudaFree(shear_rate);
                cudaFree(solution_vis);
                cudaFree(ct_vis_increase_exp);
                cudaFree(blocking_factor);

                cudaFree(vis);

                if (cudaGetLastError() == cudaSuccess)
                    isInit = false;
            }
        }
    };
}

#endif //SOSIM_IMMCT_PARAMETERS_HPP
