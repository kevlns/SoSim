//
// Created by ADMIN on 2024/3/26.
//
//

#include "imm_cuda_api.cuh"
#include "imm_macro.hpp"
#include "libs/SPHKernelL/kernels.cuh"
#include "libs/AnalysisL/statistic_util.hpp"

/**
 * cuda impl
 */

namespace SoSim {
    __global__ void
    init_data_cuda(IMMConstantParams *d_const,
                   IMMDynamicParams *d_data,
                   NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // TODO
        DATA_VALUE(flag_negative_vol_frac, p_i) = 0;
        DATA_VALUE(volume, p_i) = CONST_VALUE(rest_volume);
        DATA_VALUE(kappa_div, p_i) = 0;
        DATA_VALUE(color, p_i) = DATA_VALUE(vol_frac, p_i).x * CONST_VALUE(phase1_color) +
                                 DATA_VALUE(vol_frac, p_i).y * CONST_VALUE(phase2_color);
        DATA_VALUE(acc_phase_1, p_i) *= 0;
        DATA_VALUE(Cd, p_i) = CONST_VALUE(Cd0);
        DATA_VALUE(CT, p_i) = Mat33f::eye();
        DATA_VALUE(viscoelastic_stress, p_i) *= 0;
        DATA_VALUE(solution_vis, p_i) = CONST_VALUE(solution_vis_base);
        DATA_VALUE(ct_thinning_exp, p_i) = CONST_VALUE(ct_thinning_exp0);
        DATA_VALUE(ct_vis_increase_exp, p_i) = 1;
    }

    __global__ void
    update_rest_density_and_mass_cuda(IMMConstantParams *d_const,
                                      IMMDynamicParams *d_data,
                                      NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON && DATA_VALUE(mat, p_i) != Emitter_Particle)
            return;

        DATA_VALUE(rest_density, p_i) = dot(DATA_VALUE(vol_frac, p_i), CONST_VALUE(rest_density));
        DATA_VALUE(mass, p_i) = DATA_VALUE(rest_density, p_i) * DATA_VALUE(volume, p_i);
    }

    __global__ void
    update_color_cuda(IMMConstantParams *d_const,
                      IMMDynamicParams *d_data,
                      NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(color, p_i) = DATA_VALUE(vol_frac, p_i).x * CONST_VALUE(phase1_color) +
                                 DATA_VALUE(vol_frac, p_i).y * CONST_VALUE(phase2_color);
    }

    __global__ void
    recover_phase_vel_from_mixture_cuda(IMMConstantParams *d_const,
                                        IMMDynamicParams *d_data,
                                        NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON && DATA_VALUE(mat, p_i) != Emitter_Particle)
            return;

        DATA_VALUE(vel_phase_1, p_i) = DATA_VALUE(vel, p_i);
        DATA_VALUE(vel_phase_2, p_i) = DATA_VALUE(vel, p_i);
    }

    __global__ void
    compute_rigid_volume(IMMConstantParams *d_const,
                         IMMDynamicParams *d_data,
                         NeighborSearchUGConfig *d_nsConfig,
                         NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != FIXED_BOUND && DATA_VALUE(mat, p_i) != DYNAMIC_RIGID &&
            DATA_VALUE(mat, p_i) != STIR_FAN)
            return;

        auto pos_i = DATA_VALUE(pos, p_i);
        float delta = 0;
        FOR_EACH_NEIGHBOR_Pj() {
            auto pos_j = DATA_VALUE(pos, p_j);

            if (DATA_VALUE(mat, p_j) == DATA_VALUE(mat, p_i))
                delta += CUBIC_KERNEL_VALUE();
        }

        DATA_VALUE(volume, p_i) = 1.f / delta;
        DATA_VALUE(rest_density, p_i) = DATA_VALUE(volume, p_i) * CONST_VALUE(rest_bound_density);

        if (DATA_VALUE(mat, p_i) == DYNAMIC_RIGID)
            DATA_VALUE(rest_density, p_i) = DATA_VALUE(volume, p_i) * CONST_VALUE(rest_rigid_density);
    }

    __global__ void
    compute_compression_ratio_cuda(IMMConstantParams *d_const,
                                   IMMDynamicParams *d_data,
                                   NeighborSearchUGConfig *d_nsConfig,
                                   NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        DATA_VALUE(compression_ratio, p_i) *= 0;

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        auto pos_i = DATA_VALUE(pos, p_i);
        FOR_EACH_NEIGHBOR_Pj() {
            if (DATA_VALUE(mat, p_j) == Emitter_Particle)
                continue;

            auto pos_j = DATA_VALUE(pos, p_j);

            DATA_VALUE(compression_ratio, p_i) += DATA_VALUE(volume, p_j) * CUBIC_KERNEL_VALUE();
        }
    }

    __global__ void
    compute_df_beta_cuda(IMMConstantParams *d_const,
                         IMMDynamicParams *d_data,
                         NeighborSearchUGConfig *d_nsConfig,
                         NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        DATA_VALUE(df_alpha_1, p_i) *= 0;
        DATA_VALUE(df_alpha_2, p_i) = 1e-6;

        auto pos_i = DATA_VALUE(pos, p_i);
        FOR_EACH_NEIGHBOR_Pj() {
            if (p_j == p_i || DATA_VALUE(mat, p_j) == Emitter_Particle)
                continue;

            auto pos_j = DATA_VALUE(pos, p_j);
            auto wGrad = CUBIC_KERNEL_GRAD();

            // applied to all dynamic objects
            if (DATA_VALUE(mat, p_i) == COMMON_NEWTON)
                DATA_VALUE(df_alpha_1, p_i) += DATA_VALUE(volume, p_j) * CUBIC_KERNEL_GRAD();

            // applied to all dynamic objects
            if (DATA_VALUE(mat, p_j) == COMMON_NEWTON)
                DATA_VALUE(df_alpha_2, p_i) += dot(wGrad, wGrad) * DATA_VALUE(volume, p_j) * DATA_VALUE(volume, p_j)
                                               / DATA_VALUE(mass, p_j);
        }

        DATA_VALUE(df_alpha, p_i) =
                dot(DATA_VALUE(df_alpha_1, p_i), DATA_VALUE(df_alpha_1, p_i)) / DATA_VALUE(mass, p_i)
                + DATA_VALUE(df_alpha_2, p_i);

        if (DATA_VALUE(df_alpha, p_i) < 1e-6)
            DATA_VALUE(df_alpha, p_i) = 1e-6;
    }

    __global__ void
    compute_delta_compression_ratio_cuda(IMMConstantParams *d_const,
                                         IMMDynamicParams *d_data,
                                         NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(delta_compression_ratio, p_i) = DATA_VALUE(compression_ratio, p_i) - 1.f;
    }

    __global__ void
    update_delta_compression_ratio_from_vel_adv_cuda(IMMConstantParams *d_const,
                                                     IMMDynamicParams *d_data,
                                                     NeighborSearchUGConfig *d_nsConfig,
                                                     NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        auto pos_i = DATA_VALUE(pos, p_i);
        auto vel_adv_i = DATA_VALUE(vel_adv, p_i);
        FOR_EACH_NEIGHBOR_Pj() {
            if (p_j == p_i || DATA_VALUE(mat, p_j) == Emitter_Particle)
                continue;

            auto pos_j = DATA_VALUE(pos, p_j);
            auto vel_adv_j = DATA_VALUE(vel_adv, p_j);
            auto wGrad = CUBIC_KERNEL_GRAD();

            DATA_VALUE(delta_compression_ratio, p_i) += dot(wGrad, vel_adv_i - vel_adv_j) *
                                                        DATA_VALUE(volume, p_j) * CONST_VALUE(dt);
        }

        if (DATA_VALUE(delta_compression_ratio, p_i) < 0)
            DATA_VALUE(delta_compression_ratio, p_i) = 0;
    }

    __global__ void
    compute_kappa_div_from_delta_compression_ratio_cuda(IMMConstantParams *d_const,
                                                        IMMDynamicParams *d_data,
                                                        NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        DATA_VALUE(kappa_div, p_i) *= 0;

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(kappa_div, p_i) = DATA_VALUE(delta_compression_ratio, p_i) / DATA_VALUE(df_alpha, p_i) *
                                     CONST_VALUE(inv_dt2) / DATA_VALUE(volume, p_i);
        DATA_VALUE(df_alpha_2, p_i) += DATA_VALUE(kappa_div, p_i);
    }

    __global__ void
    vf_update_vel_adv_from_kappa_div_cuda(IMMConstantParams *d_const,
                                          IMMDynamicParams *d_data,
                                          NeighborSearchUGConfig *d_nsConfig,
                                          NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        auto pos_i = DATA_VALUE(pos, p_i);
        FOR_EACH_NEIGHBOR_Pj() {
            if (p_j == p_i || DATA_VALUE(mat, p_j) == Emitter_Particle)
                continue;

            auto pos_j = DATA_VALUE(pos, p_j);
            auto wGrad = CUBIC_KERNEL_GRAD();

            DATA_VALUE(vel_adv, p_i) -= CONST_VALUE(dt) * DATA_VALUE(volume, p_i) * DATA_VALUE(volume, p_j) /
                                        DATA_VALUE(mass, p_i) *
                                        (DATA_VALUE(kappa_div, p_i) + DATA_VALUE(kappa_div, p_j)) * wGrad;
        }
    }

    __global__ void
    log_kappa_div_cuda(IMMConstantParams *d_const,
                       IMMDynamicParams *d_data,
                       NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(kappa_div, p_i) = DATA_VALUE(df_alpha_2, p_i);
    }

    __global__ void
    clear_phase_acc_cuda(IMMConstantParams *d_const,
                         IMMDynamicParams *d_data,
                         NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        DATA_VALUE(acc_phase_1, p_i) *= 0;
        DATA_VALUE(acc_phase_2, p_i) *= 0;
    }

    __global__ void
    add_phase_acc_gravity_cuda(IMMConstantParams *d_const,
                               IMMDynamicParams *d_data,
                               NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(acc_phase_1, p_i) += CONST_VALUE(gravity);
        DATA_VALUE(acc_phase_2, p_i) += CONST_VALUE(gravity);
    }

    __global__ void
    add_phase_acc_vis_cuda(IMMConstantParams *d_const,
                           IMMDynamicParams *d_data,
                           NeighborSearchUGConfig *d_nsConfig,
                           NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        Vec3f acc = {0, 0, 0};
        float h2_001 = 0.001f * pow(d_const->sph_h, 2);
        auto pos_i = DATA_VALUE(pos, p_i);
        auto vel_i = DATA_VALUE(vel, p_i);
        FOR_EACH_NEIGHBOR_Pj() {
            if (DATA_VALUE(mat, p_j) != DATA_VALUE(mat, p_i) || p_j == p_i)
                continue;

            auto pos_j = DATA_VALUE(pos, p_j);
            auto x_ij = pos_i - pos_j;
            auto vel_j = DATA_VALUE(vel, p_j);
            auto v_ij = vel_i - vel_j;
            auto wGrad = CUBIC_KERNEL_GRAD();
            auto mass_j = DATA_VALUE(mass, p_j);
            auto vis = CONST_VALUE(rest_viscosity);

            auto pi = -vis * min(0.f, dot(v_ij, pos_i - pos_j)) /
                      (x_ij.length() * x_ij.length() + h2_001);

            acc += -DATA_VALUE(mass, p_i) * mass_j * pi * wGrad;

            if (DATA_VALUE(mat, p_j) == DATA_VALUE(mat, p_i)) {
                auto v_k1_mj = DATA_VALUE(vel_phase_1, p_i) - vel_j;
                DATA_VALUE(acc_phase_1, p_i) += 10 * DATA_VALUE(volume, p_j) *
                                                dot(CONST_VALUE(phase1_vis) * (1 - CONST_VALUE(Cd)) * v_k1_mj +
                                                    (CONST_VALUE(phase1_vis) * CONST_VALUE(Cd) * v_ij),
                                                    x_ij) * wGrad / dot(x_ij, x_ij);

                auto v_k2_mj = DATA_VALUE(vel_phase_2, p_i) - vel_j;
                DATA_VALUE(acc_phase_2, p_i) += 10 * DATA_VALUE(volume, p_j) *
                                                dot(CONST_VALUE(phase2_vis) * (1 - CONST_VALUE(Cd)) * v_k2_mj +
                                                    (CONST_VALUE(phase2_vis) * CONST_VALUE(Cd) * v_ij),
                                                    x_ij) * wGrad / dot(x_ij, x_ij);
            }
        }
    }

    __global__ void
    add_phase_acc_dfsph_vis_cuda(IMMConstantParams *d_const,
                                 IMMDynamicParams *d_data,
                                 NeighborSearchUGConfig *d_nsConfig,
                                 NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        Vec3f acc = {0, 0, 0};
        float h2_001 = 0.001f * pow(CONST_VALUE(sph_h), 2);
        auto pos_i = DATA_VALUE(pos, p_i);
        auto vel_i = DATA_VALUE(vel, p_i);
        auto vis_i = CONST_VALUE(phase1_vis) * DATA_VALUE(vol_frac, p_i).x +
                     CONST_VALUE(phase2_vis) * DATA_VALUE(vol_frac, p_i).y;
        FOR_EACH_NEIGHBOR_Pj() {
            if (DATA_VALUE(mat, p_j) != DATA_VALUE(mat, p_i))
                continue;

            auto pos_j = DATA_VALUE(pos, p_j);
            auto x_ij = pos_i - pos_j;
            auto vel_j = DATA_VALUE(vel, p_j);
            auto v_ij = vel_i - vel_j;
            auto vis_j = CONST_VALUE(phase1_vis) * DATA_VALUE(vol_frac, p_j).x +
                         CONST_VALUE(phase2_vis) * DATA_VALUE(vol_frac, p_j).y;

            auto vis = (vis_i + vis_j) / 2;

            auto pi = vis * DATA_VALUE(mass, p_j) / DATA_VALUE(rest_density, p_j) * dot(v_ij, x_ij) /
                      (x_ij.length() * x_ij.length() + h2_001);

            acc += 10 * pi * CUBIC_KERNEL_GRAD();

        }

        DATA_VALUE(acc_phase_1, p_i) += acc;
        DATA_VALUE(acc_phase_2, p_i) += acc;
    }

    __global__ void
    compute_surface_normal_cuda(IMMConstantParams *d_const,
                                IMMDynamicParams *d_data,
                                NeighborSearchUGConfig *d_nsConfig,
                                NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        auto pos_i = DATA_VALUE(pos, p_i);
        Vec3f normal;
        FOR_EACH_NEIGHBOR_Pj() {
            if (DATA_VALUE(mat, p_j) != DATA_VALUE(mat, p_i))
                continue;

            auto pos_j = DATA_VALUE(pos, p_j);

            normal += CONST_VALUE(sph_h) * DATA_VALUE(mass, p_j) / DATA_VALUE(rest_density, p_j) *
                      cubic_gradient(pos_i - pos_j, CONST_VALUE(sph_h));
        }

        DATA_VALUE(surface_normal, p_i) = normal;
    }

    __global__ void
    add_phase_acc_surface_tension_cuda(IMMConstantParams *d_const,
                                       IMMDynamicParams *d_data,
                                       NeighborSearchUGConfig *d_nsConfig,
                                       NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        Vec3f acc = {0, 0, 0};
        float gamma = 0.005;
        auto pos_i = DATA_VALUE(pos, p_i);
        FOR_EACH_NEIGHBOR_Pj() {
            if (DATA_VALUE(mat, p_j) != DATA_VALUE(mat, p_i) || p_j == p_i)
                continue;

            auto pos_j = DATA_VALUE(pos, p_j);
            auto k =
                    2 * DATA_VALUE(rest_density, p_i) / (DATA_VALUE(rest_density, p_i) + DATA_VALUE(rest_density, p_j));

            auto acc_1 = -gamma * DATA_VALUE(mass, p_i) * DATA_VALUE(mass, p_j) *
                         surface_tension_C((pos_i - pos_j).length(), CONST_VALUE(sph_h)) * (pos_i - pos_j) /
                         (pos_i - pos_j).length();
            auto acc_2 = -gamma * DATA_VALUE(mass, p_i) *
                         (DATA_VALUE(surface_normal, p_i) - DATA_VALUE(surface_normal, p_j));

            acc += k * (acc_1 + acc_2);
        }

        DATA_VALUE(acc_phase_1, p_i) += acc;
        DATA_VALUE(acc_phase_2, p_i) += acc;
    }

    __global__ void
    add_phase_acc_adhesion_cuda(IMMConstantParams *d_const,
                                IMMDynamicParams *d_data,
                                NeighborSearchUGConfig *d_nsConfig,
                                NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        Vec3f acc = {0, 0, 0};
        float beta = 0.0005;
        auto pos_i = DATA_VALUE(pos, p_i);
        FOR_EACH_NEIGHBOR_Pj() {
            if (DATA_VALUE(mat, p_j) != FIXED_BOUND || p_j == p_i)
                continue;

            auto pos_j = DATA_VALUE(pos, p_j);
            auto x_ij = pos_i - pos_j;

            acc += -beta * DATA_VALUE(mass, p_i) * DATA_VALUE(volume, p_j) * adhesion_kernel_value(x_ij,
                                                                                                   CONST_VALUE(sph_h)) *
                   x_ij / x_ij.length();
        }

        DATA_VALUE(acc_phase_1, p_i) += acc;
        DATA_VALUE(acc_phase_2, p_i) += acc;
    }

    __global__ void
    correct_phase_vel_by_vis_bound_cuda(IMMConstantParams *d_const,
                                        IMMDynamicParams *d_data,
                                        NeighborSearchUGConfig *d_nsConfig,
                                        NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        int cnt = 0;
        FOR_EACH_NEIGHBOR_Pj() {
            if (DATA_VALUE(mat, p_j) != FIXED_BOUND || p_j == p_i)
                continue;

            cnt++;
        }

        float f1 = 1.f;
        if (cnt > 80)
            f1 = 0.05;

        DATA_VALUE(vel_phase_1, p_i) = f1 * DATA_VALUE(vel_phase_1, p_i);
        DATA_VALUE(vel_phase_2, p_i) = f1 * DATA_VALUE(vel_phase_2, p_i);
    }

    __global__ void
    phase_acc_2_phase_vel_cuda(IMMConstantParams *d_const,
                               IMMDynamicParams *d_data,
                               NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(vel_phase_1, p_i) += DATA_VALUE(acc_phase_1, p_i) * CONST_VALUE(dt);
        DATA_VALUE(vel_phase_2, p_i) += DATA_VALUE(acc_phase_2, p_i) * CONST_VALUE(dt);
    }

    __global__ void
    update_vel_from_phase_vel_cuda(IMMConstantParams *d_const,
                                   IMMDynamicParams *d_data,
                                   NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(vel, p_i) = DATA_VALUE(vel_phase_1, p_i) * DATA_VALUE(vol_frac, p_i).x
                               + DATA_VALUE(vel_phase_2, p_i) * DATA_VALUE(vol_frac, p_i).y;
        DATA_VALUE(vel_adv, p_i) = DATA_VALUE(vel, p_i);

        DATA_VALUE(vel_drift_phase_1, p_i) = DATA_VALUE(vel_phase_1, p_i) - DATA_VALUE(vel, p_i);
        DATA_VALUE(vel_drift_phase_2, p_i) = DATA_VALUE(vel_phase_2, p_i) - DATA_VALUE(vel, p_i);
    }

    __global__ void
    get_acc_pressure_cuda(IMMConstantParams *d_const,
                          IMMDynamicParams *d_data,
                          NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(acc, p_i) = (DATA_VALUE(vel_adv, p_i) - DATA_VALUE(vel, p_i)) * CONST_VALUE(inv_dt);
    }

    __global__ void
    distribute_acc_pressure_2_phase_cuda(IMMConstantParams *d_const,
                                         IMMDynamicParams *d_data,
                                         NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(acc_phase_1, p_i) += DATA_VALUE(acc, p_i) * (DATA_VALUE(Cd, p_i) + (1 - DATA_VALUE(Cd, p_i)) *
                                                                                      (DATA_VALUE(rest_density, p_i) /
                                                                                       CONST_VALUE(rest_density).x));
        DATA_VALUE(acc_phase_2, p_i) += DATA_VALUE(acc, p_i) * (DATA_VALUE(Cd, p_i) + (1 - DATA_VALUE(Cd, p_i)) *
                                                                                      (DATA_VALUE(rest_density, p_i) /
                                                                                       CONST_VALUE(rest_density).y));
    }

    __global__ void
    compute_kappa_incomp_from_delta_compression_ratio_cuda(IMMConstantParams *d_const,
                                                           IMMDynamicParams *d_data,
                                                           NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(kappa_incomp, p_i) = DATA_VALUE(delta_compression_ratio, p_i) / DATA_VALUE(df_alpha, p_i) *
                                        CONST_VALUE(inv_dt2) / DATA_VALUE(volume, p_i);
        DATA_VALUE(df_alpha_2, p_i) += DATA_VALUE(kappa_incomp, p_i);
    }

    __global__ void
    vf_update_vel_adv_from_kappa_incomp_cuda(IMMConstantParams *d_const,
                                             IMMDynamicParams *d_data,
                                             NeighborSearchUGConfig *d_nsConfig,
                                             NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        auto pos_i = DATA_VALUE(pos, p_i);
        FOR_EACH_NEIGHBOR_Pj() {
            if (p_j == p_i || DATA_VALUE(mat, p_j) == Emitter_Particle)
                continue;

            auto pos_j = DATA_VALUE(pos, p_j);
            auto wGrad = CUBIC_KERNEL_GRAD();

            DATA_VALUE(vel_adv, p_i) -= CONST_VALUE(dt) * DATA_VALUE(volume, p_i) * DATA_VALUE(volume, p_j) /
                                        DATA_VALUE(mass, p_i) *
                                        (DATA_VALUE(kappa_incomp, p_i) + DATA_VALUE(kappa_incomp, p_j)) * wGrad;
        }
    }

    __global__ void
    update_pos_cuda(IMMConstantParams *d_const,
                    IMMDynamicParams *d_data,
                    NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(pos, p_i) += DATA_VALUE(vel, p_i) * CONST_VALUE(dt);
        DATA_VALUE(pos_adv, p_i) = DATA_VALUE(pos, p_i);
    }

    __global__ void
    clear_val_frac_tmp_cuda(IMMConstantParams *d_const,
                            IMMDynamicParams *d_data,
                            NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(vol_frac_in, p_i) *= 0;
        DATA_VALUE(vol_frac_out, p_i) *= 0;
    }

    __global__ void
    update_phase_change_from_drift_cuda(IMMConstantParams *d_const,
                                        IMMDynamicParams *d_data,
                                        NeighborSearchUGConfig *d_nsConfig,
                                        NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        if (DATA_VALUE(flag_negative_vol_frac, p_i) != 0)
            return;

        auto pos_i = DATA_VALUE(pos, p_i);
        FOR_EACH_NEIGHBOR_Pj() {
            if (p_j == p_i || DATA_VALUE(mat, p_j) != DATA_VALUE(mat, p_i) ||
                DATA_VALUE(flag_negative_vol_frac, p_j) != 0)
                continue;

            auto pos_j = DATA_VALUE(pos, p_j);
            auto wGrad = CUBIC_KERNEL_GRAD();

            float vol_frac_change_1 = -CONST_VALUE(dt) * DATA_VALUE(volume, p_j) * dot(DATA_VALUE(vol_frac, p_i).x *
                                                                                       DATA_VALUE(vel_drift_phase_1,
                                                                                                  p_i) +
                                                                                       DATA_VALUE(vol_frac, p_j).x *
                                                                                       DATA_VALUE(vel_drift_phase_1,
                                                                                                  p_j),
                                                                                       wGrad);
            float vol_frac_change_2 = -CONST_VALUE(dt) * DATA_VALUE(volume, p_j) * dot(DATA_VALUE(vol_frac, p_i).y *
                                                                                       DATA_VALUE(vel_drift_phase_2,
                                                                                                  p_i) +
                                                                                       DATA_VALUE(vol_frac, p_j).y *
                                                                                       DATA_VALUE(vel_drift_phase_2,
                                                                                                  p_j),
                                                                                       wGrad);
            if (vol_frac_change_1 < 0)
                DATA_VALUE(vol_frac_out, p_i).x += vol_frac_change_1;
            else
                DATA_VALUE(vol_frac_in, p_i).x += vol_frac_change_1;

            if (vol_frac_change_2 < 0)
                DATA_VALUE(vol_frac_out, p_i).y += vol_frac_change_2;
            else
                DATA_VALUE(vol_frac_in, p_i).y += vol_frac_change_2;
        }
    }

    __global__ void
    update_phase_change_from_diffuse_cuda(IMMConstantParams *d_const,
                                          IMMDynamicParams *d_data,
                                          NeighborSearchUGConfig *d_nsConfig,
                                          NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        if (DATA_VALUE(flag_negative_vol_frac, p_i) != 0)
            return;

        auto pos_i = DATA_VALUE(pos, p_i);
        FOR_EACH_NEIGHBOR_Pj() {
            if (p_j == p_i || DATA_VALUE(mat, p_j) != DATA_VALUE(mat, p_i) ||
                DATA_VALUE(flag_negative_vol_frac, p_j) != 0)
                continue;

            auto pos_j = DATA_VALUE(pos, p_j);
            auto x_ij = pos_i - pos_j;
            auto wGrad = CUBIC_KERNEL_GRAD();
            auto factor = dot(wGrad, x_ij) / dot(x_ij, x_ij);

            float vol_frac_ij_1 = DATA_VALUE(vol_frac, p_i).x - DATA_VALUE(vol_frac, p_j).x;
            float vol_frac_change_1 = CONST_VALUE(dt) * CONST_VALUE(Cf) * vol_frac_ij_1 * DATA_VALUE(volume, p_j) *
                                      factor;

            float vol_frac_ij_2 = DATA_VALUE(vol_frac, p_i).y - DATA_VALUE(vol_frac, p_j).y;
            float vol_frac_change_2 = CONST_VALUE(dt) * CONST_VALUE(Cf) * vol_frac_ij_2 * DATA_VALUE(volume, p_j) *
                                      factor;

            if (vol_frac_change_1 < 0)
                DATA_VALUE(vol_frac_out, p_i).x += vol_frac_change_1;
            else
                DATA_VALUE(vol_frac_in, p_i).x += vol_frac_change_1;

            if (vol_frac_change_2 < 0)
                DATA_VALUE(vol_frac_out, p_i).y += vol_frac_change_2;
            else
                DATA_VALUE(vol_frac_in, p_i).y += vol_frac_change_2;
        }
    }

    __device__ float g_all_positive;

    __global__ void
    check_negative_cuda(IMMConstantParams *d_const,
                        IMMDynamicParams *d_data,
                        NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (p_i == 0)
            g_all_positive = 1;
        __syncthreads();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        if (DATA_VALUE(flag_negative_vol_frac, p_i) != 0)
            return;

        auto vol_frac_tmp = DATA_VALUE(vol_frac, p_i) + DATA_VALUE(vol_frac_out, p_i) + DATA_VALUE(vol_frac_in, p_i);
        if (vol_frac_tmp.x < 0 || vol_frac_tmp.y < 0) {
            DATA_VALUE(flag_negative_vol_frac, p_i) = 1;
            atomicAdd(&g_all_positive, -1);
        }
    }

    __global__ void
    update_phase_change_cuda(IMMConstantParams *d_const,
                             IMMDynamicParams *d_data,
                             NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(vol_frac, p_i) += (DATA_VALUE(vol_frac_in, p_i) + DATA_VALUE(vol_frac_out, p_i));
    }

    __global__ void
    release_unused_drift_vel_cuda(IMMConstantParams *d_const,
                                  IMMDynamicParams *d_data,
                                  NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        if (DATA_VALUE(flag_negative_vol_frac, p_i) != 0) {
            DATA_VALUE(vel_phase_1, p_i) = DATA_VALUE(vel, p_i);
            DATA_VALUE(vel_phase_2, p_i) = DATA_VALUE(vel, p_i);
        }
    }

    __global__ void
    release_negative_cuda(IMMConstantParams *d_const,
                          IMMDynamicParams *d_data,
                          NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(flag_negative_vol_frac, p_i) = 0;
    }

    __global__ void
    regularize_val_frac_cuda(IMMConstantParams *d_const,
                             IMMDynamicParams *d_data,
                             NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        float frac_sum = DATA_VALUE(vol_frac, p_i).x + DATA_VALUE(vol_frac, p_i).y;
        DATA_VALUE(vol_frac, p_i) /= frac_sum;
    }

    __global__ void
    estimate_density_cuda(IMMConstantParams *d_const,
                          IMMDynamicParams *d_data,
                          NeighborSearchUGConfig *d_nsConfig,
                          NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(density_sph, p_i) *= 0;
        auto pos_i = DATA_VALUE(pos, p_i);
        auto rest_dens_i = DATA_VALUE(rest_density, p_i);

        FOR_EACH_NEIGHBOR_Pj() {
            auto pos_j = DATA_VALUE(pos, p_j);
            auto rest_dens_j = DATA_VALUE(rest_density, p_j);
            if (DATA_VALUE(mat, p_j) != COMMON_NEWTON)
                rest_dens_j = rest_dens_i;

            DATA_VALUE(density_sph, p_i) += rest_dens_j * CONST_VALUE(rest_volume) * CUBIC_KERNEL_VALUE();
        }

        DATA_VALUE(density_sph, p_i) = fmax(DATA_VALUE(density_sph, p_i), rest_dens_i);
    }

    __global__ void
    compute_vel_grad_cuda(IMMConstantParams *d_const,
                          IMMDynamicParams *d_data,
                          NeighborSearchUGConfig *d_nsConfig,
                          NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        // TODO
        auto pos_i = DATA_VALUE(pos, p_i);
        auto vel_i = DATA_VALUE(vel, p_i);
        Mat33f vGrad_sum;

        FOR_EACH_NEIGHBOR_Pj() {

            if (DATA_VALUE(mat, p_j) != COMMON_NEWTON)
                continue;

            auto pos_j = DATA_VALUE(pos, p_j);
            auto vel_j = DATA_VALUE(vel, p_j);
            auto vel_ji = vel_j - vel_i;
            auto wGrad = CUBIC_KERNEL_GRAD();
            auto volume_j = DATA_VALUE(volume, p_j);

            vGrad_sum += volume_j * vel_ji * wGrad;
        }

        DATA_VALUE(vel_grad, p_i) = vGrad_sum;
    }

    __global__ void
    update_conformation_tensor_cuda(IMMConstantParams *d_const,
                                    IMMDynamicParams *d_data,
                                    NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        auto CT_i = DATA_VALUE(CT, p_i);

        Mat33f dQ = (CT_i * DATA_VALUE(vel_grad, p_i) +
                     DATA_VALUE(vel_grad, p_i).transpose() * CT_i - 1.f /
                                                                    (CONST_VALUE(ct_relaxation_time) +
                                                                     1e-5f) *
                                                                    (CT_i -
                                                                     Mat33f::eye())) * CONST_VALUE(dt) -
                    DATA_VALUE(ct_thinning_exp, p_i) * (CT_i - Mat33f::eye()) * CT_i;
        DATA_VALUE(CT, p_i) += dQ;
        DATA_VALUE(solution_vis, p_i) = DATA_VALUE(ct_vis_increase_exp, p_i) * CONST_VALUE(solution_vis_base);
        DATA_VALUE(solution_vis, p_i) = fmin(DATA_VALUE(solution_vis, p_i), CONST_VALUE(solution_vis_max));

        auto D = 0.5 * (DATA_VALUE(vel_grad, p_i) + DATA_VALUE(vel_grad, p_i).transpose());
        auto shearRate = sqrtf(0.5f * D.trace() * D.trace());
        DATA_VALUE(solution_vis, p_i) =
                DATA_VALUE(solution_vis, p_i) + (CONST_VALUE(solution_vis_max) - DATA_VALUE(solution_vis, p_i)) /
                                                (1 + pow(10 * shearRate, 5)) * fmin(1.f, DATA_VALUE(vol_frac, p_i).y /
                                                                                         CONST_VALUE(
                                                                                                 polymer_vol_frac0));

        DATA_VALUE(viscoelastic_stress, p_i) =
                DATA_VALUE(solution_vis, p_i) * (DATA_VALUE(CT, p_i) - Mat33f::eye());
    }

    __global__ void
    add_viscoelastic_acc_cuda(IMMConstantParams *d_const,
                              IMMDynamicParams *d_data,
                              NeighborSearchUGConfig *d_nsConfig,
                              NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        auto pos_i = DATA_VALUE(pos, p_i);
        auto dens_i = DATA_VALUE(density_sph, p_i);
        Vec3f acc;

        FOR_EACH_NEIGHBOR_Pj() {

            if (DATA_VALUE(mat, p_j) != COMMON_NEWTON)
                continue;

            auto pos_j = DATA_VALUE(pos, p_j);
            auto wGrad = CUBIC_KERNEL_GRAD();
            auto dens_j = DATA_VALUE(density_sph, p_j);

            acc += (DATA_VALUE(viscoelastic_stress, p_i) / powf(dens_i, 2) +
                    DATA_VALUE(viscoelastic_stress, p_j) / powf(dens_j, 2)) * wGrad;
        }

        DATA_VALUE(acc_phase_1, p_i) += acc * DATA_VALUE(rest_density, p_i) * DATA_VALUE(vol_frac, p_i).y;
        DATA_VALUE(acc_phase_2, p_i) += acc * DATA_VALUE(rest_density, p_i) * DATA_VALUE(vol_frac, p_i).y;
    }

    __global__ void
    update_Cd_cuda(IMMConstantParams *d_const,
                   IMMDynamicParams *d_data,
                   NeighborSearchUGConfig *d_nsConfig,
                   NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        auto pos_i = DATA_VALUE(pos, p_i);
        DATA_VALUE(blocking_factor, p_i) *= 0;
        FOR_EACH_NEIGHBOR_Pj() {
            if (DATA_VALUE(mat, p_j) != DATA_VALUE(mat, p_i) || p_i == p_j)
                continue;

            auto pos_j = DATA_VALUE(pos, p_j);

            DATA_VALUE(blocking_factor, p_i) +=
                    DATA_VALUE(volume, p_j) * DATA_VALUE(vol_frac, p_j).y * CUBIC_KERNEL_VALUE();
        }

        DATA_VALUE(Cd, p_i) = CONST_VALUE(Cd0) + (1 - CONST_VALUE(Cd0)) * DATA_VALUE(blocking_factor, p_i);
    }

    __global__ void
    compute_shear_exp_cuda(IMMConstantParams *d_const,
                           IMMDynamicParams *d_data,
                           NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        // thinning coefficient: beta
        DATA_VALUE(ct_thinning_exp, p_i) = CONST_VALUE(ct_thinning_exp0) * fmax(0.f, CONST_VALUE(polymer_vol_frac0) -
                                                                                     DATA_VALUE(vol_frac, p_i).y) /
                                           (CONST_VALUE(polymer_vol_frac0) + 1e-6);

        // vis increase coefficient: kappa
        DATA_VALUE(ct_vis_increase_exp, p_i) = powf(1.f - DATA_VALUE(vol_frac, p_i).y, -1.2);
    }

}

namespace SoSim { // extra func cuda impl
    __global__ void
    stir_fan_cuda(IMMConstantParams *d_const,
                  IMMDynamicParams *d_data,
                  NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != STIR_FAN)
            return;

        const float M_PI = 3.1415926;
        float angleRadians = -0.004f * (M_PI / 180.0f);// 将角度转换为弧度
        float cosAngle = cos(angleRadians);
        float sinAngle = sin(angleRadians);

        Vec3f center_offset = {0, 0, 0};

        auto pos = DATA_VALUE(pos, p_i) - center_offset;
        DATA_VALUE(pos, p_i).x = pos.x * cosAngle - pos.z * sinAngle + center_offset.x;
        DATA_VALUE(pos, p_i).z = pos.x * sinAngle + pos.z * cosAngle + center_offset.z;
    }

    __global__ void
    buckling_move_cuda(IMMConstantParams *d_const,
                       IMMDynamicParams *d_data,
                       NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != MOVING_TUBE && DATA_VALUE(mat, p_i) != MOVING_COVER)
            return;

        DATA_VALUE(pos, p_i) += CONST_VALUE(dt) * DATA_VALUE(vel, p_i);
    }

    __global__ void
    rotate_bowl_cuda(IMMConstantParams *d_const,
                     IMMDynamicParams *d_data,
                     NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != MOVING_BOWL && DATA_VALUE(mat, p_i) != STIR_FAN)
            return;

        const float M_PI = 3.1415926;
        float angleRadians = -0.0001;// 将角度转换为弧度
        float cosAngle = cos(angleRadians);
        float sinAngle = sin(angleRadians);

        Vec3f offset = {0, -6.8165, 0};

        auto pos = DATA_VALUE(pos, p_i) - offset;
        DATA_VALUE(pos, p_i).y = pos.y * cosAngle - pos.z * sinAngle + offset.y;
        DATA_VALUE(pos, p_i).z = pos.y * sinAngle + pos.z * cosAngle + offset.z;
    }

    __global__ void
    correct_phase_vel_by_artificial_vis_bound_cuda(IMMConstantParams *d_const,
                                                   IMMDynamicParams *d_data,
                                                   NeighborSearchUGConfig *d_nsConfig,
                                                   NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        int cnt = 0;
        FOR_EACH_NEIGHBOR_Pj() {
            if (DATA_VALUE(mat, p_j) == DATA_VALUE(mat, p_i) || DATA_VALUE(mat, p_j) == Emitter_Particle)
                continue;

            cnt++;
        }

        float f1 = 1;
        if (cnt > 15)
            f1 = 1;

        DATA_VALUE(vel_phase_1, p_i) *= f1;
        DATA_VALUE(vel_phase_2, p_i) *= f1;
    }
}


/**
 * host invoke impl
 */

namespace SoSim {
    __host__ void
    init_data(IMMConstantParams &h_const,
              IMMConstantParams *d_const,
              IMMDynamicParams *d_data,
              NeighborSearchUGParams *d_nsParams) {
        init_data_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    prepare_ims(IMMConstantParams &h_const,
                IMMConstantParams *d_const,
                IMMDynamicParams *d_data,
                NeighborSearchUGConfig *d_nsConfig,
                NeighborSearchUGParams *d_nsParams) {
        // ims update_rest_density_and_mass()
        update_rest_density_and_mass_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // ims recover_phase_vel_from_mixture()
        recover_phase_vel_from_mixture_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // compute_rigid_volume()
        compute_rigid_volume<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);
    }

    __host__ void
    sph_precompute(IMMConstantParams &h_const,
                   IMMConstantParams *d_const,
                   IMMDynamicParams *d_data,
                   NeighborSearchUGConfig *d_nsConfig,
                   NeighborSearchUGParams *d_nsParams) {
        // compute_compression_ratio(), AKA step_sph_compute_compression_ratio()
        compute_compression_ratio_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);

        // compute_df_beta(), AKA step_df_compute_beta()
        compute_df_beta_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);
    }

    __host__ void
    vfsph_div(IMMConstantParams &h_const,
              IMMDynamicParams &h_data,
              Vec3ui &obj_part_index,
              IMMConstantParams *d_const,
              IMMDynamicParams *d_data,
              NeighborSearchUGConfig *d_nsConfig,
              NeighborSearchUGParams *d_nsParams,
              bool &crash) {
        int iter = 0;
        while (true) {
            iter++;

            // compute_delta_compression_ratio()
            compute_delta_compression_ratio_cuda<<<h_const.block_num, h_const.thread_num>>>(
                    d_const, d_data, d_nsParams);

            std::vector<float> de(h_const.particle_num);
            cudaMemcpy(de.data(), h_data.delta_compression_ratio, h_const.particle_num * sizeof(float),
                       cudaMemcpyDeviceToHost);

            // update_delta_compression_ratio_from_vel_adv()
            update_delta_compression_ratio_from_vel_adv_cuda<<<h_const.block_num, h_const.thread_num>>>(
                    d_const, d_data, d_nsConfig, d_nsParams);

            cudaMemcpy(de.data(), h_data.delta_compression_ratio, h_const.particle_num * sizeof(float),
                       cudaMemcpyDeviceToHost);

            // update_vf_compressible_ratio()
            auto compressible_ratio = cal_mean(h_data.delta_compression_ratio,
                                               h_const.particle_num, obj_part_index.y);

            // compute_kappa_div_from_delta_compression_ratio()
            compute_kappa_div_from_delta_compression_ratio_cuda<<<h_const.block_num, h_const.thread_num>>>(
                    d_const, d_data, d_nsParams);

            std::vector<float> kd(h_const.particle_num);
            cudaMemcpy(kd.data(), h_data.kappa_div, h_const.particle_num * sizeof(float),
                       cudaMemcpyDeviceToHost);

            std::vector<float> mass(h_const.particle_num);
            cudaMemcpy(mass.data(), h_data.mass, h_const.particle_num * sizeof(float),
                       cudaMemcpyDeviceToHost);

            // vf_update_vel_adv_from_kappa_div()
            vf_update_vel_adv_from_kappa_div_cuda<<<h_const.block_num, h_const.thread_num>>>(
                    d_const, d_data, d_nsConfig, d_nsParams);

            // check compressible_ratio
            if (compressible_ratio < h_const.div_free_threshold || iter > 100)
                break;
        }

        // log_kappa_div()
//        log_kappa_div_cuda<<<h_const.block_num, h_const.thread_num>>>(
//                d_const, d_data, d_nsParams);

        std::cout << "div-iter: " << iter << '\n';

        if (iter == 101)
            crash = true;

        // vel = vel_adv
        //        cudaMemcpy(h_data.vel, h_data.vel_adv, h_const.particle_num * sizeof(Vec3f), cudaMemcpyDeviceToDevice);
    }

    __host__ void
    apply_pressure_acc(IMMConstantParams &h_const,
                       IMMConstantParams *d_const,
                       IMMDynamicParams *d_data,
                       NeighborSearchUGParams *d_nsParams) {
        // get_acc_pressure()
        get_acc_pressure_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // clear_phase_acc()
        clear_phase_acc_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // distribute_acc_pressure_2_phase()
        distribute_acc_pressure_2_phase_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // phase_acc_2_phase_vel()
        phase_acc_2_phase_vel_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // update_vel_from_phase_vel()
        update_vel_from_phase_vel_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    ism_gravity_vis_surface(IMMConstantParams &h_const,
                            IMMConstantParams *d_const,
                            IMMDynamicParams *d_data,
                            NeighborSearchUGConfig *d_nsConfig,
                            NeighborSearchUGParams *d_nsParams) {
        // clear_phase_acc()
        clear_phase_acc_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // add_phase_acc_gravity()
        add_phase_acc_gravity_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // add_phase_acc_vis()
        add_phase_acc_vis_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);

        // compute_surface_normal
        compute_surface_normal_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);

        // add_phase_acc_surface_tension_cuda()
        add_phase_acc_surface_tension_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);

        // phase_acc_2_phase_vel()
        phase_acc_2_phase_vel_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // update_vel_from_phase_vel()
        update_vel_from_phase_vel_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    dfsph_gravity_vis_surface(IMMConstantParams &h_const,
                              IMMConstantParams *d_const,
                              IMMDynamicParams *d_data,
                              NeighborSearchUGConfig *d_nsConfig,
                              NeighborSearchUGParams *d_nsParams) {
        // clear_phase_acc()
        clear_phase_acc_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // add_phase_acc_gravity()
        add_phase_acc_gravity_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // add_phase_acc_vis()
        add_phase_acc_dfsph_vis_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);

        // compute_surface_normal
        compute_surface_normal_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);

        // add_phase_acc_surface_tension_cuda()
        add_phase_acc_surface_tension_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);

        // phase_acc_2_phase_vel()
        phase_acc_2_phase_vel_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // update_vel_from_phase_vel()
        update_vel_from_phase_vel_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    vfsph_incomp(IMMConstantParams &h_const,
                 IMMDynamicParams &h_data,
                 Vec3ui &obj_part_index,
                 IMMConstantParams *d_const,
                 IMMDynamicParams *d_data,
                 NeighborSearchUGConfig *d_nsConfig,
                 NeighborSearchUGParams *d_nsParams,
                 bool &crash) {
        int iter = 0;
        while (true) {
            iter++;

            // compute_delta_compression_ratio()
            compute_delta_compression_ratio_cuda<<<h_const.block_num, h_const.thread_num>>>(
                    d_const, d_data, d_nsParams);

            // update_delta_compression_ratio_from_vel_adv()
            update_delta_compression_ratio_from_vel_adv_cuda<<<h_const.block_num, h_const.thread_num>>>(
                    d_const, d_data, d_nsConfig, d_nsParams);

            // update_vf_compressible_ratio()
            auto compressible_ratio = cal_mean(h_data.delta_compression_ratio,
                                               h_const.particle_num, obj_part_index.y);

            // compute_kappa_incomp_from_delta_compression_ratio()
            compute_kappa_incomp_from_delta_compression_ratio_cuda<<<h_const.block_num, h_const.thread_num>>>(
                    d_const, d_data, d_nsParams);

            // vf_update_vel_adv_from_kappa_incomp()
            vf_update_vel_adv_from_kappa_incomp_cuda<<<h_const.block_num, h_const.thread_num>>>(
                    d_const, d_data, d_nsConfig, d_nsParams);

            // check compressible_ratio
            if (compressible_ratio < h_const.incompressible_threshold || iter > 100)
                break;
        }

        std::cout << "incomp-iter: " << iter << '\n';

        if (iter == 101)
            crash = true;
    }

    __host__ void
    update_pos(IMMConstantParams &h_const,
               IMMConstantParams *d_const,
               IMMDynamicParams *d_data,
               NeighborSearchUGParams *d_nsParams) {
        update_pos_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    phase_transport_ism(IMMConstantParams &h_const,
                        IMMConstantParams *d_const,
                        IMMDynamicParams *d_data,
                        NeighborSearchUGConfig *d_nsConfig,
                        NeighborSearchUGParams *d_nsParams,
                        bool &crash) {
        // clear_val_frac_tmp()
        clear_val_frac_tmp_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // update_phase_change_from_drift()
        update_phase_change_from_drift_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);

        // update_phase_change_from_diffuse()
        update_phase_change_from_diffuse_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);

        // while: check_negative(), update_phase_change_from_drift(), update_phase_change_from_diffuse()
        float all_positive = 0;
        int iter = 1;
        while (true) {
            // check
            check_negative_cuda<<<h_const.block_num, h_const.thread_num>>>(
                    d_const, d_data, d_nsParams);
            cudaMemcpyFromSymbol(&all_positive, g_all_positive, sizeof(float), 0, cudaMemcpyDeviceToHost);
            if (all_positive == 1 || iter > 100)
                break;

            // clear_val_frac_tmp()
            clear_val_frac_tmp_cuda<<<h_const.block_num, h_const.thread_num>>>(
                    d_const, d_data, d_nsParams);

            // update_phase_change_from_drift()
            update_phase_change_from_drift_cuda<<<h_const.block_num, h_const.thread_num>>>(
                    d_const, d_data, d_nsConfig, d_nsParams);

            // update_phase_change_from_diffuse()
            update_phase_change_from_diffuse_cuda<<<h_const.block_num, h_const.thread_num>>>(
                    d_const, d_data, d_nsConfig, d_nsParams);

            iter++;
        }

        std::cout << "phase-trans-iter: " << iter << "\n";

        if (iter == 101)
            crash = true;

        // update_phase_change()
        update_phase_change_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // release_unused_drift_vel()
        release_unused_drift_vel_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // release_negative()
        release_negative_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // regularize_val_frac()
        regularize_val_frac_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // update_rest_density_and_mass()
        update_rest_density_and_mass_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // update_vel_from_phase_vel()
        update_vel_from_phase_vel_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    update_mass_and_vel(IMMConstantParams &h_const,
                        IMMConstantParams *d_const,
                        IMMDynamicParams *d_data,
                        NeighborSearchUGParams *d_nsParams) {
        // regularize_val_frac()
        regularize_val_frac_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // update_rest_density_and_mass()
        update_rest_density_and_mass_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // update_vel_from_phase_vel()
        update_vel_from_phase_vel_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    update_color(IMMConstantParams &h_const,
                 IMMConstantParams *d_const,
                 IMMDynamicParams *d_data,
                 NeighborSearchUGParams *d_nsParams) {
        // get_acc_pressure()
        update_color_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    ism_viscoelastic(IMMConstantParams &h_const,
                     IMMConstantParams *d_const,
                     IMMDynamicParams *d_data,
                     NeighborSearchUGConfig *d_nsConfig,
                     NeighborSearchUGParams *d_nsParams) {
        // estimate_density()
        estimate_density_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);

        // compute_vel_grad()
        compute_vel_grad_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);

        // update_conformation_tensor()
        update_conformation_tensor_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // clear_phase_acc()
        clear_phase_acc_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        add_viscoelastic_acc_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);

        // phase_acc_2_phase_vel()
        phase_acc_2_phase_vel_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // update_vel_from_phase_vel()
        update_vel_from_phase_vel_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    update_CT_parameters(IMMConstantParams &h_const,
                         IMMConstantParams *d_const,
                         IMMDynamicParams *d_data,
                         NeighborSearchUGConfig *d_nsConfig,
                         NeighborSearchUGParams *d_nsParams) {
        // set_Cd()
//        update_Cd_cuda<<<h_const.block_num, h_const.thread_num>>>(
//                d_const, d_data, d_nsConfig, d_nsParams);

        // compute_shear_exp()
        compute_shear_exp_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    artificial_vis_bound(IMMConstantParams &h_const,
                         IMMConstantParams *d_const,
                         IMMDynamicParams *d_data,
                         NeighborSearchUGConfig *d_nsConfig,
                         NeighborSearchUGParams *d_nsParams) {
        // correct_phase_vel_by_artificial_vis_bound()
        correct_phase_vel_by_artificial_vis_bound_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);

        // update_vel_from_phase_vel()
        update_vel_from_phase_vel_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }
}

namespace SoSim { // extra func host invoke
    __host__ void
    stirring(IMMConstantParams &h_const,
             IMMConstantParams *d_const,
             IMMDynamicParams *d_data,
             NeighborSearchUGParams *d_nsParams) {
        stir_fan_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    rotate_bowl(IMMConstantParams &h_const,
                IMMConstantParams *d_const,
                IMMDynamicParams *d_data,
                NeighborSearchUGParams *d_nsParams) {
        rotate_bowl_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    buckling(IMMConstantParams &h_const,
             IMMConstantParams *d_const,
             IMMDynamicParams *d_data,
             NeighborSearchUGParams *d_nsParams) {
        buckling_move_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }
}