//
// Created by ADMIN on 2024/3/26.
//
//

#include "ism_ct_cuda_api.cuh"
#include "ism_ct_macro.hpp"
#include "libs/SPHKernelL/kernels.cuh"
#include "libs/AnalysisL/statistic_util.hpp"

/**
 * cuda impl
 */

namespace SoSim {
    __global__ void
    init_data_cuda(IMSCTConstantParams *d_const,
                   IMSCTDynamicParams *d_data,
                   NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // TODO
        DATA_VALUE(flag_negative_vol_frac, p_i) = 0;
        DATA_VALUE(volume, p_i) = CONST_VALUE(rest_volume);
        DATA_VALUE(kappa_div, p_i) = 0;
        DATA_VALUE(color, p_i) = DATA_VALUE(vol_frac, p_i).x * CONST_VALUE(phase1_color) +
                                 DATA_VALUE(vol_frac, p_i).y * CONST_VALUE(phase2_color);
        DATA_VALUE(acc_phase_1, p_i) *= 0;
//        DATA_VALUE(Cd, p_i) = CONST_VALUE(Cd0);
//        DATA_VALUE(CT, p_i) = Mat33f::eye();
//        DATA_VALUE(viscoelastic_stress, p_i) *= 0;
//        DATA_VALUE(solution_vis, p_i) = CONST_VALUE(solution_vis0);
    }

    __global__ void
    update_rest_density_and_mass_cuda(IMSCTConstantParams *d_const,
                                      IMSCTDynamicParams *d_data,
                                      NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
            return;

        DATA_VALUE(rest_density, p_i) = dot(DATA_VALUE(vol_frac, p_i), CONST_VALUE(rest_density));
        DATA_VALUE(mass, p_i) = DATA_VALUE(rest_density, p_i) * DATA_VALUE(volume, p_i);
    }

    __global__ void
    update_color_cuda(IMSCTConstantParams *d_const,
                      IMSCTDynamicParams *d_data,
                      NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
            return;

        DATA_VALUE(color, p_i) = DATA_VALUE(vol_frac, p_i).x * CONST_VALUE(phase1_color) +
                                 DATA_VALUE(vol_frac, p_i).y * CONST_VALUE(phase2_color);
    }

    __global__ void
    recover_phase_vel_from_mixture_cuda(IMSCTConstantParams *d_const,
                                        IMSCTDynamicParams *d_data,
                                        NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
            return;

        DATA_VALUE(vel_phase_1, p_i) = DATA_VALUE(vel, p_i);
        DATA_VALUE(vel_phase_2, p_i) = DATA_VALUE(vel, p_i);
    }

    __global__ void
    compute_rigid_volume(IMSCTConstantParams *d_const,
                         IMSCTDynamicParams *d_data,
                         NeighborSearchUGConfig *d_nsConfig,
                         NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != FIXED_BOUND && DATA_VALUE(mat, p_i) != DYNAMIC_RIGID)
            return;

        auto pos_i = DATA_VALUE(pos, p_i);
        float delta = 0;
        FOR_EACH_NEIGHBOR_Pj() {
            auto pos_j = DATA_VALUE(pos, p_j);

            if (DATA_VALUE(mat, p_j) == DATA_VALUE(mat, p_i))
                delta += CUBIC_KERNEL_VALUE();
        }

        DATA_VALUE(volume, p_i) = 1.f / delta;
    }

    __global__ void
    compute_compression_ratio_cuda(IMSCTConstantParams *d_const,
                                   IMSCTDynamicParams *d_data,
                                   NeighborSearchUGConfig *d_nsConfig,
                                   NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
            return;

        DATA_VALUE(compression_ratio, p_i) = 0;

        auto pos_i = DATA_VALUE(pos, p_i);
        int neighbor_cnt = 0;
        FOR_EACH_NEIGHBOR_Pj() {
            auto pos_j = DATA_VALUE(pos, p_j);

            DATA_VALUE(compression_ratio, p_i) += DATA_VALUE(volume, p_j) * CUBIC_KERNEL_VALUE();
            neighbor_cnt++;
        }

        if (neighbor_cnt < 35)
            DATA_VALUE(compression_ratio, p_i) = fmax(DATA_VALUE(compression_ratio, p_i), DATA_VALUE(volume, p_i));
    }

    __global__ void
    compute_df_beta_cuda(IMSCTConstantParams *d_const,
                         IMSCTDynamicParams *d_data,
                         NeighborSearchUGConfig *d_nsConfig,
                         NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        DATA_VALUE(df_alpha_1, p_i) *= 0;
        DATA_VALUE(df_alpha_2, p_i) = 1e-6;

        auto pos_i = DATA_VALUE(pos, p_i);
        FOR_EACH_NEIGHBOR_Pj() {
            if (p_j == p_i)
                continue;

            auto pos_j = DATA_VALUE(pos, p_j);
            auto wGrad = CUBIC_KERNEL_GRAD();

            // applied to all dynamic objects
            if (DATA_VALUE(mat, p_i) == IMSCT_NONNEWTON)
                DATA_VALUE(df_alpha_1, p_i) += DATA_VALUE(volume, p_j) * CUBIC_KERNEL_GRAD();

            // applied to all dynamic objects
            if (DATA_VALUE(mat, p_j) == IMSCT_NONNEWTON)
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
    compute_delta_compression_ratio_cuda(IMSCTConstantParams *d_const,
                                         IMSCTDynamicParams *d_data,
                                         NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
            return;

        DATA_VALUE(delta_compression_ratio, p_i) = DATA_VALUE(compression_ratio, p_i) - 1.f;
    }

    __global__ void
    update_delta_compression_ratio_from_vel_adv_cuda(IMSCTConstantParams *d_const,
                                                     IMSCTDynamicParams *d_data,
                                                     NeighborSearchUGConfig *d_nsConfig,
                                                     NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
            return;

        auto pos_i = DATA_VALUE(pos, p_i);
        auto vel_adv_i = DATA_VALUE(vel_adv, p_i);
        FOR_EACH_NEIGHBOR_Pj() {
            if (p_j == p_i)
                continue;

            auto pos_j = DATA_VALUE(pos, p_j);
            auto vel_adv_j = DATA_VALUE(vel_adv, p_j);
            auto wGrad = CUBIC_KERNEL_GRAD();

            DATA_VALUE(delta_compression_ratio, p_i) += dot(wGrad, vel_adv_i - vel_adv_j) *
                                                        DATA_VALUE(volume, p_j) * CONST_VALUE(dt);
        }

        if (DATA_VALUE(delta_compression_ratio, p_i) < 1e-6)
            DATA_VALUE(delta_compression_ratio, p_i) = 1e-6;
    }

    __global__ void
    compute_kappa_div_from_delta_compression_ratio_cuda(IMSCTConstantParams *d_const,
                                                        IMSCTDynamicParams *d_data,
                                                        NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
            return;

        DATA_VALUE(kappa_div, p_i) = DATA_VALUE(delta_compression_ratio, p_i) / DATA_VALUE(df_alpha, p_i) *
                                     CONST_VALUE(inv_dt2) / DATA_VALUE(volume, p_i);
        DATA_VALUE(df_alpha_2, p_i) += DATA_VALUE(kappa_div, p_i);
    }

    __global__ void
    vf_update_vel_adv_from_kappa_div_cuda(IMSCTConstantParams *d_const,
                                          IMSCTDynamicParams *d_data,
                                          NeighborSearchUGConfig *d_nsConfig,
                                          NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
            return;

        auto pos_i = DATA_VALUE(pos, p_i);
        FOR_EACH_NEIGHBOR_Pj() {
            if (p_j == p_i)
                continue;

            auto pos_j = DATA_VALUE(pos, p_j);
            auto wGrad = CUBIC_KERNEL_GRAD();

            DATA_VALUE(vel_adv, p_i) -= CONST_VALUE(dt) * DATA_VALUE(volume, p_i) * DATA_VALUE(volume, p_j) /
                                        DATA_VALUE(mass, p_i) *
                                        (DATA_VALUE(kappa_div, p_i) + DATA_VALUE(kappa_div, p_j)) * wGrad;
        }
    }

    __global__ void
    log_kappa_div_cuda(IMSCTConstantParams *d_const,
                       IMSCTDynamicParams *d_data,
                       NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
            return;

        DATA_VALUE(kappa_div, p_i) = DATA_VALUE(df_alpha_2, p_i);
    }

    __global__ void
    clear_phase_acc_cuda(IMSCTConstantParams *d_const,
                         IMSCTDynamicParams *d_data,
                         NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        DATA_VALUE(acc_phase_1, p_i) *= 0;
        DATA_VALUE(acc_phase_2, p_i) *= 0;
    }

    __global__ void
    add_phase_acc_gravity_cuda(IMSCTConstantParams *d_const,
                               IMSCTDynamicParams *d_data,
                               NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
            return;

        DATA_VALUE(acc_phase_1, p_i) += CONST_VALUE(gravity);
        DATA_VALUE(acc_phase_2, p_i) += CONST_VALUE(gravity);
    }

    __global__ void
    add_phase_acc_vis_cuda(IMSCTConstantParams *d_const,
                           IMSCTDynamicParams *d_data,
                           NeighborSearchUGConfig *d_nsConfig,
                           NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
            return;

//        Vec3f acc = {0, 0, 0};
//        float h2_001 = 0.001f * pow(d_const->sph_h, 2);
        auto pos_i = DATA_VALUE(pos, p_i);
        auto vel_i = DATA_VALUE(vel, p_i);
        FOR_EACH_NEIGHBOR_Pj() {
            if (DATA_VALUE(mat, p_j) != DATA_VALUE(mat, p_i) || p_j == p_i)
                continue;

            auto pos_j = DATA_VALUE(pos, p_j);
            auto x_ij = pos_i = pos_j;
            auto vel_j = DATA_VALUE(vel, p_j);
            auto v_ij = vel_i - vel_j;
            auto wGrad = CUBIC_KERNEL_GRAD();
//            auto pi = -CONST_VALUE(rest_viscosity) * min(0.f, dot(vel_i - vel_j, pos_i - pos_j)) /
//                      ((pos_i - pos_j).length() * (pos_i - pos_j).length() + h2_001);

//            acc += -d_data->mass[p_i] * d_data->mass[p_j] * pi * CUBIC_KERNEL_GRAD();

//            if (DATA_VALUE(mat, p_j) == DATA_VALUE(mat, p_i)) {
            auto v_k1_mj = DATA_VALUE(vel_phase_1, p_i) - vel_j;
            DATA_VALUE(acc_phase_1, p_i) += 10 * DATA_VALUE(volume, p_j) *
                                            dot(CONST_VALUE(rest_viscosity) * (1 - CONST_VALUE(Cd) * v_k1_mj) +
                                                (CONST_VALUE(rest_viscosity) * CONST_VALUE(Cd) * v_ij),
                                                x_ij) * wGrad / dot(x_ij, x_ij);

            auto v_k2_mj = DATA_VALUE(vel_phase_2, p_i) - vel_j;
            DATA_VALUE(acc_phase_2, p_i) += 10 * DATA_VALUE(volume, p_j) *
                                            dot(CONST_VALUE(rest_viscosity) * (1 - CONST_VALUE(Cd) * v_k2_mj) +
                                                (CONST_VALUE(rest_viscosity) * CONST_VALUE(Cd) * v_ij),
                                                x_ij) * wGrad / dot(x_ij, x_ij);
//            }
        }

//        DATA_VALUE(acc_phase_1, p_i) += acc;
//        DATA_VALUE(acc_phase_2, p_i) += acc;
    }

    __global__ void
    add_phase_acc_surface_tension_cuda(IMSCTConstantParams *d_const,
                                       IMSCTDynamicParams *d_data,
                                       NeighborSearchUGConfig *d_nsConfig,
                                       NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
            return;

        Vec3f acc = {0, 0, 0};
        float gamma = 0.001;
        auto pos_i = DATA_VALUE(pos, p_i);
        FOR_EACH_NEIGHBOR_Pj() {
            if (DATA_VALUE(mat, p_j) != DATA_VALUE(mat, p_i))
                continue;

            auto pos_j = DATA_VALUE(pos, p_j);
            if ((pos_i - pos_j).length() < 1e-6)
                continue;


            acc += -gamma * DATA_VALUE(mass, p_i) * DATA_VALUE(mass, p_j) *
                   surface_tension_C((pos_i - pos_j).length(), CONST_VALUE(sph_h)) * (pos_i - pos_j) /
                   (pos_i - pos_j).length();
        }

        DATA_VALUE(acc_phase_1, p_i) += acc;
        DATA_VALUE(acc_phase_2, p_i) += acc;
    }

    __global__ void
    phase_acc_2_phase_vel_cuda(IMSCTConstantParams *d_const,
                               IMSCTDynamicParams *d_data,
                               NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
            return;

        DATA_VALUE(vel_phase_1, p_i) += DATA_VALUE(acc_phase_1, p_i) * CONST_VALUE(dt);
        DATA_VALUE(vel_phase_2, p_i) += DATA_VALUE(acc_phase_2, p_i) * CONST_VALUE(dt);
    }

    __global__ void
    update_vel_from_phase_vel_cuda(IMSCTConstantParams *d_const,
                                   IMSCTDynamicParams *d_data,
                                   NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
            return;

        DATA_VALUE(vel, p_i) = DATA_VALUE(vel_phase_1, p_i) * DATA_VALUE(vol_frac, p_i).x
                               + DATA_VALUE(vel_phase_2, p_i) * DATA_VALUE(vol_frac, p_i).y;
        DATA_VALUE(vel_adv, p_i) = DATA_VALUE(vel, p_i);

        DATA_VALUE(vel_drift_phase_1, p_i) = DATA_VALUE(vel_phase_1, p_i) - DATA_VALUE(vel, p_i);
        DATA_VALUE(vel_drift_phase_2, p_i) = DATA_VALUE(vel_phase_2, p_i) - DATA_VALUE(vel, p_i);
    }

    __global__ void
    get_acc_pressure_cuda(IMSCTConstantParams *d_const,
                          IMSCTDynamicParams *d_data,
                          NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
            return;

        DATA_VALUE(acc, p_i) = (DATA_VALUE(vel_adv, p_i) - DATA_VALUE(vel, p_i)) * CONST_VALUE(inv_dt);
    }

    __global__ void
    distribute_acc_pressure_2_phase_cuda(IMSCTConstantParams *d_const,
                                         IMSCTDynamicParams *d_data,
                                         NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
            return;

        DATA_VALUE(acc_phase_1, p_i) += DATA_VALUE(acc, p_i) * (CONST_VALUE(Cd) + (1 - CONST_VALUE(Cd)) *
                                                                                  (DATA_VALUE(rest_density, p_i) /
                                                                                   CONST_VALUE(rest_density).x));
        DATA_VALUE(acc_phase_2, p_i) += DATA_VALUE(acc, p_i) * (CONST_VALUE(Cd) + (1 - CONST_VALUE(Cd)) *
                                                                                  (DATA_VALUE(rest_density, p_i) /
                                                                                   CONST_VALUE(rest_density).y));
    }

    __global__ void
    compute_kappa_incomp_from_delta_compression_ratio_cuda(IMSCTConstantParams *d_const,
                                                           IMSCTDynamicParams *d_data,
                                                           NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
            return;

        DATA_VALUE(kappa_incomp, p_i) = DATA_VALUE(delta_compression_ratio, p_i) / DATA_VALUE(df_alpha, p_i) *
                                        CONST_VALUE(inv_dt2) / DATA_VALUE(volume, p_i);
        DATA_VALUE(df_alpha_2, p_i) += DATA_VALUE(kappa_incomp, p_i);
    }

    __global__ void
    vf_update_vel_adv_from_kappa_incomp_cuda(IMSCTConstantParams *d_const,
                                             IMSCTDynamicParams *d_data,
                                             NeighborSearchUGConfig *d_nsConfig,
                                             NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
            return;

        auto pos_i = DATA_VALUE(pos, p_i);
        FOR_EACH_NEIGHBOR_Pj() {
            if (p_j == p_i)
                continue;

            auto pos_j = DATA_VALUE(pos, p_j);
            auto wGrad = CUBIC_KERNEL_GRAD();

            DATA_VALUE(vel_adv, p_i) -= CONST_VALUE(dt) * DATA_VALUE(volume, p_i) * DATA_VALUE(volume, p_j) /
                                        DATA_VALUE(mass, p_i) *
                                        (DATA_VALUE(kappa_incomp, p_i) + DATA_VALUE(kappa_incomp, p_j)) * wGrad;
        }
    }

    __global__ void
    update_pos_cuda(IMSCTConstantParams *d_const,
                    IMSCTDynamicParams *d_data,
                    NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
            return;

        DATA_VALUE(pos, p_i) += DATA_VALUE(vel, p_i) * CONST_VALUE(dt);
        DATA_VALUE(pos_adv, p_i) = DATA_VALUE(pos, p_i);
    }

    __global__ void
    clear_val_frac_tmp_cuda(IMSCTConstantParams *d_const,
                            IMSCTDynamicParams *d_data,
                            NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
            return;

        DATA_VALUE(vol_frac_in, p_i) *= 0;
        DATA_VALUE(vol_frac_out, p_i) *= 0;
    }

    __global__ void
    update_phase_change_from_drift_cuda(IMSCTConstantParams *d_const,
                                        IMSCTDynamicParams *d_data,
                                        NeighborSearchUGConfig *d_nsConfig,
                                        NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
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
    update_phase_change_from_diffuse_cuda(IMSCTConstantParams *d_const,
                                          IMSCTDynamicParams *d_data,
                                          NeighborSearchUGConfig *d_nsConfig,
                                          NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
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
    check_negative_cuda(IMSCTConstantParams *d_const,
                        IMSCTDynamicParams *d_data,
                        NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (p_i == 0)
            g_all_positive = 1;
        __syncthreads();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
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
    update_phase_change_cuda(IMSCTConstantParams *d_const,
                             IMSCTDynamicParams *d_data,
                             NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
            return;

        DATA_VALUE(vol_frac, p_i) += (DATA_VALUE(vol_frac_in, p_i) + DATA_VALUE(vol_frac_out, p_i));
    }

    __global__ void
    release_unused_drift_vel_cuda(IMSCTConstantParams *d_const,
                                  IMSCTDynamicParams *d_data,
                                  NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
            return;

        if (DATA_VALUE(flag_negative_vol_frac, p_i) != 0) {
            DATA_VALUE(vel_phase_1, p_i) = DATA_VALUE(vel, p_i);
            DATA_VALUE(vel_phase_2, p_i) = DATA_VALUE(vel, p_i);
        }
    }

    __global__ void
    release_negative_cuda(IMSCTConstantParams *d_const,
                          IMSCTDynamicParams *d_data,
                          NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
            return;

        DATA_VALUE(flag_negative_vol_frac, p_i) = 0;
    }

    __global__ void
    regularize_val_frac_cuda(IMSCTConstantParams *d_const,
                             IMSCTDynamicParams *d_data,
                             NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
            return;

        float frac_sum = DATA_VALUE(vol_frac, p_i).x + DATA_VALUE(vol_frac, p_i).y;
        DATA_VALUE(vol_frac, p_i) /= frac_sum;
    }

//    __global__ void
//    compute_vel_grad_cuda(IMSCTConstantParams *d_const,
//                          IMSCTDynamicParams *d_data,
//                          NeighborSearchUGConfig *d_nsConfig,
//                          NeighborSearchUGParams *d_nsParams) {
//        CHECK_THREAD();
//
//        // applied to all dynamic objects
//        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
//            return;
//
//        DATA_VALUE(vel_grad, p_i) *= 0;
//        auto pos_i = DATA_VALUE(pos, p_i);
//        auto vel_i = DATA_VALUE(vel, p_i);
//
//        Mat33f A = {1, 2, 3,
//                    4, 5, 6, 7, 8, 9};
//        Mat33f B = {1, 2, 3,
//                    4, 5, 6, 7, 8, 9};
//
//        Vec3f C = {1, 1, 1};
//        Vec3f D = {4, 5, 6};
//
//        FOR_EACH_NEIGHBOR_Pj() {
//            if (DATA_VALUE(mat, p_j) != DATA_VALUE(mat, p_i))
//                continue;
//
//            auto pos_j = DATA_VALUE(pos, p_j);
//            auto vel_j = DATA_VALUE(vel, p_j);
//            auto vel_ji = vel_j - vel_i;
//
////            DATA_VALUE(vel_grad, p_i) += DATA_VALUE(compression_ratio, p_j) * vel_ji * CUBIC_KERNEL_GRAD();
//            DATA_VALUE(vel_grad, p_i) = vel_i * C;
//        }
//    }
//
//    __global__ void
//    update_conformation_tensor_cuda(IMSCTConstantParams *d_const,
//                                    IMSCTDynamicParams *d_data,
//                                    NeighborSearchUGParams *d_nsParams) {
//        CHECK_THREAD();
//
//        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
//            return;
//
//        auto CT_i = DATA_VALUE(CT, p_i);
//        auto d_CT = CT_i * DATA_VALUE(vel_grad, p_i) +
//                    DATA_VALUE(vel_grad, p_i).transpose() * CT_i -
//                    1.f / (CONST_VALUE(ct_relaxation_time) + 1e-5f) * (CT_i - Mat33f::eye());
//        DATA_VALUE(CT, p_i) += d_CT * CONST_VALUE(dt);
//    }
//
//    __global__ void
//    add_viscoelastic_acc_cuda(IMSCTConstantParams *d_const,
//                              IMSCTDynamicParams *d_data,
//                              NeighborSearchUGConfig *d_nsConfig,
//                              NeighborSearchUGParams *d_nsParams) {
//        CHECK_THREAD();
//
//        if (DATA_VALUE(mat, p_i) != IMSCT_NONNEWTON)
//            return;
//
//        auto CT_i = DATA_VALUE(CT, p_i);
//        DATA_VALUE(viscoelastic_stress, p_i) = DATA_VALUE(solution_vis, p_i) * (CT_i - Mat33f::eye())
//                                               - DATA_VALUE(ct_thinning_exp, p_i) * DATA_VALUE(solution_vis, p_i) *
//                                                 (CT_i * CT_i - CT_i);
//        __syncthreads();
//
//        DATA_VALUE(viscoelastic_stress, p_i) *= 0;
//        auto pos_i = DATA_VALUE(pos, p_i);
//        FOR_EACH_NEIGHBOR_Pj() {
//
//            if (DATA_VALUE(mat, p_j) != DATA_VALUE(mat, p_i))
//                continue;
//
//            auto pos_j = DATA_VALUE(pos, p_j);
//
//            DATA_VALUE(acc_phase_2, p_i) +=
//                    (DATA_VALUE(viscoelastic_stress, p_i) / powf(DATA_VALUE(rest_density, p_i), 2) +
//                     DATA_VALUE(viscoelastic_stress, p_j) / powf(DATA_VALUE(rest_density, p_j), 2)) *
//                    CUBIC_KERNEL_GRAD();
//        }
//    }
}


/**
 * host invoke impl
 */

namespace SoSim {
    __host__ void
    init_data(IMSCTConstantParams &h_const,
              IMSCTConstantParams *d_const,
              IMSCTDynamicParams *d_data,
              NeighborSearchUGParams *d_nsParams) {
        init_data_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    prepare_ims(IMSCTConstantParams &h_const,
                IMSCTConstantParams *d_const,
                IMSCTDynamicParams *d_data,
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
    sph_precompute(IMSCTConstantParams &h_const,
                   IMSCTConstantParams *d_const,
                   IMSCTDynamicParams *d_data,
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
    vfsph_div(IMSCTConstantParams &h_const,
              IMSCTDynamicParams &h_data,
              Vec3ui &obj_part_index,
              IMSCTConstantParams *d_const,
              IMSCTDynamicParams *d_data,
              NeighborSearchUGConfig *d_nsConfig,
              NeighborSearchUGParams *d_nsParams) {
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

            // compute_kappa_div_from_delta_compression_ratio()
            compute_kappa_div_from_delta_compression_ratio_cuda<<<h_const.block_num, h_const.thread_num>>>(
                    d_const, d_data, d_nsParams);

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

        // vel = vel_adv
        //        cudaMemcpy(h_data.vel, h_data.vel_adv, h_const.particle_num * sizeof(Vec3f), cudaMemcpyDeviceToDevice);
    }

    __host__ void
    apply_pressure_acc(IMSCTConstantParams &h_const,
                       IMSCTConstantParams *d_const,
                       IMSCTDynamicParams *d_data,
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
    ism_gravity_vis_surface(IMSCTConstantParams &h_const,
                            IMSCTConstantParams *d_const,
                            IMSCTDynamicParams *d_data,
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
    vfsph_incomp(IMSCTConstantParams &h_const,
                 IMSCTDynamicParams &h_data,
                 Vec3ui &obj_part_index,
                 IMSCTConstantParams *d_const,
                 IMSCTDynamicParams *d_data,
                 NeighborSearchUGConfig *d_nsConfig,
                 NeighborSearchUGParams *d_nsParams) {
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
    }

    __host__ void
    update_pos(IMSCTConstantParams &h_const,
               IMSCTConstantParams *d_const,
               IMSCTDynamicParams *d_data,
               NeighborSearchUGParams *d_nsParams) {
        update_pos_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    phase_transport_ism(IMSCTConstantParams &h_const,
                        IMSCTConstantParams *d_const,
                        IMSCTDynamicParams *d_data,
                        NeighborSearchUGConfig *d_nsConfig,
                        NeighborSearchUGParams *d_nsParams) {
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
    update_mass_and_vel(IMSCTConstantParams &h_const,
                        IMSCTConstantParams *d_const,
                        IMSCTDynamicParams *d_data,
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
    update_color(IMSCTConstantParams &h_const,
                 IMSCTConstantParams *d_const,
                 IMSCTDynamicParams *d_data,
                 NeighborSearchUGParams *d_nsParams) {
        // get_acc_pressure()
        update_color_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

//    __host__ void
//    ism_viscoelastic(IMSCTConstantParams &h_const,
//                     IMSCTConstantParams *d_const,
//                     IMSCTDynamicParams *d_data,
//                     NeighborSearchUGConfig *d_nsConfig,
//                     NeighborSearchUGParams *d_nsParams) {
//        // compute_vel_grad()
//        compute_vel_grad_cuda<<<h_const.block_num, h_const.thread_num>>>(
//                d_const, d_data, d_nsConfig, d_nsParams);
//
//        // update_conformation_tensor()
//        update_conformation_tensor_cuda<<<h_const.block_num, h_const.thread_num>>>(
//                d_const, d_data, d_nsParams);
//
//        // clear_phase_acc()
//        clear_phase_acc_cuda<<<h_const.block_num, h_const.thread_num>>>(
//                d_const, d_data, d_nsParams);
//
//        add_viscoelastic_acc_cuda<<<h_const.block_num, h_const.thread_num>>>(
//                d_const, d_data, d_nsConfig, d_nsParams);
//
//        // phase_acc_2_phase_vel()
//        phase_acc_2_phase_vel_cuda<<<h_const.block_num, h_const.thread_num>>>(
//                d_const, d_data, d_nsParams);
//
//        // update_vel_from_phase_vel()
//        update_vel_from_phase_vel_cuda<<<h_const.block_num, h_const.thread_num>>>(
//                d_const, d_data, d_nsParams);
//    }
}