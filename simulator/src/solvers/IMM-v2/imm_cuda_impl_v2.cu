//
// Created by ADMIN on 2024/3/26.
//
//

#include "imm_cuda_api_v2.cuh"
#include "imm_macro_v2.hpp"
#include "libs/SPHKernelL/kernels.cuh"
#include "libs/AnalysisL/statistic_util.hpp"

/**
 * cuda impl
 */

namespace SoSim {
    __global__ void
    update_const_vec_ptr_cuda(IMMConstantParams_v2 *d_const,
                              float *d_phase_density,
                              Vec3f *d_phase_color,
                              float *d_phase_vis) {
        d_const->phase_rest_density = d_phase_density;
        d_const->phase_color = d_phase_color;
        d_const->phase_vis = d_phase_vis;
    }

    __global__ void
    init_data_cuda(IMMConstantParams_v2 *d_const,
                   IMMDynamicParams_v2 *d_data) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->particle_num)
            return;

        DATA_VALUE(volume, i) = CONST_VALUE(rest_volume);
        DATA_VALUE(kappa_div, i) = 0;
        DATA_VALUE(acc, i) *= 0;
        DATA_VALUE(vel_adv, i) = DATA_VALUE(vel, i);
        DATA_VALUE(rest_density, i) = 0;
        DATA_VALUE(flag_negative_vol_frac, i) = 0;

        FOR_EACH_PHASE_k() {
            DATA_VALUE_PHASE(vol_frac_in, i, k) = 0;
            DATA_VALUE_PHASE(vol_frac_out, i, k) = 0;
            DATA_VALUE_PHASE(vel_phase, i, k) = DATA_VALUE(vel, i);
            DATA_VALUE_PHASE(acc_phase, i, k) *= 0;
            DATA_VALUE_PHASE(vel_drift_phase, i, k) *= 0;
            DATA_VALUE(rest_density, i) += DATA_VALUE_PHASE(vol_frac, i, k) * CONST_VALUE_PHASE(phase_rest_density, k);
        }
    }

    __global__ void
    update_rest_density_and_mass_cuda(IMMConstantParams_v2 *d_const,
                                      IMMDynamicParams_v2 *d_data,
                                      NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON && DATA_VALUE(mat, p_i) != Emitter_Particle)
            return;

        DATA_VALUE(rest_density, p_i) *= 0;
        FOR_EACH_PHASE_k() {
            DATA_VALUE(rest_density, p_i) += DATA_VALUE_PHASE(vol_frac, p_i, k) *
                                             CONST_VALUE_PHASE(phase_rest_density, k);
        }
        DATA_VALUE(mass, p_i) = DATA_VALUE(rest_density, p_i) * DATA_VALUE(volume, p_i);
    }

    __global__ void
    compute_rigid_volume(IMMConstantParams_v2 *d_const,
                         IMMDynamicParams_v2 *d_data,
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
    compute_compression_ratio_cuda(IMMConstantParams_v2 *d_const,
                                   IMMDynamicParams_v2 *d_data,
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
    compute_df_beta_cuda(IMMConstantParams_v2 *d_const,
                         IMMDynamicParams_v2 *d_data,
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
    compute_delta_compression_ratio_cuda(IMMConstantParams_v2 *d_const,
                                         IMMDynamicParams_v2 *d_data,
                                         NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(delta_compression_ratio, p_i) = DATA_VALUE(compression_ratio, p_i) - 1.f;
    }

    __global__ void
    update_delta_compression_ratio_from_vel_adv_cuda(IMMConstantParams_v2 *d_const,
                                                     IMMDynamicParams_v2 *d_data,
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
    compute_kappa_div_from_delta_compression_ratio_cuda(IMMConstantParams_v2 *d_const,
                                                        IMMDynamicParams_v2 *d_data,
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
    vf_update_vel_adv_from_kappa_div_cuda(IMMConstantParams_v2 *d_const,
                                          IMMDynamicParams_v2 *d_data,
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
    clear_acc_cuda(IMMConstantParams_v2 *d_const,
                   IMMDynamicParams_v2 *d_data,
                   NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        DATA_VALUE(acc, p_i) *= 0;
    }

    __global__ void
    get_acc_pressure_cuda(IMMConstantParams_v2 *d_const,
                          IMMDynamicParams_v2 *d_data,
                          NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        DATA_VALUE(acc, p_i) = (DATA_VALUE(vel_adv, p_i) - DATA_VALUE(vel, p_i)) * CONST_VALUE(inv_dt);
    }

    __global__ void
    clear_phase_acc_cuda(IMMConstantParams_v2 *d_const,
                         IMMDynamicParams_v2 *d_data,
                         NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        FOR_EACH_PHASE_k() {
            DATA_VALUE_PHASE(acc_phase, p_i, k) *= 0;
        }
    }

    __global__ void
    add_acc_gravity_cuda(IMMConstantParams_v2 *d_const,
                         IMMDynamicParams_v2 *d_data,
                         NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(acc, p_i) += CONST_VALUE(gravity);
    }

    __global__ void
    add_phase_acc_gravity_cuda(IMMConstantParams_v2 *d_const,
                               IMMDynamicParams_v2 *d_data,
                               NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        FOR_EACH_PHASE_k() {
            DATA_VALUE_PHASE(acc_phase, p_i, k) += CONST_VALUE(gravity);
        }
    }

    __global__ void
    add_phase_acc_vis_cuda(IMMConstantParams_v2 *d_const,
                           IMMDynamicParams_v2 *d_data,
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

            FOR_EACH_PHASE_k() {
                auto v_k_mj = DATA_VALUE_PHASE(vel_phase, p_i, k) - vel_j;
                DATA_VALUE_PHASE(acc_phase, p_i, k) += 10 * DATA_VALUE(volume, p_j) *
                                                       dot(CONST_VALUE_PHASE(phase_vis, k) * (1 - CONST_VALUE(Cd)) *
                                                           v_k_mj +
                                                           (CONST_VALUE_PHASE(phase_vis, k) * CONST_VALUE(Cd) * v_ij),
                                                           x_ij) * wGrad / dot(x_ij, x_ij);
            }
        }
    }

    __global__ void
    compute_surface_normal_cuda(IMMConstantParams_v2 *d_const,
                                IMMDynamicParams_v2 *d_data,
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
    add_acc_surface_tension_cuda(IMMConstantParams_v2 *d_const,
                                 IMMDynamicParams_v2 *d_data,
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

        DATA_VALUE(acc, p_i) += acc;
    }

    __global__ void
    add_phase_acc_surface_tension_cuda(IMMConstantParams_v2 *d_const,
                                       IMMDynamicParams_v2 *d_data,
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

        FOR_EACH_PHASE_k() {
            DATA_VALUE_PHASE(acc_phase, p_i, k) += acc;
        }
    }

    __global__ void
    phase_acc_2_phase_vel_cuda(IMMConstantParams_v2 *d_const,
                               IMMDynamicParams_v2 *d_data,
                               NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        FOR_EACH_PHASE_k() {
            DATA_VALUE_PHASE(vel_phase, p_i, k) += DATA_VALUE_PHASE(acc_phase, p_i, k) * CONST_VALUE(dt);
        }
    }

    __global__ void
    update_vel_from_phase_vel_cuda(IMMConstantParams_v2 *d_const,
                                   IMMDynamicParams_v2 *d_data,
                                   NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(vel, p_i) *= 0;
        FOR_EACH_PHASE_k() {
            DATA_VALUE(vel, p_i) += DATA_VALUE_PHASE(vel_phase, p_i, k) * DATA_VALUE_PHASE(vol_frac, p_i, k);
        }

        DATA_VALUE(vel_adv, p_i) = DATA_VALUE(vel, p_i);

        FOR_EACH_PHASE_k() {
            DATA_VALUE_PHASE(vel_drift_phase, p_i, k) =
                    DATA_VALUE_PHASE(vel_phase, p_i, k) - DATA_VALUE(vel, p_i);
        }
    }

    __global__ void
    acc_2_vel(IMMConstantParams_v2 *d_const,
              IMMDynamicParams_v2 *d_data,
              NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(vel, p_i) += DATA_VALUE(acc, p_i) * CONST_VALUE(dt);
        DATA_VALUE(vel_adv, p_i) = DATA_VALUE(vel, p_i);
    }

    __global__ void
    distribute_acc_pressure_2_phase_cuda(IMMConstantParams_v2 *d_const,
                                         IMMDynamicParams_v2 *d_data,
                                         NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        FOR_EACH_PHASE_k() {
            DATA_VALUE_PHASE(acc_phase, p_i, k) +=
                    DATA_VALUE(acc, p_i) * (CONST_VALUE(Cd) + (1 - CONST_VALUE(Cd)) *
                                                              (DATA_VALUE(rest_density, p_i) /
                                                               CONST_VALUE_PHASE(phase_rest_density, k)));
        }
    }

    __global__ void
    compute_kappa_incomp_from_delta_compression_ratio_cuda(IMMConstantParams_v2 *d_const,
                                                           IMMDynamicParams_v2 *d_data,
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
    vf_update_vel_adv_from_kappa_incomp_cuda(IMMConstantParams_v2 *d_const,
                                             IMMDynamicParams_v2 *d_data,
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
    update_pos_cuda(IMMConstantParams_v2 *d_const,
                    IMMDynamicParams_v2 *d_data,
                    NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(pos, p_i) += DATA_VALUE(vel, p_i) * CONST_VALUE(dt);
        DATA_VALUE(pos_adv, p_i) = DATA_VALUE(pos, p_i);
    }

    __global__ void
    regularize_vol_frac_cuda(IMMConstantParams_v2 *d_const,
                             IMMDynamicParams_v2 *d_data,
                             NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        float frac_sum = 0;
        FOR_EACH_PHASE_k() {
            frac_sum += DATA_VALUE_PHASE(vol_frac, p_i, k);
        }

        FOR_EACH_PHASE_k() {
            DATA_VALUE_PHASE(vol_frac, p_i, k) /= frac_sum;
        }
    }

    __global__ void
    update_color_cuda(IMMConstantParams_v2 *d_const,
                      IMMDynamicParams_v2 *d_data,
                      NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(color, p_i) *= 0;
        FOR_EACH_PHASE_k() {
            DATA_VALUE(color, p_i) += DATA_VALUE_PHASE(vol_frac, p_i, k) * CONST_VALUE_PHASE(phase_color, k);
        }
    }

    __global__ void
    estimate_density_cuda(IMMConstantParams_v2 *d_const,
                          IMMDynamicParams_v2 *d_data,
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
    clear_vol_frac_tmp_cuda(IMMConstantParams_v2 *d_const,
                            IMMDynamicParams_v2 *d_data,
                            NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        FOR_EACH_PHASE_k() {
            DATA_VALUE_PHASE(vol_frac_in, p_i, k) = 0;
            DATA_VALUE_PHASE(vol_frac_out, p_i, k) = 0;
        }
    }

    __global__ void
    update_phase_change_from_drift_cuda(IMMConstantParams_v2 *d_const,
                                        IMMDynamicParams_v2 *d_data,
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

            FOR_EACH_PHASE_k() {
                float vol_frac_change =
                        -CONST_VALUE(dt) * DATA_VALUE(volume, p_j) * dot(DATA_VALUE_PHASE(vol_frac, p_i, k) *
                                                                         DATA_VALUE_PHASE(vel_drift_phase, p_i, k) +
                                                                         DATA_VALUE_PHASE(vol_frac, p_j, k) *
                                                                         DATA_VALUE_PHASE(vel_drift_phase, p_j, k),
                                                                         wGrad);

                if (vol_frac_change < 0)
                    DATA_VALUE_PHASE(vol_frac_out, p_i, k) += vol_frac_change;
                else
                    DATA_VALUE_PHASE(vol_frac_in, p_i, k) += vol_frac_change;
            }
        }
    }

    __global__ void
    update_phase_change_from_diffuse_cuda(IMMConstantParams_v2 *d_const,
                                          IMMDynamicParams_v2 *d_data,
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

            FOR_EACH_PHASE_k() {
                float vol_frac_ij = DATA_VALUE_PHASE(vol_frac, p_i, k) - DATA_VALUE_PHASE(vol_frac, p_j, k);
                float vol_frac_change = CONST_VALUE(dt) * CONST_VALUE(Cf) * vol_frac_ij * DATA_VALUE(volume, p_j) *
                                        factor;

                if (vol_frac_change < 0)
                    DATA_VALUE_PHASE(vol_frac_out, p_i, k) += vol_frac_change;
                else
                    DATA_VALUE_PHASE(vol_frac_in, p_i, k) += vol_frac_change;
            }
        }
    }

    __device__ float g_all_positive;

    __global__ void
    check_negative_cuda(IMMConstantParams_v2 *d_const,
                        IMMDynamicParams_v2 *d_data,
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

        FOR_EACH_PHASE_k() {
            auto vol_frac_tmp = DATA_VALUE_PHASE(vol_frac, p_i, k) + DATA_VALUE_PHASE(vol_frac_out, p_i, k)
                                + DATA_VALUE_PHASE(vol_frac_in, p_i, k);
            if (vol_frac_tmp < 0) {
                DATA_VALUE(flag_negative_vol_frac, p_i) = 1;
                atomicAdd(&g_all_positive, -1);
            }
        }
    }

    __global__ void
    update_phase_change_cuda(IMMConstantParams_v2 *d_const,
                             IMMDynamicParams_v2 *d_data,
                             NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        FOR_EACH_PHASE_k() {
            DATA_VALUE_PHASE(vol_frac, p_i, k) += DATA_VALUE_PHASE(vol_frac_out, p_i, k)
                                                  + DATA_VALUE_PHASE(vol_frac_in, p_i, k);
        }
    }

    __global__ void
    release_unused_drift_vel_cuda(IMMConstantParams_v2 *d_const,
                                  IMMDynamicParams_v2 *d_data,
                                  NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        if (DATA_VALUE(flag_negative_vol_frac, p_i) != 0) {
            FOR_EACH_PHASE_k() {
                DATA_VALUE_PHASE(vel_phase, p_i, k) = DATA_VALUE(vel, p_i);
            }
        }
    }

    __global__ void
    release_negative_cuda(IMMConstantParams_v2 *d_const,
                          IMMDynamicParams_v2 *d_data,
                          NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(flag_negative_vol_frac, p_i) = 0;
    }
}

namespace SoSim { // extra func cuda impl
    __global__ void
    stir_fan_cuda(IMMConstantParams_v2 *d_const,
                  IMMDynamicParams_v2 *d_data,
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
    buckling_move_cuda(IMMConstantParams_v2 *d_const,
                       IMMDynamicParams_v2 *d_data,
                       NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != MOVING_TUBE && DATA_VALUE(mat, p_i) != MOVING_COVER)
            return;

        DATA_VALUE(pos, p_i) += CONST_VALUE(dt) * DATA_VALUE(vel, p_i);
    }

    __global__ void
    rotate_bowl_cuda(IMMConstantParams_v2 *d_const,
                     IMMDynamicParams_v2 *d_data,
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
}


/**
 * host invoke impl
 */

namespace SoSim {
    __host__ void
    init_data(IMMConstantParams_v2 &h_const,
              IMMConstantParams_v2 *d_const,
              IMMDynamicParams_v2 *d_data,
              float *d_phase_density,
              Vec3f *d_phase_color,
              float *d_phase_vis) {
        update_const_vec_ptr_cuda<<<1, 1>>>(
                d_const, d_phase_density, d_phase_color, d_phase_vis);

        init_data_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data);
    }

    __host__ void
    prepare_imm(IMMConstantParams_v2 &h_const,
                IMMConstantParams_v2 *d_const,
                IMMDynamicParams_v2 *d_data,
                NeighborSearchUGConfig *d_nsConfig,
                NeighborSearchUGParams *d_nsParams) {

        update_rest_density_and_mass_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        compute_rigid_volume<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);
    }

    __host__ void
    sph_precompute(IMMConstantParams_v2 &h_const,
                   IMMConstantParams_v2 *d_const,
                   IMMDynamicParams_v2 *d_data,
                   NeighborSearchUGConfig *d_nsConfig,
                   NeighborSearchUGParams *d_nsParams) {
        compute_compression_ratio_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);

        compute_df_beta_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);
    }

    __host__ void
    vfsph_div(IMMConstantParams_v2 &h_const,
              IMMDynamicParams_v2 &h_data,
              Vec3ui &obj_part_index,
              IMMConstantParams_v2 *d_const,
              IMMDynamicParams_v2 *d_data,
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
    apply_pressure_acc(IMMConstantParams_v2 &h_const,
                       IMMConstantParams_v2 *d_const,
                       IMMDynamicParams_v2 *d_data,
                       NeighborSearchUGParams *d_nsParams) {
        get_acc_pressure_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        clear_phase_acc_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        distribute_acc_pressure_2_phase_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        phase_acc_2_phase_vel_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        update_vel_from_phase_vel_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    imm_gravity_vis_surface(IMMConstantParams_v2 &h_const,
                            IMMConstantParams_v2 *d_const,
                            IMMDynamicParams_v2 *d_data,
                            NeighborSearchUGConfig *d_nsConfig,
                            NeighborSearchUGParams *d_nsParams) {
        clear_phase_acc_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        add_phase_acc_gravity_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        add_phase_acc_vis_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);

        compute_surface_normal_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);

        add_phase_acc_surface_tension_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);

        phase_acc_2_phase_vel_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        update_vel_from_phase_vel_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    vfsph_incomp(IMMConstantParams_v2 &h_const,
                 IMMDynamicParams_v2 &h_data,
                 Vec3ui &obj_part_index,
                 IMMConstantParams_v2 *d_const,
                 IMMDynamicParams_v2 *d_data,
                 NeighborSearchUGConfig *d_nsConfig,
                 NeighborSearchUGParams *d_nsParams,
                 bool &crash) {
        int iter = 0;
        while (true) {
            iter++;

            compute_delta_compression_ratio_cuda<<<h_const.block_num, h_const.thread_num>>>(
                    d_const, d_data, d_nsParams);

            update_delta_compression_ratio_from_vel_adv_cuda<<<h_const.block_num, h_const.thread_num>>>(
                    d_const, d_data, d_nsConfig, d_nsParams);

            auto compressible_ratio = cal_mean(h_data.delta_compression_ratio,
                                               h_const.particle_num, obj_part_index.y);

            compute_kappa_incomp_from_delta_compression_ratio_cuda<<<h_const.block_num, h_const.thread_num>>>(
                    d_const, d_data, d_nsParams);

            vf_update_vel_adv_from_kappa_incomp_cuda<<<h_const.block_num, h_const.thread_num>>>(
                    d_const, d_data, d_nsConfig, d_nsParams);

            if (compressible_ratio < h_const.incompressible_threshold || iter > 100)
                break;
        }

        std::cout << "incomp-iter: " << iter << '\n';

        if (iter == 101)
            crash = true;
    }

    __host__ void
    update_pos(IMMConstantParams_v2 &h_const,
               IMMConstantParams_v2 *d_const,
               IMMDynamicParams_v2 *d_data,
               NeighborSearchUGParams *d_nsParams) {
        update_pos_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    phase_transport_ism(IMMConstantParams_v2 &h_const,
                        IMMConstantParams_v2 *d_const,
                        IMMDynamicParams_v2 *d_data,
                        NeighborSearchUGConfig *d_nsConfig,
                        NeighborSearchUGParams *d_nsParams,
                        bool &crash) {
        // clear_val_frac_tmp()
        clear_vol_frac_tmp_cuda<<<h_const.block_num, h_const.thread_num>>>(
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
            clear_vol_frac_tmp_cuda<<<h_const.block_num, h_const.thread_num>>>(
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
        regularize_vol_frac_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // update_rest_density_and_mass()
        update_rest_density_and_mass_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // update_vel_from_phase_vel()
        update_vel_from_phase_vel_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    update_mass_and_vel(IMMConstantParams_v2 &h_const,
                        IMMConstantParams_v2 *d_const,
                        IMMDynamicParams_v2 *d_data,
                        NeighborSearchUGParams *d_nsParams) {
        regularize_vol_frac_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        update_rest_density_and_mass_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        update_vel_from_phase_vel_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    update_color(IMMConstantParams_v2 &h_const,
                 IMMConstantParams_v2 *d_const,
                 IMMDynamicParams_v2 *d_data,
                 NeighborSearchUGParams *d_nsParams) {
        update_color_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }
}

namespace SoSim { // extra func host invoke
    __host__ void
    stirring(IMMConstantParams_v2 &h_const,
             IMMConstantParams_v2 *d_const,
             IMMDynamicParams_v2 *d_data,
             NeighborSearchUGParams *d_nsParams) {
        stir_fan_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    rotate_bowl(IMMConstantParams_v2 &h_const,
                IMMConstantParams_v2 *d_const,
                IMMDynamicParams_v2 *d_data,
                NeighborSearchUGParams *d_nsParams) {
        rotate_bowl_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    buckling(IMMConstantParams_v2 &h_const,
             IMMConstantParams_v2 *d_const,
             IMMDynamicParams_v2 *d_data,
             NeighborSearchUGParams *d_nsParams) {
        buckling_move_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }
}