//
// Created by ADMIN on 2024/3/26.
//
//

#include "dfsph_cuda_api.cuh"
#include "dfsph_macro.hpp"
#include "libs/SPHKernelL/kernels.cuh"
#include "libs/AnalysisL/statistic_util.hpp"

/**
 * cuda impl
 */

namespace SoSim {
    __global__ void
    init_data_cuda(DFSPHConstantParams *d_const,
                   DFSPHDynamicParams *d_data,
                   NeighborSearchUGParams *d_nsParams) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->particle_num)
            return;

        DATA_VALUE(volume, i) = CONST_VALUE(rest_volume);
        DATA_VALUE(kappa_div, i) = 0;
        DATA_VALUE(acc, i) *= 0;
        DATA_VALUE(vel_adv, i) = DATA_VALUE(vel, i);
        DATA_VALUE(vis, i) = CONST_VALUE(rest_viscosity);
        DATA_VALUE(rest_density, i) = CONST_VALUE(rest_density);
    }

    __global__ void
    update_mass_cuda(DFSPHConstantParams *d_const,
                     DFSPHDynamicParams *d_data,
                     NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON && DATA_VALUE(mat, p_i) != Emitter_Particle)
            return;

        DATA_VALUE(mass, p_i) = DATA_VALUE(rest_density, p_i) * DATA_VALUE(volume, p_i);
    }

    __global__ void
    compute_rigid_volume(DFSPHConstantParams *d_const,
                         DFSPHDynamicParams *d_data,
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
    compute_compression_ratio_cuda(DFSPHConstantParams *d_const,
                                   DFSPHDynamicParams *d_data,
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
    compute_df_beta_cuda(DFSPHConstantParams *d_const,
                         DFSPHDynamicParams *d_data,
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
    compute_delta_compression_ratio_cuda(DFSPHConstantParams *d_const,
                                         DFSPHDynamicParams *d_data,
                                         NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(delta_compression_ratio, p_i) = DATA_VALUE(compression_ratio, p_i) - 1.f;
    }

    __global__ void
    update_delta_compression_ratio_from_vel_adv_cuda(DFSPHConstantParams *d_const,
                                                     DFSPHDynamicParams *d_data,
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
    compute_kappa_div_from_delta_compression_ratio_cuda(DFSPHConstantParams *d_const,
                                                        DFSPHDynamicParams *d_data,
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
    vf_update_vel_adv_from_kappa_div_cuda(DFSPHConstantParams *d_const,
                                          DFSPHDynamicParams *d_data,
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
    clear_acc_cuda(DFSPHConstantParams *d_const,
                   DFSPHDynamicParams *d_data,
                   NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        DATA_VALUE(acc, p_i) *= 0;
    }

    __global__ void
    add_acc_gravity_cuda(DFSPHConstantParams *d_const,
                         DFSPHDynamicParams *d_data,
                         NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(acc, p_i) += CONST_VALUE(gravity);
    }

    __global__ void
    add_acc_explicit_lap_vis_cuda(DFSPHConstantParams *d_const,
                                  DFSPHDynamicParams *d_data,
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
        FOR_EACH_NEIGHBOR_Pj() {
            if (DATA_VALUE(mat, p_j) != DATA_VALUE(mat, p_i))
                continue;

            auto pos_j = DATA_VALUE(pos, p_j);
            auto x_ij = pos_i - pos_j;
            auto vel_j = DATA_VALUE(vel, p_j);
            auto v_ij = vel_i - vel_j;

            auto vis = (DATA_VALUE(vis, p_i) + DATA_VALUE(vis, p_j)) / 2;

            auto pi = vis * DATA_VALUE(mass, p_j) / DATA_VALUE(rest_density, p_j) * dot(v_ij, x_ij) /
                      (x_ij.length() * x_ij.length() + h2_001);

            acc += 10 * pi * CUBIC_KERNEL_GRAD();
        }

        DATA_VALUE(acc, p_i) += acc;
    }

    __global__ void
    compute_surface_normal_cuda(DFSPHConstantParams *d_const,
                                DFSPHDynamicParams *d_data,
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
    add_acc_surface_tension_cuda(DFSPHConstantParams *d_const,
                                 DFSPHDynamicParams *d_data,
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
    acc_2_vel(DFSPHConstantParams *d_const,
              DFSPHDynamicParams *d_data,
              NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(vel, p_i) += DATA_VALUE(acc, p_i) * CONST_VALUE(dt);
        DATA_VALUE(vel_adv, p_i) = DATA_VALUE(vel, p_i);
    }

    __global__ void
    add_acc_pressure_cuda(DFSPHConstantParams *d_const,
                          DFSPHDynamicParams *d_data,
                          NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(acc, p_i) += (DATA_VALUE(vel_adv, p_i) - DATA_VALUE(vel, p_i)) * CONST_VALUE(inv_dt);
    }

    __global__ void
    compute_kappa_incomp_from_delta_compression_ratio_cuda(DFSPHConstantParams *d_const,
                                                           DFSPHDynamicParams *d_data,
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
    vf_update_vel_adv_from_kappa_incomp_cuda(DFSPHConstantParams *d_const,
                                             DFSPHDynamicParams *d_data,
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
    update_pos_cuda(DFSPHConstantParams *d_const,
                    DFSPHDynamicParams *d_data,
                    NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        // applied to all dynamic objects
        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(pos, p_i) += DATA_VALUE(vel, p_i) * CONST_VALUE(dt);
        DATA_VALUE(pos_adv, p_i) = DATA_VALUE(pos, p_i);
    }

    __global__ void
    estimate_density_cuda(DFSPHConstantParams *d_const,
                          DFSPHDynamicParams *d_data,
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
}

namespace SoSim { // extra func cuda impl
    __global__ void
    stir_fan_cuda(DFSPHConstantParams *d_const,
                  DFSPHDynamicParams *d_data,
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
    buckling_move_cuda(DFSPHConstantParams *d_const,
                       DFSPHDynamicParams *d_data,
                       NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != MOVING_TUBE && DATA_VALUE(mat, p_i) != MOVING_COVER)
            return;

        DATA_VALUE(pos, p_i) += CONST_VALUE(dt) * DATA_VALUE(vel, p_i);
    }

    __global__ void
    rotate_bowl_cuda(DFSPHConstantParams *d_const,
                     DFSPHDynamicParams *d_data,
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
    correct_vel_by_artificial_vis_bound_cuda(DFSPHConstantParams *d_const,
                                             DFSPHDynamicParams *d_data,
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

        DATA_VALUE(vel, p_i) *= f1;
    }
}


/**
 * host invoke impl
 */

namespace SoSim {
    __host__ void
    init_data(DFSPHConstantParams &h_const,
              DFSPHConstantParams *d_const,
              DFSPHDynamicParams *d_data,
              NeighborSearchUGParams *d_nsParams) {
        init_data_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    prepare_dfsph(DFSPHConstantParams &h_const,
                  DFSPHConstantParams *d_const,
                  DFSPHDynamicParams *d_data,
                  NeighborSearchUGConfig *d_nsConfig,
                  NeighborSearchUGParams *d_nsParams) {
        // ims update_rest_density_and_mass()
        update_mass_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // compute_rigid_volume()
        compute_rigid_volume<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);
    }

    __host__ void
    sph_precompute(DFSPHConstantParams &h_const,
                   DFSPHConstantParams *d_const,
                   DFSPHDynamicParams *d_data,
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
    vfsph_div(DFSPHConstantParams &h_const,
              DFSPHDynamicParams &h_data,
              Vec3ui &obj_part_index,
              DFSPHConstantParams *d_const,
              DFSPHDynamicParams *d_data,
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
    apply_pressure_acc(DFSPHConstantParams &h_const,
                       DFSPHConstantParams *d_const,
                       DFSPHDynamicParams *d_data,
                       NeighborSearchUGParams *d_nsParams) {
        // clear_phase_acc()
        clear_acc_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // get_acc_pressure()
        add_acc_pressure_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // phase_acc_2_phase_vel()
        acc_2_vel<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    dfsph_gravity_vis_surface(DFSPHConstantParams &h_const,
                              DFSPHConstantParams *d_const,
                              DFSPHDynamicParams *d_data,
                              NeighborSearchUGConfig *d_nsConfig,
                              NeighborSearchUGParams *d_nsParams) {
        // clear_phase_acc()
        clear_acc_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // add_phase_acc_gravity()
        add_acc_gravity_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);

        // add_acc_explicit_lap_vis_cuda()
        add_acc_explicit_lap_vis_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);

        // compute_surface_normal
        compute_surface_normal_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);

        // add_acc_surface_tension_cuda()
        add_acc_surface_tension_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);

        // acc_2_vel()
        acc_2_vel<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    vfsph_incomp(DFSPHConstantParams &h_const,
                 DFSPHDynamicParams &h_data,
                 Vec3ui &obj_part_index,
                 DFSPHConstantParams *d_const,
                 DFSPHDynamicParams *d_data,
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
    update_pos(DFSPHConstantParams &h_const,
               DFSPHConstantParams *d_const,
               DFSPHDynamicParams *d_data,
               NeighborSearchUGParams *d_nsParams) {
        update_pos_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    artificial_vis_bound(DFSPHConstantParams &h_const,
                         DFSPHConstantParams *d_const,
                         DFSPHDynamicParams *d_data,
                         NeighborSearchUGConfig *d_nsConfig,
                         NeighborSearchUGParams *d_nsParams) {
        // correct_phase_vel_by_artificial_vis_bound()
        correct_vel_by_artificial_vis_bound_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsConfig, d_nsParams);
    }
}

namespace SoSim { // extra func host invoke
    __host__ void
    stirring(DFSPHConstantParams &h_const,
             DFSPHConstantParams *d_const,
             DFSPHDynamicParams *d_data,
             NeighborSearchUGParams *d_nsParams) {
        stir_fan_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    rotate_bowl(DFSPHConstantParams &h_const,
                DFSPHConstantParams *d_const,
                DFSPHDynamicParams *d_data,
                NeighborSearchUGParams *d_nsParams) {
        rotate_bowl_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }

    __host__ void
    buckling(DFSPHConstantParams &h_const,
             DFSPHConstantParams *d_const,
             DFSPHDynamicParams *d_data,
             NeighborSearchUGParams *d_nsParams) {
        buckling_move_cuda<<<h_const.block_num, h_const.thread_num>>>(
                d_const, d_data, d_nsParams);
    }
}