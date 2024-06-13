//
// Created by ADMIN on 2024/6/13.
//

#include "pbf_cuda_api.cuh"

#include "pbf_macro.hpp"
#include "libs/SPHKernelL/kernels.cuh"

namespace SoSim { // cuda kernels
    __global__ void
    init_cuda(PBFConstantParams *d_const,
              PBFDynamicParams *d_data) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->particle_num)
            return;

        DATA_VALUE(volume, i) = CONST_VALUE(rest_volume);
        DATA_VALUE(dx, i) *= 0;
        DATA_VALUE(mass, i) = CONST_VALUE(rest_density) * CONST_VALUE(rest_volume);
        DATA_VALUE(error, i) = 0;
        DATA_VALUE(error_grad, i) *= 0;
    }

    __global__ void
    compute_rigid_particle_volume_cuda(PBFConstantParams *d_const,
                                       PBFDynamicParams *d_data,
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

        if (DATA_VALUE(mat, p_i) == FIXED_BOUND)
            DATA_VALUE(mass, p_i) = CONST_VALUE(rest_bound_density) * DATA_VALUE(volume, p_i);
        else if (DATA_VALUE(mat, p_i) == DYNAMIC_RIGID)
            DATA_VALUE(mass, p_i) = CONST_VALUE(rest_rigid_density) * DATA_VALUE(volume, p_i);
    }

    __global__ void
    compute_sph_density_and_error_cuda(PBFConstantParams *d_const,
                                       PBFDynamicParams *d_data,
                                       NeighborSearchUGConfig *d_nsConfig,
                                       NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        auto pos_i = DATA_VALUE(pos, p_i);
        DATA_VALUE(density_sph, p_i) *= 0;

        FOR_EACH_NEIGHBOR_Pj() {
            auto pos_j = DATA_VALUE(pos, p_j);
            auto m_j = DATA_VALUE(mass, p_j);

            DATA_VALUE(density_sph, p_i) += m_j * CUBIC_KERNEL_VALUE();
        }

        if (DATA_VALUE(density_sph, p_i) < CONST_VALUE(rest_density))
            DATA_VALUE(density_sph, p_i) = CONST_VALUE(rest_density);

        DATA_VALUE(error, p_i) = DATA_VALUE(density_sph, p_i) / CONST_VALUE(rest_density) - 1.f;
    }

    __global__ void
    update_lamb_cuda(PBFConstantParams *d_const,
                     PBFDynamicParams *d_data,
                     NeighborSearchUGConfig *d_nsConfig,
                     NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        auto pos_i = DATA_VALUE(pos, p_i);
        DATA_VALUE(error_grad, p_i) *= 0;
        DATA_VALUE(lamb, p_i) = 0;

        FOR_EACH_NEIGHBOR_Pj() {
            auto pos_j = DATA_VALUE(pos, p_j);
            auto wGrad = CUBIC_KERNEL_GRAD();

            auto error_grad_j = -1 / CONST_VALUE(rest_density) * wGrad;
            DATA_VALUE(error_grad, p_i) += wGrad;
            DATA_VALUE(lamb, p_i) += dot(error_grad_j, error_grad_j);
        }

        DATA_VALUE(error_grad, p_i) /= CONST_VALUE(rest_density);

        DATA_VALUE(lamb, p_i) += dot(DATA_VALUE(error_grad, p_i), DATA_VALUE(error_grad, p_i));

        DATA_VALUE(lamb, p_i) = -DATA_VALUE(error, p_i) / (DATA_VALUE(lamb, p_i) + 1e-6f);
    }

    __global__ void
    compute_dx_cuda(PBFConstantParams *d_const,
                    PBFDynamicParams *d_data,
                    NeighborSearchUGConfig *d_nsConfig,
                    NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        auto pos_i = DATA_VALUE(pos, p_i);
        DATA_VALUE(dx, p_i) *= 0;

        FOR_EACH_NEIGHBOR_Pj() {
            auto pos_j = DATA_VALUE(pos, p_j);
            auto wGrad = CUBIC_KERNEL_GRAD();

            DATA_VALUE(dx, p_i) += (DATA_VALUE(lamb, p_i) + DATA_VALUE(lamb, p_j)) * wGrad;
        }

        DATA_VALUE(dx, p_i) /= CONST_VALUE(rest_density);
    }

    __global__ void
    apply_ext_force_cuda(PBFConstantParams *d_const,
                         PBFDynamicParams *d_data,
                         NeighborSearchUGConfig *d_nsConfig,
                         NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(vel, p_i) += CONST_VALUE(gravity) * CONST_VALUE(dt);
        DATA_VALUE(pos, p_i) += DATA_VALUE(vel, p_i) * CONST_VALUE(dt);
    }

    __global__ void
    apply_dx_cuda(PBFConstantParams *d_const,
                  PBFDynamicParams *d_data,
                  NeighborSearchUGConfig *d_nsConfig,
                  NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        DATA_VALUE(pos, p_i) += DATA_VALUE(dx, p_i);
        DATA_VALUE(vel, p_i) += DATA_VALUE(dx, p_i) / CONST_VALUE(dt);
    }

    __global__ void
    XSPH_cuda(PBFConstantParams *d_const,
              PBFDynamicParams *d_data,
              NeighborSearchUGConfig *d_nsConfig,
              NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD();

        if (DATA_VALUE(mat, p_i) != COMMON_NEWTON)
            return;

        auto pos_i = DATA_VALUE(pos, p_i);
        auto vel_i = DATA_VALUE(vel, p_i);
        DATA_VALUE(dx, p_i) *= 0;
        Vec3f dv;

        FOR_EACH_NEIGHBOR_Pj() {
            auto pos_j = DATA_VALUE(pos, p_j);
            auto vel_j = DATA_VALUE(vel, p_j);

            dv += CONST_VALUE(XSPH_k) * DATA_VALUE(volume, p_j) * (vel_j - vel_i) * CUBIC_KERNEL_VALUE();
        }

        DATA_VALUE(dx, p_i) += dv * CONST_VALUE(dt);
    }
}

namespace SoSim { // host invoke api
    __host__ void
    init_data(PBFConstantParams &h_const,
              PBFConstantParams *d_const,
              PBFDynamicParams *d_data,
              NeighborSearchUGConfig *d_nsConfig,
              NeighborSearchUGParams *d_nsParams) {
        // init_data_cuda
        init_cuda<<<h_const.block_num, h_const.thread_num>>>(d_const,
                                                             d_data);

        // compute_rigid_particle_volume_cuda
        compute_rigid_particle_volume_cuda<<<h_const.block_num, h_const.thread_num>>>(d_const,
                                                                                      d_data,
                                                                                      d_nsConfig,
                                                                                      d_nsParams);
    }

    __host__ void
    compute_sph_density_and_error(PBFConstantParams &h_const,
                                  PBFConstantParams *d_const,
                                  PBFDynamicParams *d_data,
                                  NeighborSearchUGConfig *d_nsConfig,
                                  NeighborSearchUGParams *d_nsParams) {
        // compute_sph_density_and_lamb_cuda
        compute_sph_density_and_error_cuda<<<h_const.block_num, h_const.thread_num>>>(d_const,
                                                                                      d_data,
                                                                                      d_nsConfig,
                                                                                      d_nsParams);
    }

    __host__ void
    update_lamb(PBFConstantParams &h_const,
                PBFConstantParams *d_const,
                PBFDynamicParams *d_data,
                NeighborSearchUGConfig *d_nsConfig,
                NeighborSearchUGParams *d_nsParams) {
        // update_lamb_cuda
        update_lamb_cuda<<<h_const.block_num, h_const.thread_num>>>(d_const,
                                                                    d_data,
                                                                    d_nsConfig,
                                                                    d_nsParams);
    }

    __host__ void
    compute_dx(PBFConstantParams &h_const,
               PBFConstantParams *d_const,
               PBFDynamicParams *d_data,
               NeighborSearchUGConfig *d_nsConfig,
               NeighborSearchUGParams *d_nsParams) {
        // compute_dx_cuda
        compute_dx_cuda<<<h_const.block_num, h_const.thread_num>>>(d_const,
                                                                   d_data,
                                                                   d_nsConfig,
                                                                   d_nsParams);
    }

    __host__ void
    apply_ext_force(PBFConstantParams &h_const,
                    PBFConstantParams *d_const,
                    PBFDynamicParams *d_data,
                    NeighborSearchUGConfig *d_nsConfig,
                    NeighborSearchUGParams *d_nsParams) {
        // apply_ext_force_cuda
        apply_ext_force_cuda<<<h_const.block_num, h_const.thread_num>>>(d_const,
                                                                        d_data,
                                                                        d_nsConfig,
                                                                        d_nsParams);
    }

    __host__ void
    apply_dx(PBFConstantParams &h_const,
             PBFConstantParams *d_const,
             PBFDynamicParams *d_data,
             NeighborSearchUGConfig *d_nsConfig,
             NeighborSearchUGParams *d_nsParams) {
        // apply_dx_cuda
        apply_dx_cuda<<<h_const.block_num, h_const.thread_num>>>(d_const,
                                                                 d_data,
                                                                 d_nsConfig,
                                                                 d_nsParams);
    }

    __host__ void
    post_correct(PBFConstantParams &h_const,
                 PBFConstantParams *d_const,
                 PBFDynamicParams *d_data,
                 NeighborSearchUGConfig *d_nsConfig,
                 NeighborSearchUGParams *d_nsParams) {
        // XSPH_cuda
        XSPH_cuda<<<h_const.block_num, h_const.thread_num>>>(d_const,
                                                             d_data,
                                                             d_nsConfig,
                                                             d_nsParams);

        // apply_dx_cuda
        apply_dx_cuda<<<h_const.block_num, h_const.thread_num>>>(d_const,
                                                                 d_data,
                                                                 d_nsConfig,
                                                                 d_nsParams);
    }
}