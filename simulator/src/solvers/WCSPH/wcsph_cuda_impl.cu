//
// Created by ADMIN on 2024/3/14.
//

#include "macro.hpp"
#include "wcsph_cuda_api.cuh"

namespace SoSim {
    /*
     *  Cuda API
     */
    __device__ inline float
    cubic_value(const Vec3f &r, float h) {
        float r_norm = r.length();
        const float PI = 3.14159265;
        const float cubicSigma = 8.f / PI / static_cast<float>(std::pow(h, 3));

        float res = 0.0;
        float invH = 1 / h;
        float q = r_norm * invH;

        if (q <= 1) {
            if (q <= 0.5) {
                auto q2 = q * q;
                auto q3 = q2 * q;
                res = static_cast<float>(cubicSigma * (6.0 * q3 - 6.0 * q2 + 1));
            } else {
                res = static_cast<float>(cubicSigma * 2 * std::pow(1 - q, 3));
            }
        }

        return res;
    }

    __device__ inline Vec3f
    cubic_gradient(const Vec3f &r, float h) {
        const float PI = 3.14159265;
        const float cubicSigma = 8.f / PI / static_cast<float>(std::pow(h, 3));

        auto res = Vec3f();
        float invH = 1 / h;
        float q = r.length() * invH;

        if (q < 1e-6 || q > 1)
            return res;

        Vec3f grad_q = r / (r.length() * h);
        if (q <= 0.5)
            res = (6 * (3 * q * q - 2 * q)) * grad_q * cubicSigma;
        else {
            auto factor = 1 - q;
            res = -6 * factor * factor * grad_q * cubicSigma;
        }

        return res;
    }

    __global__ void
    init_cuda(WCSPHConstantParams *d_const,
              WCSPHDynamicParams *d_data,
              NeighborSearchUGConfig *d_nsConfig,
              NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD()

        DATA_VALUE(density, p_i) = CONST_VALUE(rest_density);
        DATA_VALUE(volume, p_i) = CONST_VALUE(rest_volume);
        DATA_VALUE(pos_adv, p_i) = DATA_VALUE(pos, p_i);

        // compute rigid particle volume
        if (DATA_VALUE(mat, p_i) == FIXED_BOUND || DATA_VALUE(mat, p_i) == DYNAMIC_RIGID) {
            float density = 0;
            auto pos_i = DATA_VALUE(pos, p_i);
            FOR_EACH_NEIGHBOR_Pj() {
                auto pos_j = DATA_VALUE(pos, p_j);

                if (DATA_VALUE(mat, p_j) == DATA_VALUE(mat, p_i))
                    density += cubic_value(pos_i - pos_j, CONST_VALUE(h));
            }
            DATA_VALUE(volume, p_i) = max(1 / density, CONST_VALUE(rest_volume));
            DATA_VALUE(density, p_i) = CONST_VALUE(rest_rigid_density);
        }
    }

    __global__ void
    computeDensityAndPressure_cuda(WCSPHConstantParams *d_const,
                                   WCSPHDynamicParams *d_data,
                                   NeighborSearchUGConfig *d_nsConfig,
                                   NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD()

        if (DATA_VALUE(mat, p_i) != FLUID)
            return;

        // define your local var
        auto pos_i = DATA_VALUE(pos_adv, p_i);
        float density_sph = 0;

        FOR_EACH_NEIGHBOR_Pj() {
            auto pos_j = DATA_VALUE(pos_adv, p_j);

            float dens = DATA_VALUE(density, p_j);
            if (DATA_VALUE(mat, p_j) == FIXED_BOUND || DATA_VALUE(mat, p_j) == DYNAMIC_RIGID)
                dens = DATA_VALUE(density, p_i);

            density_sph += dens * cubic_value(pos_i - pos_j, CONST_VALUE(h));
        }

        DATA_VALUE(density, p_i) = max(density_sph, CONST_VALUE(rest_density));
        DATA_VALUE(pressure, p_i) =
                CONST_VALUE(stiff) * (pow(DATA_VALUE(density, p_i) / CONST_VALUE(rest_density), 7) - 1);
    }

    __global__ void
    computeGravityForce_cuda(WCSPHConstantParams *d_const,
                             WCSPHDynamicParams *d_data,
                             NeighborSearchUGConfig *d_nsConfig,
                             NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD()

        if (DATA_VALUE(mat, p_i) != FLUID)
            return;

        DATA_VALUE(vel_adv, p_i) += CONST_VALUE(dt) * CONST_VALUE(gravity);
    }

    __global__ void
    computePressureForce_cuda(WCSPHConstantParams *d_const,
                              WCSPHDynamicParams *d_data,
                              NeighborSearchUGConfig *d_nsConfig,
                              NeighborSearchUGParams *d_nsParams) {

        CHECK_THREAD()

        if (DATA_VALUE(mat, p_i) != FLUID)
            return;

        auto pos_i = DATA_VALUE(pos_adv, p_i);
        auto pressure_i = DATA_VALUE(pressure, p_i);
        Vec3f acc;
        FOR_EACH_NEIGHBOR_Pj() {
            auto pos_j = DATA_VALUE(pos_adv, p_j);

            auto pressure_j = DATA_VALUE(pressure, p_j);
            if (DATA_VALUE(mat, p_j) == FIXED_BOUND || DATA_VALUE(mat, p_j) == DYNAMIC_RIGID)
                pressure_j = pressure_i;

            acc += -DATA_VALUE(density, p_j) * DATA_VALUE(volume, p_j) *
                   (pressure_i / pow(DATA_VALUE(density, p_i), 2) + pressure_j / pow(DATA_VALUE(density, p_j), 2)) *
                   cubic_gradient(pos_i - pos_j, CONST_VALUE(h));
        }
        DATA_VALUE(vel_adv, p_i) += CONST_VALUE(dt) * acc;
    }

    __global__ void
    computeViscousForce_cuda(WCSPHConstantParams *d_const,
                             WCSPHDynamicParams *d_data,
                             NeighborSearchUGConfig *d_nsConfig,
                             NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD()

        if (DATA_VALUE(mat, p_i) != FLUID)
            return;

        auto pos_i = DATA_VALUE(pos_adv, p_i);
        auto vel_i = DATA_VALUE(vel, p_i);
        Vec3f acc;
        FOR_EACH_NEIGHBOR_Pj() {
            auto pos_j = DATA_VALUE(pos_adv, p_j);
            auto vel_j = DATA_VALUE(vel, p_j);

            acc += 10 * CONST_VALUE(rest_vis) * DATA_VALUE(volume, p_j) *
                   dot(vel_i - vel_j, pos_i - pos_j) * cubic_gradient(pos_i - pos_j, CONST_VALUE(h)) /
                   (pow((pos_i - pos_j).length(), 2) + 1e-6f);
        }
        DATA_VALUE(vel_adv, p_i) += CONST_VALUE(dt) * acc;
    }

    __global__ void
    advect_cuda(WCSPHConstantParams *d_const,
                WCSPHDynamicParams *d_data,
                NeighborSearchUGConfig *d_nsConfig,
                NeighborSearchUGParams *d_nsParams) {
        CHECK_THREAD()

        if (DATA_VALUE(mat, p_i) != FLUID)
            return;

        DATA_VALUE(pos_adv, p_i) += CONST_VALUE(dt) * DATA_VALUE(vel_adv, p_i);
        DATA_VALUE(pos, p_i) = DATA_VALUE(pos_adv, p_i);
        DATA_VALUE(vel, p_i) = DATA_VALUE(vel_adv, p_i);
    }

}