//
// Created by ADMIN on 2024/3/22.
//

#ifndef SOSIM_KERNELS_CUH
#define SOSIM_KERNELS_CUH

#include "core/math/matrix.hpp"

namespace SoSim {
    __device__ inline float
    cubic_value(const Vec3f &r, float h) {
        const float r_norm = r.length();
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

    __host__ __device__ inline float
    surface_tension_C(const float r_norm, const float h) {
        const float PI = 3.14159265;
        const float cSigma = 32.f / PI / static_cast<float>(std::pow(h, 9));

        if (r_norm * 2 > h && r_norm <= h)
            return cSigma * std::pow(h - r_norm, 3) * std::pow(r_norm, 3);
        else if (r_norm > 0 && r_norm * 2 <= h)
            return cSigma * (2 * std::pow(h - r_norm, 3) * std::pow(r_norm, 3) - std::pow(h, 6) / 64);

        return 0;
    }

    __host__ __device__ inline float
    df_viscosity_kernel_laplacian(const Vec3f &r, const float h) {
        const float PI = 3.14159265;
        const float r_norm = r.length();
        return (r_norm <= h) ? (45.0f * (h - r_norm) / (PI * powf(h, 6))) : 0.0f;
    }
}

#endif //SOSIM_KERNELS_CUH
