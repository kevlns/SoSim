//
// Created by ADMIN on 2024/3/5.
//
#include "libs/SPHL/kernels.cuh"

namespace SoSim {

    __host__ __device__ inline float
    cubic_value(const float r_norm, const float h) {
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

    __host__ __device__ inline Vec3f
    cubic_gradient(const Vec3f &r, const float h) {
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

}