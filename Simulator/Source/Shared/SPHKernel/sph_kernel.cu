//@author        : Long Shen
//@date          : 2023/10/26
//@description   :
//@version       : 1.0

#include "Public/Shared/SPHKernel/sph_kernel.cuh"
#include "Public/Shared/Math/helper_math.hpp"

namespace SoSim {

    inline __device__ float cubic_value(const float3 &x_ij, float sph_h) {
        const float PI = 3.14159265;
        const float cubicSigma = 8.f / PI / static_cast<float>(std::pow(sph_h, 3));
        const float r_norm = length(x_ij);

        float res = 0.0;
        float invH = 1 / sph_h;
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

    inline __device__ float3 cubic_grad(const float3 &x_ij, float sph_h) {
        const float PI = 3.14159265;
        const float cubicSigma = 8.f / PI / static_cast<float>(std::pow(sph_h, 3));

        auto res = make_float3(0, 0, 0);
        float invH = 1 / sph_h;
        float q = length(x_ij) * invH;

        if (q < 1e-6 || q > 1)
            return res;

        float3 grad_q = x_ij / (length(x_ij) * sph_h);
        if (q <= 0.5)
            res = (6 * (3 * q * q - 2 * q)) * grad_q * cubicSigma;
        else {
            auto factor = 1 - q;
            res = -6 * factor * factor * grad_q * cubicSigma;
        }

        return res;
    }

    inline __device__ float poly6_value(const float3 &x_ij, float sph_h) {
        const float PI = 3.14159265;
        const float sigma = 315 / (64 * PI * pow(sph_h, 9));
        const float r = length(x_ij);

        float res = 0;

        if (r < 1e-6 || r > sph_h)
            return res;

        res = pow((sph_h * sph_h - r * r), 3) * sigma;
        return res;
    }

    inline __device__ float3 spiky_grad(const float3 &x_ij, float sph_h) {
        const float PI = 3.14159265;
        const float sigma = 315 / (64 * PI * pow(sph_h, 9));
        const float r = length(x_ij);

        float3 res = {0, 0, 0};

        if (r < 1e-6 || r > sph_h)
            return res;

        res = -3 * pow(sph_h - r, 2) * sigma * normalize(x_ij);
        return res;
    }

}