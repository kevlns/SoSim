//@author        : Long Shen
//@date          : 2023/11/8
//@description   :
//@version       : 1.0

#ifndef SOSIM_MMSPH_DATA_PACK_HPP
#define SOSIM_MMSPH_DATA_PACK_HPP

#include "Public/Shared/CudaUtils/cuda_tool.hpp"
#include "Public/Framework/mat.hpp"

namespace SoSim::MMSPH {

    struct ConstParams {
        float dt;
        float rest_volume;
        float sph_h;
        uint32_t total_particle_num;
        uint32_t max_neighbor_num;
        float3 gravity;
        float2 rest_density; // two-phase particle
    };

    struct DynamicParams {
        float3 *pos;
        float3 *predictPos;
        float3 *acc_mix;
        float3 *vel_mix;
        float3 *drift_vel_k1;
        float3 *drift_vel_k2;
        float3 *M_m;
        float2 *alpha;
        float *density_mix;
        float *density_sph;
        Material *mat;

        void destroy() const {
            cudaFree(pos);
            cudaFree(predictPos);
            cudaFree(acc_mix);
            cudaFree(vel_mix);
            cudaFree(drift_vel_k1);
            cudaFree(drift_vel_k2);
            cudaFree(M_m);
            cudaFree(density_mix);
            cudaFree(density_sph);
            cudaFree(mat);

            cudaGetLastError_t("ERROR::MMSPH::DynamicParams::destroy() failed.");
        }
    };

}

#endif //SOSIM_MMSPH_DATA_PACK_HPP
