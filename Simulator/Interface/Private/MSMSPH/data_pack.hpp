//@author        : Long Shen
//@date          : 2023/10/24
//@description   :
//@version       : 1.0

#ifndef SOSIM_DATA_MSMSPH_PACK_HPP
#define SOSIM_DATA_MSMSPH_PACK_HPP

#include <cuda_runtime.h>
#include <vector_types.h>
#include <iostream>

#include "Public/Shared/CudaUtils/cuda_tool.hpp"
#include "Public/Framework/mat.hpp"

namespace SoSim::MSMSPH {

    struct ConstParams {
        float3 gravity;
        float2 rest_density; // two phase rest density_m vec
        float rest_volume;
        uint32_t total_particle_num;
        uint32_t ns_maxNeighborNum;
        float sph_h;

        float dt;
    };

    struct DynamicParams {
        float3 *pos;
        float3 *predictPos;
        float3 *v_m;
        float3 *v_mk1; // drift vel 1
        float3 *v_mk2; // drift vel 2
        float3 *acc_m;
        float3 *M_m;
        float2 *alpha_k; // two phase volume fraction vec
        float2 *delta_alpha;
        float *pressure_m; // two phase volume pressure_m
        float2 *pressure_k; // two phase volume pressure_m
        float *density_m;
        float *density_sph;
        Material *mat;
        Phase *original_phase;


        void destroy() const {
            cudaFree(pos);
            cudaFree(predictPos);
            cudaFree(v_m);
            cudaFree(v_mk1);
            cudaFree(v_mk2);
            cudaFree(acc_m);
            cudaFree(M_m);
            cudaFree(alpha_k);
            cudaFree(delta_alpha);
            cudaFree(pressure_k);
            cudaFree(pressure_m);
            cudaFree(density_m);
            cudaFree(density_sph);
            cudaFree(mat);
            cudaFree(original_phase);

            cudaGetLastError_t("ERROR::MSMSPH::DynamicParams::destroy() failed.");
        }
    };

}

#endif //SOSIM_DATA_MSMSPH_PACK_HPP