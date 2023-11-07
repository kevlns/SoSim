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
        float2 rest_density; // two phase rest density vec
        float rest_volume;  // TODO if is right?
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
        float3 *acc;
        float2 *volF_k; // two phase volume fraction vec
        float2 *d_volF;
        float2 *pressure_k; // two phase volume pressure
        float *pressure; // two phase volume pressure
        float *mass;
        float *density;
        float *density_ba;
        Material *mat;
        Phase *original_phase;


        void destroy() const{
            cudaFree(pos);
            cudaFree(predictPos);
            cudaFree(v_m);
            cudaFree(v_mk1);
            cudaFree(v_mk2);
            cudaFree(acc);
            cudaFree(volF_k);
            cudaFree(d_volF);
            cudaFree(pressure_k);
            cudaFree(pressure);
            cudaFree(mass);
            cudaFree(density);
            cudaFree(density_ba);
            cudaFree(mat);
            cudaFree(original_phase);

            cudaGetLastError_t("ERROR::MSMSPH::DynamicParams::destroy() failed.");
        }
    };

}

#endif //SOSIM_DATA_MSMSPH_PACK_HPP
