//@author        : Long Shen
//@date          : 2023/10/26
//@description   :
//@version       : 1.0

#ifndef SOSIM_SPH_KERNEL_CUH
#define SOSIM_SPH_KERNEL_CUH

namespace SoSim {

    extern inline __device__ float cubic_value(const float3 &x_ij, float sph_h);

    extern inline __device__ float3 cubic_grad(const float3 &x_ij, float sph_h);

    extern inline __device__ float poly6_value(const float3 &x_ij, float sph_h);

    extern inline __device__ float3 spiky_grad(const float3 &x_ij, float sph_h);

}

#endif //SOSIM_SPH_KERNEL_CUH
