//@author        : Long Shen
//@date          : 2023/10/27
//@description   :
//@version       : 1.0

#include <iostream>
#include <cuda_runtime.h>

#include "Public/Shared/Analyzer/params_analysis_helper.hpp"
#include "Public/Shared/Math/helper_math.hpp"

namespace SoSim {


    template<typename ParamPType>
    extern void dump_avg(unsigned int num, ParamPType d_params) {
        if (std::is_same<ParamPType, float *>::value) {
            auto *c_params = new float[num];
            cudaMemcpy(c_params, d_params, num * sizeof(float), cudaMemcpyDeviceToHost);

            float avg = 0.f;
            for (int i = 0; i < num; ++i)
                avg += c_params[i];

            if (avg < 1e-6)
                avg = 0;

            std::cout << "dump_avg():: The specified param's AVG_VALUE: " << avg / static_cast<float >(num) << "\n";

            delete[] c_params;
        } else if (std::is_same<ParamPType, unsigned int *>::value) {
            auto *c_params = new unsigned int[num];
            cudaMemcpy(c_params, d_params, num * sizeof(unsigned int), cudaMemcpyDeviceToHost);

            float avg = 0.f;
            for (int i = 0; i < num; ++i)
                avg += c_params[i];

            std::cout << "dump_avg():: The specified param's AVG_VALUE: " << avg / static_cast<float >(num) << "\n";

            delete[] c_params;

        } else if (std::is_same<ParamPType, float3 *>::value) {
            auto *c_params = new float3[num];
            cudaMemcpy(c_params, d_params, num * sizeof(float3), cudaMemcpyDeviceToHost);

            float3 avg = {0, 0, 0};
            for (int i = 0; i < num; ++i)
                avg += c_params[i] / num;

            auto res = avg;
            if (res.x < 1e-6)
                res.x = 0;
            if (res.y < 1e-6)
                res.y = 0;
            if (res.z < 1e-6)
                res.z = 0;

            std::cout << "dump_avg():: The specified param's AVG_VALUE: " << "[" << res.x << "," << res.y << ","
                      << res.z << "]" << "\n";

            delete[] c_params;

        } else {
            throw std::runtime_error("dump_avg():: Param type not supported!");
        }
    }

    template void dump_avg<float *>(unsigned int num, float *d_params);

    template void dump_avg<float3 *>(unsigned int num, float3 *d_params);

    template void dump_avg<float2 *>(unsigned int num, float2 *d_params);

    template void dump_avg<uint32_t *>(unsigned int num, uint32_t *d_params);

    template<typename ParamPType>
    void dump_max(unsigned int num, ParamPType d_params) {
        if (std::is_same<ParamPType, float *>::value) {
            auto *c_params = new float[num];
            cudaMemcpy(c_params, d_params, num * sizeof(float), cudaMemcpyDeviceToHost);

            float max_value = 0.f;
            for (int i = 0; i < num; ++i)
                max_value = std::max(max_value, c_params[i]);

            std::cout << "dump_max():: The specified param's MAX_VALUE: " << max_value << "\n";

            delete[] c_params;
        } else if (std::is_same<ParamPType, unsigned int *>::value) {
            auto *c_params = new unsigned int[num];
            cudaMemcpy(c_params, d_params, num * sizeof(unsigned int), cudaMemcpyDeviceToHost);

            unsigned max_value = 0.f;
            for (int i = 0; i < num; ++i)
                max_value = std::max(max_value, c_params[i]);

            std::cout << "dump_max():: The specified param's MAX_VALUE: " << max_value << "\n";

            delete[] c_params;

        } else if (std::is_same<ParamPType, float3 *>::value) {
            auto *c_params = new float3[num];
            cudaMemcpy(c_params, d_params, num * sizeof(float3), cudaMemcpyDeviceToHost);

            float3 max_value = {0, 0, 0};
            for (int i = 0; i < num; ++i) {
                if (length(c_params[i]) > length(max_value))
                    max_value = c_params[i];
            }

            auto res = max_value;
            std::cout << "dump_max():: The specified param's MAX_VALUE: " << "[" << res.x << "," << res.y << ","
                      << res.z << "]" << "\n";

            delete[] c_params;

        } else {
            throw std::runtime_error("dump_max():: Param type not supported!");
        }
    }

    template void dump_max<float *>(unsigned int num, float *d_params);

    template void dump_max<float3 *>(unsigned int num, float3 *d_params);

    template void dump_max<float2 *>(unsigned int num, float2 *d_params);

    template void dump_max<uint32_t *>(unsigned int num, uint32_t *d_params);

    template<typename ParamPType>
    void dump_min(unsigned int num, ParamPType d_params) {
        if (std::is_same<ParamPType, float *>::value) {
            auto *c_params = new float[num];
            cudaMemcpy(c_params, d_params, num * sizeof(float), cudaMemcpyDeviceToHost);

            float min_value = 0xffffff;
            for (int i = 0; i < num; ++i)
                min_value = std::min(min_value, c_params[i]);

            std::cout << "dump_min():: The specified param's MIN_VALUE: " << min_value << "\n";

            delete[] c_params;
        } else if (std::is_same<ParamPType, unsigned int *>::value) {
            auto *c_params = new unsigned int[num];
            cudaMemcpy(c_params, d_params, num * sizeof(unsigned int), cudaMemcpyDeviceToHost);

            unsigned min_value = 0xffffff;
            for (int i = 0; i < num; ++i)
                min_value = std::min(min_value, c_params[i]);

            std::cout << "dump_min():: The specified param's MIN_VALUE: " << min_value << "\n";

            delete[] c_params;

        } else {
            throw std::runtime_error("dump_min():: Param type not supported!");
        }
    }

}
