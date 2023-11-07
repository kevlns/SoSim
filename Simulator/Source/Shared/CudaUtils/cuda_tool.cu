//@author        : Long Shen
//@date          : 2023/10/27
//@description   :
//@version       : 1.0

#include "Public/Shared/CudaUtils/cuda_tool.hpp"

namespace SoSim {

    extern bool cudaGetLastError_t(const std::string &msg) {
        auto c_error = cudaGetLastError();
        if (cudaSuccess != c_error) {
            std::cerr << "error:: cudaGetLastError:: "
                      << msg << " : "
                      << "(" << static_cast<int>(c_error) << ") "
                      << cudaGetErrorString(c_error) << ".\n";
            return false;
        }
        return true;
    }

    extern void cudaMalloc_t(void **pPtr, size_t size, double &mem) {
        cudaMalloc(pPtr, size);
        mem += static_cast<double>(size) / 1000000;
    }

}




