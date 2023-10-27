//@author        : Long Shen
//@date          : 2023/10/7
//@description   :
//@version       : 1.0

#ifndef PHYSIKA_CUDA_TOOL_CUH
#define PHYSIKA_CUDA_TOOL_CUH

#include <iostream>
#include <cuda_runtime.h>

/**
 *
 * @param msg user-defined error message
 * @return if cuda function failed, such as cudaMalloc() and cudaFree(), it will output the error msg and return false; otherwise, return true.
 */
static bool cudaGetLastError_t(const std::string& msg)
{
    auto c_error = cudaGetLastError();
    if (cudaSuccess != c_error)
    {
        std::cerr << "error:: cudaGetLastError:: "
                  << msg << " : "
                  << "(" << static_cast<int>(c_error) << ") "
                  << cudaGetErrorString(c_error) << ".\n";
        return false;
    }
    return true;
}

/**
 *
 * @param pPtr device pointer to free
 * @param size size of space that pointer occupies
 * @param mem a counter to restore allocated mem
 */
static void cudaMalloc_t(void** pPtr, size_t size, double& mem)
{
    cudaMalloc(pPtr, size);
    mem += static_cast<double>(size) / 1000000;
}

#endif  // PHYSIKA_CUDA_TOOL_CUH
