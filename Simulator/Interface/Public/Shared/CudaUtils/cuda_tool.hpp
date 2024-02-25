//@author        : Long Shen
//@date          : 2023/10/7
//@description   :
//@version       : 1.0

#ifndef SOSIM_CUDA_TOOL_CUH
#define SOSIM_CUDA_TOOL_CUH

#include <iostream>
#include <cuda_runtime.h>

#include "Private/BuildSystem/macro_helper.hpp"

namespace SoSim {
/**
 *
 * @param msg user-defined error message
 * @return if cuda function failed, such as cudaMalloc() and cudaFree(), it will output the error msg and return false; otherwise, return true.
 */
    SOSIM_API
    extern bool cudaGetLastError_t(const std::string &msg);

/**
 *
 * @param pPtr device pointer to free
 * @param size size of space that pointer occupies
 * @param mem a counter to restore allocated mem
 */
    SOSIM_API
    extern void cudaMalloc_t(void **pPtr, size_t size, double &mem);

}
#endif  // SOSIM_CUDA_TOOL_CUH
