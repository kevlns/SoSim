//
// Created by ADMIN on 2024/3/13.
//

#include <vector>
#include <cuda_runtime.h>

#include "libs/AnalysisL/statistic_util.hpp"

namespace SoSim {

    template<typename T>
    T dump_mean(T *d_ptr, unsigned raw_size, unsigned target_size, unsigned target_start) {
        std::vector<T> data(raw_size);
        cudaMemcpy(data.data(), d_ptr, raw_size * sizeof(T), cudaMemcpyDeviceToHost);

        T avg;
        if (std::is_same<float, T>::value || std::is_same<double, T>::value ||
            std::is_same<int, T>::value || std::is_same<unsigned, T>::value)
            avg = 0;
        for (int i = 0; i < raw_size; ++i) {
            if (i >= target_start && i < target_start + target_size)
                avg += data[i];
        }

        // TODO log system
        std::cout << "Target Mean: " << avg / target_size << "\n";

        return avg;
    }

    template<typename T>
    T dump_max(T *d_ptr, unsigned raw_size, unsigned target_size, unsigned target_start) {
        std::vector<T> data(raw_size);
        cudaMemcpy(data.data(), d_ptr, raw_size * sizeof(T), cudaMemcpyDeviceToHost);

        T max_v;
        if (std::is_same<float, T>::value || std::is_same<double, T>::value ||
            std::is_same<int, T>::value || std::is_same<unsigned, T>::value)
            max_v = 0;
        for (int i = 0; i < raw_size; ++i) {
            if (i >= target_start && i < target_start + target_size)
                if (max_v < data[i])
                    max_v = data[i];
        }

        // TODO log system
        std::cout << "Target Max: " << max_v << "\n";

        return max_v;
    }

    template<typename T>
    T cal_mean(T *d_ptr, unsigned raw_size, unsigned target_size, unsigned target_start) {
        std::vector<T> data(raw_size);
        cudaMemcpy(data.data(), d_ptr, raw_size * sizeof(T), cudaMemcpyDeviceToHost);

        T mean;
        if (std::is_same<float, T>::value || std::is_same<double, T>::value ||
            std::is_same<int, T>::value || std::is_same<unsigned, T>::value)
            mean = 0;
        for (int i = 0; i < raw_size; ++i) {
            if (i >= target_start && i < target_start + target_size)
                mean += data[i];
        }

        return mean / data.size();
    }

    // explicit template instance
    template Vec3f dump_mean<Vec3f>(Vec3f *d_ptr, unsigned raw_size, unsigned target_size, unsigned target_start);

    template float dump_mean<float>(float *d_ptr, unsigned raw_size, unsigned target_size, unsigned target_start);

    template unsigned
    dump_mean<unsigned>(unsigned *d_ptr, unsigned raw_size, unsigned target_size, unsigned target_start);

    template Vec3f dump_max<Vec3f>(Vec3f *d_ptr, unsigned raw_size, unsigned target_size, unsigned target_start);

    template float dump_max<float>(float *d_ptr, unsigned raw_size, unsigned target_size, unsigned target_start);

    template unsigned
    dump_max<unsigned>(unsigned *d_ptr, unsigned raw_size, unsigned target_size, unsigned target_start);

    template float cal_mean<float>(float *d_ptr, unsigned raw_size, unsigned target_size, unsigned target_start);
}