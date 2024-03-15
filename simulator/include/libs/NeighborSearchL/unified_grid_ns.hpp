//
// Created by ADMIN on 2024/3/5.
//

#ifndef SOSIM_UNIFIED_GRID_NS_HPP
#define SOSIM_UNIFIED_GRID_NS_HPP

#include <vector>

#include "core/math/matrix.hpp"

namespace SoSim {

    struct NeighborSearchUGConfig {
        Vec3f sceneLB;
        Vec3f sceneSize;
        Vec3ui gridSize;
        float cellLength;
        unsigned cellNum;
        unsigned particle_num;
        unsigned maxNeighborNum;

        unsigned kernel_blocks{0};
        unsigned kernel_threads{0};
    };

    struct NeighborSearchUGParams {
    public:
        Vec3i *cellOffsets_cuData{nullptr};
        unsigned *particleIndices_cuData{nullptr};
        unsigned *cellIndices_cuData{nullptr};
        unsigned *cellStart_cuData{nullptr};
        unsigned *cellEnd_cuData{nullptr};
        unsigned *neighborNum_cuData{nullptr};
        unsigned *neighbors_cuData{nullptr};

    private:
        bool isInit{false};

    public:
        void malloc(const NeighborSearchUGConfig &config);

        void freeMemory();
    };

    class NeighborSearchUG {
    public:
        void setConfig(NeighborSearchUGConfig config);

        void malloc();

        void update(Vec3f *pos_cuData);

        void dump() const;

        void freeMemory();

    public:
        NeighborSearchUGConfig h_config;
        NeighborSearchUGConfig *d_config{nullptr};
        NeighborSearchUGParams h_params;
        NeighborSearchUGParams *d_params{nullptr};

    private:
        double m_mem{0};

    };

}

#endif //SOSIM_UNIFIED_GRID_NS_HPP
