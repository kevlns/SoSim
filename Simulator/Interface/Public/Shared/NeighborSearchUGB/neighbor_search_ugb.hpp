//@author        : Long Shen
//@date          : 2023/10/11
//@description   :
//@version       : 1.0

#ifndef SOSIM_NEIGHBOR_SEARCH_UGB_HPP
#define SOSIM_NEIGHBOR_SEARCH_UGB_HPP

#include "Public/Shared/NeighborSearchUGB/neighbor_search_config.hpp"
#include "Private/Shared/NeighborSearchUGB/data_pack.hpp"

namespace SoSim::NSUGB {

    class NeighborSearcher {

    public:
        NeighborSearcher() = default;

        ~NeighborSearcher() = default;

        void initialize(float3 scene_lb, float3 scene_size, unsigned total_particle_num, float sph_support_radius);

        void update(float3 *device_pos);

        void dumpInfo() const;

        void setConfig(const NeighborSearchConfig *config);

        uint32_t *getPartIndexDevicePtr() const;

        uint32_t *getNeighborsDevicePtr() const;

        void destroy();

    private:
        double m_mem{0};
        bool m_isInit{false};
        NeighborSearchConfig *m_config{nullptr};
        ConstParams m_host_cp{};
        DynamicParams m_host_dp{};

    private:
        ConstParams *m_device_cp{nullptr};
        DynamicParams *m_device_dp{nullptr};
    };

} // namespace SoSim::NSUGB

#endif // SOSIM_NEIGHBOR_SEARCH_UGB_HPP
