//@author        : Long Shen
//@date          : 2023/10/11
//@description   :
//@version       : 1.0

#ifndef SOSIM_NEIGHBOR_SEARCH_UGB_HPP
#define SOSIM_NEIGHBOR_SEARCH_UGB_HPP

#include "Public/Shared/NeighborSearchUGB/neighbor_search_config.hpp"
#include "Private/Shared/NeighborSearchUGB/data_pack.hpp"
#include "Private/BuildSystem/macro_helper.hpp"

namespace SoSim::NSUGB {

    class NeighborSearcher {

    public:
        SOSIM_API
        NeighborSearcher() = default;

        SOSIM_API
        ~NeighborSearcher() = default;

        SOSIM_API
        void initialize();

        SOSIM_API
        void update(float3 *device_pos);

        SOSIM_API
        void dumpInfo() const;

        SOSIM_API
        void setConfig(const NeighborSearchConfig *config);

        SOSIM_API
        uint32_t *getPartIndexDevicePtr() const;

        SOSIM_API
        uint32_t *getNeighborsDevicePtr() const;

        SOSIM_API
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
