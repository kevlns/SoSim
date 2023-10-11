//@author        : Long Shen
//@date          : 2023/10/11
//@description   : 
//@version       : 1.0

#ifndef SOSIM_NEIGHBOR_SEARCH_UGB_HPP
#define SOSIM_NEIGHBOR_SEARCH_UGB_HPP

namespace SoSim::NSUGB {

    class NeighborSearchUGB {

    public:
        NeighborSearchUGB() = default;

        ~NeighborSearchUGB() = default;

        void initialize();

//        void update(float3 *device_pos);

        void destroy();

    };

}

#endif //SOSIM_NEIGHBOR_SEARCH_UGB_HPP
