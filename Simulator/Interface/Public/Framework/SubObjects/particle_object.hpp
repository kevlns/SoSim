//@author        : Long Shen
//@date          : 2023/11/22
//@description   :
//@version       : 1.0

#ifndef SOSIM_PARTICLE_OBJECT_HPP
#define SOSIM_PARTICLE_OBJECT_HPP

#include <vector>

#include "Public/Framework/object.hpp"
#include "Public/Framework/framework_config.hpp"

namespace SoSim {

    class ParticleObject : Object {
    public:

        explicit ParticleObject(ParticleObjectConfig &config);

        ~ParticleObject() override = default;

        void destroy() override;

        void refresh() override;

    private:

        void gen();

    public:
        ParticleObjectConfig m_config{};
        float3 m_vel_start{0, 0, 0};

    private:
        std::vector<float3> m_pos;
    };

}

#endif //SOSIM_PARTICLE_OBJECT_HPP
