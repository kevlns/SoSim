//
// Created by ADMIN on 2024/5/31.
//

#ifndef SOSIM_PARTICLE_EMITTER_HPP
#define SOSIM_PARTICLE_EMITTER_HPP

#include <optional>
#include <string>

#include "core/math/matrix.hpp"
#include "framework/MetaFramework/object.hpp"

namespace SoSim {

    struct ParticleEmitterConfig {
        // required settings
        std::optional<float> particle_radius;
        std::optional<bool> use_unified_buffer;
        std::optional<bool> use_emit_agent;
        std::optional<Material> emit_mat;
        unsigned max_particle_num{10000};
        float emit_vel{1};

        // if multiphase
        std::vector<float> phases;

        // if use_unified_buffer = true, then set these three terms
        Vec3f *unified_pos_buffer{nullptr};
        Vec3f *unified_vel_buffer{nullptr};
        Material *unified_mat_buffer{nullptr};
        unsigned insert_index{0};

        // if use_agent_file = true, then set this term
        std::string agent_file;
        std::string agent_normal_file;
        // else
        ParticleObjectConfig agent_obj_config;
    };

    class ParticleEmitter {
    public:
        ParticleEmitter();

        ~ParticleEmitter();

        std::shared_ptr<ParticleEmitterConfig> getConfig();

        Vec3f *getCudaTemplatePosBuffer();

        Vec3f *getCudaTemplateVelBuffer();

        unsigned getTemplatePartNum() const;

        void update();

        bool emit(float t);

    private:
        std::shared_ptr<ParticleEmitterConfig> m_config;
        unsigned m_emit_times{0};
        float m_emit_gap{1e-6};
        unsigned m_template_parts_num;

        Vec3f *m_device_cuda_pos_buffer{nullptr};
        Vec3f *m_device_cuda_vel_buffer{nullptr};
        Material *m_device_cuda_mat_buffer{nullptr};
        Vec3f *m_device_cuda_T_pos{nullptr};
        Vec3f *m_device_cuda_T_vel{nullptr};
        Material *m_device_cuda_T_mat{nullptr};
        unsigned m_insert_index{0};
    };

}

#endif //SOSIM_PARTICLE_EMITTER_HPP
