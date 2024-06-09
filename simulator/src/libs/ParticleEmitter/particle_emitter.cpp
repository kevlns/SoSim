//
// Created by ADMIN on 2024/5/31.
//

#include "libs/ParticleEmitter/particle_emitter.hpp"

#include <cuda_runtime.h>

#include "libs/ModelL/model_helper.hpp"

namespace SoSim {

    ParticleEmitter::ParticleEmitter() {
        m_config = std::make_shared<ParticleEmitterConfig>();
    }

    ParticleEmitter::~ParticleEmitter() {
        if (!m_config->use_unified_buffer) {
            if (m_device_cuda_pos_buffer)
                cudaFree(m_device_cuda_pos_buffer);

            if (m_device_cuda_vel_buffer)
                cudaFree(m_device_cuda_vel_buffer);

            if (m_device_cuda_mat_buffer)
                cudaFree(m_device_cuda_mat_buffer);

            if (m_device_cuda_T_pos)
                cudaFree(m_device_cuda_T_pos);

            if (m_device_cuda_T_vel)
                cudaFree(m_device_cuda_T_vel);

            if (m_device_cuda_T_mat)
                cudaFree(m_device_cuda_T_mat);
        }
    }

    std::shared_ptr<ParticleEmitterConfig> ParticleEmitter::getConfig() {
        return m_config;
    }

    void ParticleEmitter::update() {
        if (m_config->use_unified_buffer) {
            m_device_cuda_pos_buffer = m_config->unified_pos_buffer;
            m_device_cuda_vel_buffer = m_config->unified_vel_buffer;
            m_device_cuda_mat_buffer = m_config->unified_mat_buffer;
        } else {
            if (m_device_cuda_pos_buffer)
                cudaFree(m_device_cuda_pos_buffer);

            if (m_device_cuda_vel_buffer)
                cudaFree(m_device_cuda_vel_buffer);

            size_t s = m_config->max_particle_num * sizeof(Vec3f);
            cudaMalloc((void **) &m_device_cuda_pos_buffer, s);
            cudaMalloc((void **) &m_device_cuda_vel_buffer, s);
        }

        if (m_device_cuda_T_pos)
            cudaFree(m_device_cuda_T_pos);

        if (m_device_cuda_T_vel)
            cudaFree(m_device_cuda_T_vel);

        if (m_device_cuda_T_mat)
            cudaFree(m_device_cuda_T_mat);

        m_insert_index = m_config->insert_index;
        m_emit_gap = (m_config->particle_radius.value() * 2.5f) / m_config->emit_vel;

        // TODO load template particles
        Vec3f dir;
        m_config->agent_obj_config.particle_radius = m_config->particle_radius;
        if (m_config->use_emit_agent) {
            m_config->agent_obj_config.model_file = m_config->agent_file;
            dir = ModelHelper::loadEmitterAgentNormal(m_config->agent_normal_file);
        }

        auto pos = ModelHelper::create3DParticleModel(m_config->agent_obj_config);
        m_template_parts_num = pos.size();
        Vec3f v = m_config->emit_vel * dir;
        std::vector<Vec3f> vel(pos.size(), v);

        if (!m_config->emit_mat.has_value())
            throw std::runtime_error("Error:: ParticleEmitter need to set emit material.");
        Material m = m_config->emit_mat.value();
        std::vector<Material> mat(pos.size(), m);

        cudaMalloc((void **) &m_device_cuda_T_pos, pos.size() * sizeof(Vec3f));
        cudaMalloc((void **) &m_device_cuda_T_vel, pos.size() * sizeof(Vec3f));
        cudaMalloc((void **) &m_device_cuda_T_mat, pos.size() * sizeof(Material));

        cudaMemcpy(m_device_cuda_T_pos, pos.data(), pos.size() * sizeof(Vec3f), cudaMemcpyHostToDevice);
        cudaMemcpy(m_device_cuda_T_vel, vel.data(), pos.size() * sizeof(Vec3f), cudaMemcpyHostToDevice);
        cudaMemcpy(m_device_cuda_T_mat, mat.data(), pos.size() * sizeof(Material), cudaMemcpyHostToDevice);

        cudaGetLastError();
    }

    bool ParticleEmitter::emit(float t) {
        static size_t size = m_template_parts_num * sizeof(Vec3f);
        static size_t size_mat = m_template_parts_num * sizeof(Material);
        if (t >= m_emit_gap * m_emit_times && m_emit_times * m_template_parts_num < m_config->max_particle_num) {
            cudaMemcpy(m_device_cuda_pos_buffer + m_insert_index, m_device_cuda_T_pos, size, cudaMemcpyDeviceToDevice);
            cudaMemcpy(m_device_cuda_vel_buffer + m_insert_index, m_device_cuda_T_vel, size, cudaMemcpyDeviceToDevice);
            cudaMemcpy(m_device_cuda_mat_buffer + m_insert_index, m_device_cuda_T_mat, size_mat,
                       cudaMemcpyDeviceToDevice);

            cudaGetLastError();

            m_insert_index += m_template_parts_num;
            m_emit_times += 1;

            return true;
        }

        return false;
    }

    Vec3f *ParticleEmitter::getCudaTemplatePosBuffer() {
        return m_device_cuda_T_pos;
    }

    Vec3f *ParticleEmitter::getCudaTemplateVelBuffer() {
        return m_device_cuda_T_vel;
    }

    unsigned ParticleEmitter::getTemplatePartNum() const {
        return m_template_parts_num;
    }

}