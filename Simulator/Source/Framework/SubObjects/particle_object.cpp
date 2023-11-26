//@author        : Long Shen
//@date          : 2023/11/22
//@description   :
//@version       : 1.0

#include "Public/Framework/SubObjects/particle_object.hpp"
#include "Public/Shared/ModelUtils/model_tool.hpp"

namespace SoSim {

    ParticleObject::ParticleObject(ParticleObjectConfig &config) {
        m_config = std::move(config);
        gen();
    }

    void ParticleObject::gen() {
        if (!m_config.model_file_path.empty()) {

        } else {
            if (m_config.shape == "cube") {
                m_pos = generate_cube(m_config.lb, m_config.size, m_config.particle_radius);
            } else if (m_config.shape == "box") {
                m_pos = generate_box(m_config.lb, m_config.size, m_config.particle_radius);
            } else if (m_config.shape == "plane-x") {
                m_pos = generate_plane_X(m_config.lb, m_config.size, m_config.particle_radius);
            } else if (m_config.shape == "plane-z") {
                m_pos = generate_plane_Z(m_config.lb, m_config.size, m_config.particle_radius);
            } else if (m_config.shape == "cylinder") {
                m_pos = generate_cylinder(m_config.top_center, m_config.height, m_config.area_radius,
                                          m_config.particle_radius);
            }
        }
    }

    void ParticleObject::refresh() {
        gen();
    }

    void ParticleObject::destroy() {

    }

}