//
// Created by ADMIN on 2024/3/30.
//
#include "so_render.hpp"

namespace SoSim {

    SRenderer::SRenderer(std::shared_ptr<ObjectManager> object_manager) {
        m_object_manager_ref = object_manager;

        prepare_render_objs();
    }

    void SRenderer::tick() {

    }

    void SRenderer::prepare_render_objs() {
        auto objs = m_object_manager_ref->getObjects();
        for (const auto &obj: objs) {

        }
    }


}
