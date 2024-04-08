//
// Created by ADMIN on 2024/3/30.
//

#ifndef SOSIM_SO_RENDER_HPP
#define SOSIM_SO_RENDER_HPP

#include "render_component.hpp"

#include "core/math/matrix.hpp"
#include "framework/object_manager.hpp"

namespace SoSim {
    class SRenderer {
    public:
        explicit SRenderer(std::shared_ptr<ObjectManager> object_manager);

        void tick();

    private:
        void prepare_render_objs();

        void particle_obj_render_system();

        void mesh_obj_render_system();

    private:
        std::shared_ptr<ObjectManager> m_object_manager_ref;
        std::vector<std::shared_ptr<Object>> m_particle_render_targets;
        std::vector<std::shared_ptr<Object>> m_mesh_render_targets;

    };
}

#endif //SOSIM_SO_RENDER_HPP
