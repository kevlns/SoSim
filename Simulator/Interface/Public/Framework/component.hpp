//@author        : Long Shen
//@date          : 2023/10/27
//@description   :
//@version       : 1.0

#ifndef SOSIM_COMPONENT_HPP
#define SOSIM_COMPONENT_HPP

#include <vector_types.h>
#include <iostream>
#include <cuda_runtime.h>

#include "enum.hpp"

namespace SoSim {
    struct Component {
        virtual void destroy() = 0;

        virtual ~Component() {};
    };

    struct BaseMoveComponent : Component {
        float3 pos;
        float3 acc;
        float3 vel;

        ~BaseMoveComponent() override = default;

        void destroy() override {
            std::cout << "BaseMoveComponent destructed.\n";
        }
    };

    struct BaseParticleComponent : Component {
        float3 *pos_cuPtr;
        std::vector<float3> pos;
        std::vector<float3> vel;
        std::vector<float3> acc;

        ~BaseParticleComponent() override = default;

        void destroy() override {
            if (!pos.empty())
                cudaFree(pos_cuPtr);

            std::cout << "BaseParticleComponent destructed.\n";
        }
    };

    struct Particle_GLRenderableComponent : public Component {
        // use opengl-cuda interop
        // note: if object isn't attached by any solver, then assign cuPtr itself; if later be attached, then free the cuPtr and assign solver pos ptr.
        float3 *render_pos_cuPtr{nullptr};
        float3 *render_color_cuPtr{nullptr};
        float render_part_radius{0};
        unsigned render_part_num{0};
        unsigned VAO{0};
        unsigned VBO{0};
        unsigned EBO{0};

        ~Particle_GLRenderableComponent() override = default;

        void destroy() override {
            render_pos_cuPtr = nullptr;
            render_color_cuPtr = nullptr;
            std::cout << "Particle_GLRenderableComponent destructed.\n";
        }
    };


    static inline Component *createComponentInstance(const ComponentType &component) {
        Component *componentInstance = nullptr;
        switch (component) {
            case ComponentType::BASE_MOVE_COMPONENT:
                componentInstance = new BaseMoveComponent;
                break;
            case ComponentType::BASE_PARTICLE_COMPONENT:
                componentInstance = new BaseParticleComponent;
                break;
            case ComponentType::PARTICLE_GL_RRENDERABLE_COMPONENT:
                componentInstance = new Particle_GLRenderableComponent;
                break;
                // TODO add more Component below
                // case ...
            default:
                // TODO log error: no such component type
                break;
        }

        return componentInstance;
    }

}

#endif //SOSIM_COMPONENT_HPP
