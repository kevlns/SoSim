//@author        : Long Shen
//@date          : 2023/10/11
//@description   :
//@version       : 1.0

#include "Public/Framework/object.hpp"
#include "Public/Shared/ModelUtils/model_builder.hpp"

namespace SoSim {

    Object::Object() {
        m_config = new ObjectConfig;

        std::cout << "Object created.\n";
    }

    Object::~Object() {
        delete m_config;
        m_config = nullptr;
    }

    bool Object::hasComponent(const ComponentType &component) const {
        return m_components.count(component) > 0;
    }

    void Object::destroy() {
        for (auto &component: m_components) {
            component.second->destroy();
            delete component.second;
            component.second = nullptr;
        }

        m_components.clear();

        std::cout << "Object destroyed.\n";
    }

    void Object::refresh() {
        std::cout << "Object refresh.\n";

        // TODO generate particles
        // TODO only particle component supported yet!
        if (hasComponent(BASE_PARTICLE_COMPONENT)) {
            BaseParticleComponent *particleComponent = static_cast<BaseParticleComponent *>(getComponent(
                    BASE_PARTICLE_COMPONENT));

            if (particleComponent->pos_cuPtr)
                cudaFree(particleComponent->pos_cuPtr);
            particleComponent->pos = genFromObjectConfig(m_config);
            cudaMalloc((void **) &particleComponent->pos_cuPtr, particleComponent->pos.size() * sizeof(float3));
            cudaMemcpy(particleComponent->pos_cuPtr, particleComponent->pos.data(),
                       particleComponent->pos.size() * sizeof(float3), cudaMemcpyHostToDevice);
        }

        // if has other component
        // ...
    }

    void Object::addComponent(const ComponentType &component) {
        if (m_components.count(component) > 0) {
            // TODO log error: component already exists
            std::cout << "error: component already exists.\n";
            return;
        }

        m_components[component] = createComponentInstance(component);
        m_config->components.insert(component);
        std::cout << "Component created.\n";

        refresh();
    }

    void Object::removeComponent(const ComponentType &component) {
        if (m_components.count(component) > 0) {
            m_components[component]->destroy();
            delete m_components[component];
            m_components[component] = nullptr;
            m_components.erase(component);
            m_config->components.erase(component);

            std::cout << "Component destroyed.\n";

            refresh();
        }
    }

    Component *Object::getComponent(const ComponentType &component) {
        if (hasComponent(component))
            return m_components[component];

        // TODO log errro: component doesn't exist
        std::cout << "error: component doesn't exist.\n";
        return nullptr;
    }

    ObjectConfig *Object::getConfig() {
        return m_config;
    }

}
