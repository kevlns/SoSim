//@author        : Long Shen
//@date          : 2023/10/11
//@description   : 
//@version       : 1.0

#ifndef SOSIM_OBJECT_HPP
#define SOSIM_OBJECT_HPP

#include <vector>
#include <iostream>
#include <unordered_map>

#include "Public/Framework/component.hpp"
#include "Public/Framework/framework_config.hpp"

namespace SoSim {
    /**
     * @breif TODO
     */
    class Object {
    public:
        Object();

        ~Object();

        void addComponent(const ComponentType &component);

        void removeComponent(const ComponentType &component);

        bool hasComponent(const ComponentType &component) const;

        void destroy();

        void refresh();

        Component *getComponent(const ComponentType &component);

        ObjectConfig *getConfig();

    private:

        ObjectConfig *m_config;

        std::unordered_map<ComponentType, Component *> m_components;
    };
}


#endif //SOSIM_OBJECT_HPP
