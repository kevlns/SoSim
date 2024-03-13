//
// Created by ADMIN on 2024/3/5.
//

#ifndef SOSIM_OBJECT_MANADER_HPP
#define SOSIM_OBJECT_MANADER_HPP

#include <set>

#include "framework/object.hpp"

namespace SoSim {
    class ObjectManager {
    public:
        ~ObjectManager();

        std::shared_ptr<Object> createObject();

        void removeObject(std::shared_ptr<Object> object);

    private:
        void destroy();

    private:
        std::set<std::shared_ptr<Object>> m_objects;
        unsigned m_hash_counter{0};
    };
}

#endif //SOSIM_OBJECT_MANADER_HPP
