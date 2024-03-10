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
        Object *createObject();

        void removeObject(Object *object);

        void destroy();

    private:
        std::set<Object *> m_all_objects;
        unsigned m_hash_counter{0};
    };
}

#endif //SOSIM_OBJECT_MANADER_HPP
