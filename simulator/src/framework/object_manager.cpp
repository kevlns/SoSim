//
// Created by ADMIN on 2024/3/5.
//
#include "framework/object_manager.hpp"

namespace SoSim {
    Object *ObjectManager::createObject() {
        auto obj = new Object(m_hash_counter++);
        m_all_objects.insert(obj);
        return obj;
    }

    void ObjectManager::removeObject(Object *object) {
        if(m_all_objects.count(object) > 0) {
            object->destroy();
            m_all_objects.erase(object);
            delete object;
        }
    }

    void ObjectManager::destroy() {
        for (auto obj: m_all_objects) {
            obj->destroy();
            delete obj;
        }
        m_all_objects.clear();
    }
}