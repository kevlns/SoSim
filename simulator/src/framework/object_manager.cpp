//
// Created by ADMIN on 2024/3/5.
//
#include <memory>

#include "framework/object_manager.hpp"

namespace SoSim {

    ObjectManager::~ObjectManager() {
        destroy();
    }

    std::shared_ptr<Object> ObjectManager::createObject() {
        auto obj = std::make_shared<Object>(m_hash_counter++);
        m_objects.insert(obj);
        return obj;
    }

    void ObjectManager::removeObject(std::shared_ptr<Object> object) {
        if (m_objects.count(object) > 0) {
            m_objects.erase(object);
            object.reset();
        }
    }

    void ObjectManager::destroy() {
        for (auto obj: m_objects) {
            obj.reset();
        }
        m_objects.clear();
    }

    std::set<std::shared_ptr<Object>> &ObjectManager::getObjects() {
        return m_objects;
    }

}