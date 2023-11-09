//@author        : Long Shen
//@date          : 2023/10/11
//@description   : 
//@version       : 1.0

#ifndef SOSIM_OBJECT_HPP
#define SOSIM_OBJECT_HPP

#include <vector>
#include <iostream>

#include "Public/Framework/component.hpp"

namespace SoSim {

    /**
     * @breif TODO
     */
    class Object {
    public:
        Object() = default;

        ~Object() = default;

        void attachComponent();

        void detachComponent();

        void destroy();

    public:
        uint32_t m_id;

    };

}


#endif //SOSIM_OBJECT_HPP
