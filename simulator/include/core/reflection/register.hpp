//
// Created by ADMIN on 2024/3/17.
//

#ifndef SOSIM_REGISTER_HPP
#define SOSIM_REGISTER_HPP

#include "reflection.hpp"

namespace SoSim {

    class SoRegister {
    public:

        void classRegister(std::string& class_name);

        void funcRegister();

        void dataRegister();

    };

}

#endif //SOSIM_REGISTER_HPP
