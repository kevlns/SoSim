//@author        : Long Shen
//@date          : 2023/11/21
//@description   :
//@version       : 1.0

#ifndef SOSIM_GUI_HPP
#define SOSIM_GUI_HPP

#include "Public/GUI/camera.hpp"
#include "Public/GUI/shader.hpp"
#include "Public/GUI/light.hpp"
#include "Public/Framework/simulator.hpp"
#include "Public/Framework/scene.hpp"
#include "Public/Framework/solver.hpp"

#include "Private/GUI/gui_header_prefix.hpp"

namespace SoSim::GUI {

    class GUI {
    public:
        GUI();

        ~GUI();

        void run();

        void terminate();

    private:
        bool m_doRendering{true};
        Camera *m_camera{nullptr};

    private:
        // opengl context
        GLFWwindow *m_window{nullptr};
    };

}

#endif //SOSIM_GUI_HPP
