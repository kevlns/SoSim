//@author        : Long Shen
//@date          : 2023/11/21
//@description   :
//@version       : 1.0

#ifndef SOSIM_GUI_HPP
#define SOSIM_GUI_HPP

#include "Public/Framework/simulator.hpp"
#include "Public/GUI/render.hpp"
#include "Private/GUI/gui_header_prefix.hpp"


#ifdef GUI_KEEP_SILENT
const bool KEEP_SILENT = true;
#else
const bool KEEP_SILENT = false;
#endif

namespace SoSim {

    class GUI {
    public:
        GUI();

        ~GUI() = default;

        void run(bool keepSilent = KEEP_SILENT);

        void terminate();

    private:
        void initialize();

        void runPure();

        void runGUI();

    private:
        bool m_keep_silent{false};
        Simulator *m_simulator{nullptr};
        Renderer *m_renderer{nullptr};

    private:
        // gui resource
        ImGuiIO m_io;
        GLFWwindow *m_main_window{nullptr};

    };

}

#endif //SOSIM_GUI_HPP
