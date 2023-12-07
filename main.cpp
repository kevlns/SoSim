//@author        : Long Shen
//@date          : 2023/11/21
//@description   :
//@version       : 1.0

//#define GUI_KEEP_SILENT
#include "Public/GUI/gui.hpp"

using namespace SoSim;

int main() {
    GUI sosim_editor;

    sosim_editor.run();

    sosim_editor.terminate();
}