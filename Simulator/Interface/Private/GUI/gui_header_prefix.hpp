//@author        : Long Shen
//@date          : 2023/11/22
//@description   :
//@version       : 1.0

#ifndef SOSIM_GUI_HEADER_PREFIX_HPP
#define SOSIM_GUI_HEADER_PREFIX_HPP

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "GLFW/glfw3.h"

#if defined(_MSC_VER) && (_MSC_VER >= 1900) && !defined(IMGUI_DISABLE_WIN32_FUNCTIONS)
#pragma comment(lib, "legacy_stdio_definitions")
#endif

// This example can also compile and run with Emscripten! See 'Makefile.emscripten' for details.
#ifdef __EMSCRIPTEN__
#include "../libs/emscripten/emscripten_mainloop_stub.h"
#endif

#endif //SOSIM_GUI_HEADER_PREFIX_HPP
