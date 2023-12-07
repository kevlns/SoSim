//@author        : Long Shen
//@date          : 2023/11/27
//@description   :
//@version       : 1.0

#ifndef SOSIM_WIDGET_TOOL_HPP
#define SOSIM_WIDGET_TOOL_HPP

#include <string>
#include <vector_types.h>

namespace SoSim {

    extern void setNextItemMid(const char *text);

    extern const char *makeItemTitle(const char *title, const std::string &postfix);

    extern bool checkFloat3ItemChanged(const char *title, const std::string &itemPostfix, float3 &var);

    extern bool checkFloatItemChanged(const char *title, const std::string &itemPostfix, float &var);

    extern bool checkComboItemChanged(const char *title, const std::string &itemPostfix, int &var, char *map);

    extern bool checkColorEdit3ItemChanged(const char *title, const std::string &itemPostfix, float3 &var);

}

#endif //SOSIM_WIDGET_TOOL_HPP
