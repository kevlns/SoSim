//@author        : Long Shen
//@date          : 2023/11/27
//@description   :
//@version       : 1.0

#include "Private/GUI/widget_tool.hpp"
#include "Public/Shared/Math/helper_math.hpp"
#include "imgui.h"

namespace SoSim {

    extern void setNextItemMid(const char *text) {
        ImVec2 textSize = ImGui::CalcTextSize(text);
        float windowWidth = ImGui::GetContentRegionAvail().x;
        float buttonWidth = textSize.x + ImGui::GetStyle().FramePadding.x * 2.0f;
        float centerPos = (windowWidth - buttonWidth) * 0.5f;
        ImGui::SetCursorPosX(centerPos);
    }

    extern const char *makeItemTitle(const char *title, const std::string &postfix) {
        static std::string s;
        std::string s_n = std::string(title) + "##SOSIM_" + postfix;
        s = s_n;
        return s.c_str();
    }

    extern bool checkFloat3ItemChanged(const char *title, const std::string &itemPostfix, float3 &var) {
        ImGui::SetNextItemWidth(250.0f);
        float3 var_old = var;
        ImGui::InputFloat3(makeItemTitle(title, itemPostfix), &var.x);
        if (length(var_old - var) > 1e-6)
            return true;
        return false;
    }

    extern bool checkFloatItemChanged(const char *title, const std::string &itemPostfix, float &var) {
        ImGui::SetNextItemWidth(250.0f);
        bool flag = ImGui::InputFloat(makeItemTitle(title, itemPostfix), &var);
        return flag;
    }

    extern bool checkComboItemChanged(const char *title, const std::string &itemPostfix, int &var, char *map) {
        ImGui::SetNextItemWidth(250.0f);
        bool flag = ImGui::Combo(makeItemTitle(title, itemPostfix), &var, map, IM_ARRAYSIZE(map));
        return flag;
    }

    extern bool checkColorEdit3ItemChanged(const char *title, const std::string &itemPostfix, float3 &var) {
        ImGui::SetNextItemWidth(250.0f);
        bool flag = ImGui::ColorEdit3(makeItemTitle(title, itemPostfix),
                                      reinterpret_cast<float *>(&var));
        return flag;
    }

}

