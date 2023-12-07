//@author        : Long Shen
//@date          : 2023/11/26
//@description   :
//@version       : 1.0

#include <string>

#include "imgui.h"
#include "Private/GUI/widget_tool.hpp"
#include "Private/GUI/widgets.hpp"
#include "Public/Shared/Math/helper_math.hpp"
#include "Public/PhysicalSolvers/solver_ref.hpp"

namespace SoSim {
    extern void ShowMainWindow(Simulator *simulator) {
        ImVec4 *colors = ImGui::GetStyle().Colors;
        // 设置新的按钮颜色
        colors[ImGuiCol_Button] = ImVec4(0.0f, 0.0f, 0.0f, 1.0f);
        colors[ImGuiCol_ButtonHovered] = ImVec4(0.1f, 0.1f, 0.1f, 1.0f);
        colors[ImGuiCol_ButtonActive] = ImVec4(0.5f, 0.5f, 0.5f, 1.0f); // 深红色

        ImGui::Begin("SoSim Simulator");

        // 推送新的padding值到样式堆栈
        const ImVec2 customPadding = ImVec2(25.0f, 25.0f);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, customPadding);

        ImGui::SetNextItemAllowOverlap();
        ImGui::Text("Add Scene");;
        ImGui::SameLine();
        if (ImGui::SmallButton(makeItemTitle("+", "sim_setup")))
            simulator->addSceneDefault();

        ShowScenesConfig(simulator);

        ImGui::PopStyleVar();
        ImGui::End();
    }

    extern void ShowScenesConfig(Simulator *simulator) {
        if (simulator->getSceneMap().empty())
            return;

        if (ImGui::TreeNode("Scene List")) {
            int id = 1;
            for (auto &scene: simulator->getSceneMap()) {
                std::string itemPostfix = "scene_config_" + std::to_string(id);

                ImGui::SetNextItemAllowOverlap();
                bool showSceneDetail = ImGui::TreeNode((void *) (intptr_t) id, "SCENE_ %d", id);
                ImGui::SameLine();
                if (ImGui::SmallButton(makeItemTitle("-", itemPostfix)))
                    simulator->addToSceneRemoveList(scene.first);

                if (showSceneDetail) {

                    if (ImGui::Checkbox(makeItemTitle("preview", itemPostfix), &scene.first->preview)) {
                        if (scene.first->preview) {
                            // TODO 这里应该是跟 render 相关的设置而不是 scene->refresh()
                        }
                    }

                    ImGui::SeparatorText("Scene Config");

                    if (checkFloat3ItemChanged("scene lb", itemPostfix, scene.first->scene_lb))
                        scene.second->refresh();

                    if (checkFloat3ItemChanged("scene size", itemPostfix, scene.first->scene_size))
                        scene.second->refresh();

                    ImGui::SeparatorText("");
                    ImGui::SetNextItemAllowOverlap();
                    ImGui::Text("Add Object");
                    ImGui::SameLine();
                    if (ImGui::SmallButton(makeItemTitle("+", itemPostfix + "_1")))
                        scene.second->addObjectDefault();
                    ShowObjectDetail(scene.second);

//                    ImGui::SeparatorText("");
//                    ImGui::SetNextItemAllowOverlap();
//                    ImGui::Text("Add Solver");;
//                    ImGui::SameLine();
//                    if (ImGui::SmallButton(makeItemTitle("+", itemPostfix + "_2")))
//                        scene->addSolverDefault();
//                    ShowSolverDetail(scene);

                    ImGui::TreePop();
                }
                id++;
            }
            simulator->clearSceneRemoveList();
            ImGui::TreePop();
        }
    }

    extern void ShowObjectDetail(Scene *scene) {
        static const char *shapeItems[] = {"None", "cube", "box", "cylinder", "plane-x", "plane-z"};
        using EnumType = std::underlying_type<ComponentType>::type;
        static auto minEnum = static_cast<EnumType>(ComponentType::COMPONENT_NONE);
        static auto maxEnum = static_cast<EnumType>(ComponentType::COMPONENT_BOTTOM);
        static int numberOfComponentCases = maxEnum - minEnum;

        if (scene->getObjectMap().empty())
            return;

        if (ImGui::TreeNode("Object List")) {

            int id = 1;
            for (auto obj: scene->getObjectMap()) {
                std::string obj_name = "OBJ_" + std::to_string(id) + ": " + matItems[obj.first->mat_index];
                std::string itemPostfix = "obj_config_" + std::to_string(id);

                ImGui::SetNextItemAllowOverlap();
                bool showObjectDetail = ImGui::TreeNode((void *) (intptr_t) id, obj_name.c_str(), id);
                ImGui::SameLine();
                if (ImGui::SmallButton(makeItemTitle("-", itemPostfix)))
                    scene->addToObjectRemoveList(obj.first);

                if (showObjectDetail) {

                    // for transformer
                    if (checkFloat3ItemChanged("transform", itemPostfix, obj.first->t_config.transform))
                        obj.second->refresh();

                    if (checkFloat3ItemChanged("scale", itemPostfix, obj.first->t_config.scale)) {
                        if (obj.first->t_config.scale.x < 0)
                            obj.first->t_config.scale.x = 1.0;
                        if (obj.first->t_config.scale.y < 0)
                            obj.first->t_config.scale.y = 1.0;
                        if (obj.first->t_config.scale.z < 0)
                            obj.first->t_config.scale.z = 1.0;
                        obj.second->refresh();
                    }

                    if (checkFloat3ItemChanged("rotate", itemPostfix, obj.first->t_config.rotate))
                        obj.second->refresh();

                    // for model_file_path
                    // TODO make a selective resource manager window
                    ImGui::SeparatorText("");
                    ImGui::SetNextItemWidth(250.0f);
                    ImGui::InputText(makeItemTitle("model path", itemPostfix),
                                     obj.first->model_file_path_cStr,
                                     256);


                    // for shape
                    ImGui::SetNextItemWidth(250.0f);
                    if (ImGui::Combo(makeItemTitle("shape", itemPostfix), &obj.first->shape_index, shapeItems,
                                     IM_ARRAYSIZE(shapeItems))) {
                        obj.first->shape = shapeItems[obj.first->shape_index];
                        obj.second->refresh();
                    }

                    // for material
                    ImGui::SetNextItemWidth(250.0f);
                    if (ImGui::Combo(makeItemTitle("material", itemPostfix), &obj.first->mat_index, matItems,
                                     IM_ARRAYSIZE(matItems))) {
                        obj.first->mat = static_cast<MaterialType>(obj.first->mat_index);
                        obj.second->refresh();
                    }

                    // for phase
                    ImGui::SetNextItemWidth(250.0f);
                    if (ImGui::Combo(makeItemTitle("phase", itemPostfix), &obj.first->phase_index, phaseItems,
                                     IM_ARRAYSIZE(phaseItems))) {
                        obj.first->phase = static_cast<PhaseType>(obj.first->phase_index);
                        obj.second->refresh();
                    }

                    // for particle radius
                    ImGui::SetNextItemWidth(250.0f);
                    if (checkFloatItemChanged("particle radius", itemPostfix, obj.first->particle_radius)) {
                        if (obj.first->particle_radius < 0.01)
                            obj.first->particle_radius = 0.01;
                        if (obj.first->particle_radius > 1.0)
                            obj.first->particle_radius = 1.0;
                        obj.second->refresh();
                    }

                    // for color
                    if (checkColorEdit3ItemChanged("color", itemPostfix, obj.first->color))
                        obj.second->refresh();


                    if (obj.first->shape_index == 3) {

                        if (checkFloatItemChanged("height", itemPostfix, obj.first->height))
                            obj.second->refresh();

                        if (checkFloatItemChanged("area radius", itemPostfix, obj.first->area_radius))
                            obj.second->refresh();

                        if (checkFloat3ItemChanged("top center", itemPostfix, obj.first->top_center))
                            obj.second->refresh();

                    } else if (obj.first->shape_index != 0) {

                        if (checkFloat3ItemChanged("lb", itemPostfix, obj.first->lb))
                            obj.second->refresh();

                        if (checkFloat3ItemChanged("size", itemPostfix, obj.first->size))
                            obj.second->refresh();

//                        ImGui::SetNextItemWidth(250.0f);
//                        if (ImGui::InputInt(makeItemTitle("particle layer num", itemPostfix),
//                                            &obj.first->particle_layer_num)) {
//                            if (obj.first->particle_layer_num < 1)
//                                obj.first->particle_layer_num = 1;
//                        }
                    }

                    // for components
                    ImGui::SeparatorText("");
                    if (obj.first->component_selected.empty()) {
                        obj.first->component_selected = std::vector<bool>(numberOfComponentCases, false);
                    }
                    for (int i = 1; i < obj.first->component_selected.size(); ++i) {
                        bool selected = obj.first->component_selected[i];
                        if (ImGui::Checkbox(makeItemTitle(componentItems[i], itemPostfix), &selected)) {
                            auto type = static_cast<ComponentType>(i);
                            if (selected) {
                                if (!obj.second->hasComponent(type))
                                    obj.second->addComponent(type);
                            } else {
                                if (obj.second->hasComponent(type))
                                    obj.second->removeComponent(type);
                            }
                            obj.first->component_selected[i] = selected;
                        }
                    }

                    ImGui::TreePop();
                }
                id++;
            }
            scene->clearObjectRemoveList();
            ImGui::TreePop();
        }
    }
//
//    extern void ShowSolverDetail(Scene *scene) {
//        if (scene->getSolverConfigs().empty())
//            return;
//
//        if (ImGui::TreeNode("Solver List")) {
//
//            int id = 1;
//            for (auto solverConfig: scene->getSolverConfigs()) {
//                std::string solver_name =
//                        "SOLVER_" + std::to_string(id) + ": " + solverItems[solverConfig->solver_index];
//                std::string itemPostfix = "solver_config_" + std::to_string(id);
//
//                ImGui::SetNextItemAllowOverlap();
//                bool showSolverDetail = ImGui::TreeNode((void *) (intptr_t) id, solver_name.c_str(), id);
//                ImGui::SameLine();
//                if (ImGui::SmallButton(makeItemTitle("-", itemPostfix)))
//                    scene->addToSolverRemoveList(solverConfig);
//
//                if (showSolverDetail) {
//
//                    ImGui::SetNextItemWidth(250.0f);
//                    ImGui::InputFloat3(makeItemTitle("gravity", itemPostfix), &solverConfig->gravity.x);
//                    ImGui::SetNextItemWidth(250.0f);
//                    if (ImGui::InputFloat(makeItemTitle("dt", itemPostfix), &solverConfig->dt)) {
//                        if (solverConfig->dt < 0)
//                            solverConfig->dt = 0;
//                    }
//                    ImGui::SetNextItemWidth(250.0f);
//                    ImGui::Text(" -- keep same with obj particle radius -- ");
//                    ImGui::SetNextItemWidth(250.0f);
//                    ImGui::InputFloat(makeItemTitle("particle radius", itemPostfix),
//                                      &solverConfig->unified_particle_radius);
//
//                    // for solver instance type
//                    ImGui::SetNextItemWidth(250.0f);
//                    if (ImGui::Combo(makeItemTitle("solver instance", itemPostfix), &solverConfig->solver_index,
//                                     solverItems,
//                                     IM_ARRAYSIZE(solverItems))) {
//                        solverConfig->instanceType = static_cast<SolverInstanceType>(solverConfig->solver_index);
//                    }
//
//                    ImGui::TreePop();
//                }
//                id++;
//            }
//            scene->clearSolverRemoveList();
//            ImGui::TreePop();
//        }
//    }

}
