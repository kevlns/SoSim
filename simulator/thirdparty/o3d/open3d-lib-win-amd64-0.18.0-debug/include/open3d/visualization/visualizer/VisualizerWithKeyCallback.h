// ----------------------------------------------------------------------------
// -                        Open3D: www.open3d.org                            -
// ----------------------------------------------------------------------------
// Copyright (c) 2018-2023 www.open3d.org
// SPDX-License-Identifier: MIT
// ----------------------------------------------------------------------------

#pragma once

#include <map>

#include "open3d/visualization/visualizer/Visualizer.h"

namespace open3d {
namespace visualization {

/// \class VisualizerWithKeyCallback
///
/// \brief Visualizer with custom key callback capabilities.
class VisualizerWithKeyCallback : public Visualizer {
public:
    typedef std::pair<int, std::function<bool(Visualizer *)>> KeyCallbackPair;

public:
    /// \brief Default Constructor.
    VisualizerWithKeyCallback();
    ~VisualizerWithKeyCallback() override;
    VisualizerWithKeyCallback(const VisualizerWithKeyCallback &) = delete;
    VisualizerWithKeyCallback &operator=(const VisualizerWithKeyCallback &) =
            delete;

public:
    void PrintVisualizerHelp() override;
    void RegisterKeyCallback(int key,
                             std::function<bool(Visualizer *)> callback);
    /// Register callback function with access to GLFW key actions.
    ///
    /// \param key GLFW key value, see [GLFW key
    /// values](https://www.glfw.org/docs/latest/group__keys.html).
    ///
    /// \param callback The callback function. The callback function takes
    /// `Visualizer *`, `action` and `mods` as input and returns a boolean
    /// indicating UpdateGeometry() needs to be run. The `action` can be one of
    /// GLFW_RELEASE (0), GLFW_PRESS (1) or GLFW_REPEAT (2), see [GLFW input
    /// interface](https://www.glfw.org/docs/latest/group__input.html). The
    /// `mods` specifies the modifier key, see [GLFW modifier
    /// key](https://www.glfw.org/docs/latest/group__mods.html).
    void RegisterKeyActionCallback(
            int key, std::function<bool(Visualizer *, int, int)> callback);

protected:
    void KeyPressCallback(GLFWwindow *window,
                          int key,
                          int scancode,
                          int action,
                          int mods) override;
    std::string PrintKeyToString(int key);

protected:
    std::map<int, std::function<bool(Visualizer *)>> key_to_callback_;
    std::map<int, std::function<bool(Visualizer *, int, int)>>
            key_action_to_callback_;
};

}  // namespace visualization
}  // namespace open3d
