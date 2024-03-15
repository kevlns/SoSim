//
// Created by ADMIN on 2024/3/13.
//

#ifndef SOSIM_MODEL_EXPORTER_HPP
#define SOSIM_MODEL_EXPORTER_HPP

#include <vector>
#include <string>

#include "core/math/matrix.hpp"

namespace SoSim {

    class ModelExporter {
    public:
        static void exportVecSetAsPly(const std::string &dir,
                                      const std::string &file_name,
                                      const std::vector<Vec3f> &pos);

        static void exportVecSetAsPly(const std::string &dir,
                                      const std::string &file_name,
                                      const std::vector<Vec3f> &pos,
                                      const std::vector<Vec3f> &color);
    };

}

#endif //SOSIM_MODEL_EXPORTER_HPP
