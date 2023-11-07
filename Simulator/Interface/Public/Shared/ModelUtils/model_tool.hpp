//@author        : Long Shen
//@date          : 2023/10/26
//@description   :
//@version       : 1.0

#ifndef SOSIM_MODEL_TOOL_HPP
#define SOSIM_MODEL_TOOL_HPP

#include <vector>
#include <vector_types.h>
#include <fstream>
#include <iostream>

#include "Public/Framework/mat.hpp"

namespace SoSim {

    extern uint32_t loadObjFromFile(const std::string &file_path, std::vector<float3> &pos);

    extern void genObjFromJson(const std::string &json_path, std::vector<float3> &pos, std::vector<float3> &vel,
                               std::vector<float> &den, std::vector<Material> &mat, std::vector<Phase> &phase,
                               float &radius);

    extern std::vector<float3>
    generate_cube(float3 cubeLB, float3 cubeSize, float particleRadius);

    extern std::vector<float3>
    generate_box(float3 boxLB, float3 boxSize, float particleRadius);

    extern std::vector<float3>
    generate_plane_X(float3 planeLB, float3 planeSize, float particleRadius);

    extern std::vector<float3>
    generate_plane_Z(float3 planeLB, float3 planeSize, float particleRadius);

    extern std::vector<float3>
    generate_cylinder(float3 topCenter, float height, float areaRadius, float particleRadius);

    extern std::pair<std::vector<float3>, float>
    import_ply_model(const std::string &file_name);

    extern void write_ply(const std::string &filename, const std::vector<float3> &points);

}

#endif //SOSIM_MODEL_TOOL_HPP
