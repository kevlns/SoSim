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

namespace SoSim {

    extern void genObjFromJson();

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

}

#endif //SOSIM_MODEL_TOOL_HPP
