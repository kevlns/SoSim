//
// Created by ADMIN on 2024/6/12.
//

#include "core/math/matrix_v2.hpp"
#include "core/math/matrix.hpp"

using namespace SoSim;

int main() {

    Mat22<float> m1;
    m1 = {1., 2.,
          3., 4.};

    Vec2<float> v1;
    v1 = {1, 2};

    std::cout << v1[0] << std::endl;
}