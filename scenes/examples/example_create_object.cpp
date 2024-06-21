//
// Created by ADMIN on 2024/6/12.
//

#include "framework/object_manager.hpp"

using namespace SoSim;

int main() {

    /**  =============================================================
     * create objects
     */
    ObjectManager objectManager;

    // create host cube object; lb:left-bottom-coord; size:cube-size
    auto cube = objectManager.createObject();
    auto cube_config = cube->getParticleObjectConfig();
    cube_config->shape = ObjectShape::Cube;
    cube_config->particle_radius = 0.05;
    cube_config->particle_mat = COMMON_NEWTON;
    cube_config->lb = {-1, 0, 0};
    cube_config->size = {2, 2, 2};
    cube_config->vel_start = {0.5, 0, 0}; // if not set, then {0,0,0}
    cube->setName("cube");
    cube->update();

    // create host box object; layer:particle layer, one-layer box or two-layer box
    auto box = objectManager.createObject();
    auto box_config = box->getParticleObjectConfig();
    box_config->shape = ObjectShape::Box;
    box_config->particle_radius = 0.05;
    box_config->particle_mat = FIXED_BOUND;
    box_config->lb = {-1, 0, 0};
    box_config->size = {2, 2, 2};
    box_config->layer = 1;
    box->setName("box");
    box->update();

    // create host plane object; layer:particle layer, one-layer plane or two-layer plane
    auto plane = objectManager.createObject();
    auto plane_config = plane->getParticleObjectConfig();
    plane_config->shape = ObjectShape::Plane;
    plane_config->particle_radius = 0.05;
    plane_config->particle_mat = FIXED_BOUND;
    plane_config->lb = {-5, -2.5, -5};
    plane_config->size = {10, 10, 10};
    plane_config->layer = 2;
    plane->setName("plane");
    plane->update();

    // create host sphere object;
    auto sphere = objectManager.createObject();
    auto sphere_config = sphere->getParticleObjectConfig();
    sphere_config->shape = ObjectShape::Sphere;
    sphere_config->particle_radius = 0.05;
    sphere_config->particle_mat = FIXED_BOUND;
    sphere_config->center = {0, 0, 0};
    sphere_config->volume_radius = 5;
    sphere->setName("sphere");
    sphere->update();

    // create host cylinder object;
    auto cylinder = objectManager.createObject();
    auto cylinder_config = cylinder->getParticleObjectConfig();
    cylinder_config->shape = ObjectShape::Cylinder;
    cylinder_config->particle_radius = 0.05;
    cylinder_config->particle_mat = FIXED_BOUND;
    cylinder_config->center = {0, 0, 0};
    cylinder_config->bottom_area_radius = 2.5;
    cylinder_config->height = 5;
    cylinder->setName("cylinder");
    cylinder->update();

    // sample particles from obj mesh
    auto model = objectManager.createObject();
    auto model_config = model->getParticleObjectConfig();
    model_config->model_file = "D:\\xuyuhang\\obj\\bunny_10k_sample.obj"; // only ply model supported yet.
    model_config->shape = Surface_Sample;
    model_config->ratio = 1.9;
    model_config->particle_radius = 0.025;
    model_config->particle_mat = COMMON_NEWTON;
    model->setName("sample model");
    model->update();
    model->exportAsPly("D:\\xuyuhang\\sample_ply", "bunny_10k_sample");

    // once you use a multiphase solver, you need to assign phase fraction
    // config->phases.assign({0.5,0.5}); only two-phase fluid supported yet.
}