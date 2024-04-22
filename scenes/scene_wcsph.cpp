//
// Created by ADMIN on 2024/3/14.
//

#include "framework/object_manager.hpp"
#include "framework/solver_manager.hpp"

using namespace SoSim;

int main() {

    /**  =============================================================
     * create objects
     */
    ObjectManager objectManager;
    auto cube_1 = objectManager.createObject();
    auto cubeConfig_1 = cube_1->getParticleObjectConfig();
    cubeConfig_1->particle_radius = 0.05;
    cubeConfig_1->particle_mat = COMMON_FLUID;
    cubeConfig_1->shape = ObjectShape::Cube;
    cubeConfig_1->lb = {-1, 0.2, -1};
    cubeConfig_1->size = {2, 2, 2};
    cube_1->setName("cube_1");
    cube_1->update();

    auto cube_2 = objectManager.createObject();
    auto cubeConfig_2 = cube_2->getParticleObjectConfig();
    cubeConfig_2->particle_radius = 0.1;
    cubeConfig_2->particle_mat = COMMON_FLUID;
    cubeConfig_2->shape = ObjectShape::Cube;
    cubeConfig_2->lb = {1.5, 0.1, 0};
    cubeConfig_2->size = {2, 2, 2};
    cube_2->setName("cube_2");
    cube_2->update();

    auto plane_1 = objectManager.createObject();
    auto planeConfig = plane_1->getParticleObjectConfig();
    planeConfig->particle_radius = 0.05;
    planeConfig->particle_mat = FIXED_BOUND;
    planeConfig->shape = ObjectShape::Plane;
    planeConfig->lb = {-2, -0.4, -2};
    planeConfig->size = {4, 2, 4};
    planeConfig->layer = 2;
    plane_1->setName("plane_1");
    plane_1->update();


    /**  =============================================================
     * create solvers
     */
    SolverManager solverManager;
    auto wcsph_solver = solverManager.createSolver<WCSPHSolver>();
    auto wcsphSolverConfig = dynamic_cast<WCSPHSolverConfig *>(wcsph_solver->getConfig().get());
    wcsphSolverConfig->dt = 0.001;
    wcsphSolverConfig->gravity = {0, -9.8, 0};
    wcsphSolverConfig->rest_density = 1000;
    wcsphSolverConfig->rest_rigid_density = 1000;
    wcsphSolverConfig->cs = 5;
    wcsphSolverConfig->rest_vis = 0.001;
    wcsphSolverConfig->max_neighborNum = 40;
    wcsphSolverConfig->export_data = true;
    wcsphSolverConfig->scene_lb = {-15, -15, -15};
    wcsphSolverConfig->scene_size = {30, 30, 30};
    wcsph_solver->attachObject(cube_1);
//    wcsph_solver->attachObject(cube_2);
    wcsph_solver->attachObject(plane_1);

    /**  =============================================================
     * run simulation
     */
    wcsph_solver->run(2);

}
