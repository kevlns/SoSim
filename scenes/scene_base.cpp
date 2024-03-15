//
// Created by ADMIN on 2024/3/7.
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
    cube_1->setName("cube_1");
    cubeConfig_1->particle_radius = 0.05;
    cubeConfig_1->particle_mat = FLUID;
    cubeConfig_1->shape = "cube";
    cubeConfig_1->lb = {-1, 0, 0};
    cubeConfig_1->size = {2, 2, 2};
    cubeConfig_1->vel_start = {0.5, 0, 0};
    cube_1->update();

    auto cube_2 = objectManager.createObject();
    auto cubeConfig_2 = cube_2->getParticleObjectConfig();
    cube_2->setName("cube_2");
    cubeConfig_2->particle_radius = 0.05;
    cubeConfig_2->particle_mat = FLUID;
    cubeConfig_2->shape = "cube";
    cubeConfig_2->lb = {1.5, 0, 0};
    cubeConfig_2->size = {2, 2, 2};
    cube_2->update();

    auto plane = objectManager.createObject();
    auto planeConfig = plane->getParticleObjectConfig();
    plane->setName("plane_1");
    planeConfig->particle_radius = 0.05;
    planeConfig->particle_mat = FIXED_BOUND;
    planeConfig->shape = "plane";
    planeConfig->lb = {-5, -2.5, -5};
    planeConfig->size = {10, 10, 10};
    planeConfig->layer = 2;
    plane->update();


    /**  =============================================================
     * create solvers
     */
    SolverManager solverManager;
    auto dfsph_solver = solverManager.createSolver<DFSPHSolver>();
    auto dfsphSolverConfig = dynamic_cast<DFSPHSolverConfig *>(dfsph_solver->getConfig().get());
    dfsphSolverConfig->dt = 0.01;
    dfsphSolverConfig->gravity = {0,-0.5,0};
    dfsphSolverConfig->rest_density = 1000;
    dfsphSolverConfig->rest_rigid_density = 1000;
    dfsphSolverConfig->rest_vis = 0.001;
    dfsphSolverConfig->max_neighborNum = 35;
    dfsphSolverConfig->export_data = false;
    dfsph_solver->attachObject(cube_1);
    dfsph_solver->attachObject(cube_2);
//    dfsph_solver->attachObject(plane);

    /**  =============================================================
     * run simulation
     */
    dfsph_solver->run(2);

}