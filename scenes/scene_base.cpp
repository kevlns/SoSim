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
    auto cube = objectManager.createObject();
    cube->rename("cube_1");
    auto cubeConfig = new ParticleObjectConfig;
    cubeConfig->particle_radius = 0.05;
    cubeConfig->particle_mat = FLUID;
    cubeConfig->shape = "cube";
    cubeConfig->lb = {-2, -2, -2};
    cubeConfig->size = {4, 4, 4};
    cube->setConfig(cubeConfig);
    cube->update();

    auto box = objectManager.createObject();
    box->rename("box_1");
    auto boxConfig = new ParticleObjectConfig;
    boxConfig->particle_radius = 0.05;
    boxConfig->particle_mat = FIXED_BOUND;
    boxConfig->shape = "box";
    boxConfig->lb = {-5, -5, -5};
    boxConfig->size = {10, 10, 10};
    boxConfig->layer = 1;
    box->setConfig(boxConfig);
    box->update();


    /**  =============================================================
     * create solvers
     */
    SolverManager solverManager;
    auto dfsph_solver = solverManager.createSolver<DFSPHSolver>();
    auto dfsphSolverConfig = new DFSPHSolverConfig;
    dfsphSolverConfig->dt = 0.0015;
    dfsphSolverConfig->rest_density = 1200;
    dfsphSolverConfig->rest_vis = 0.001;
    reinterpret_cast<DFSPHSolver *>(dfsph_solver)->setConfig(dfsphSolverConfig);
    dfsph_solver->attachObject(cube);
//    dfsph_solver->attachObject(box);
    dfsph_solver->initialize();


    /**  =============================================================
     * run simulation
     */
    dfsph_solver->run(10);


    /**  =============================================================
     * terminate
     */
    solverManager.destroy();
    objectManager.destroy();
}