//
// Created by ADMIN on 2024/3/26.
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
    cubeConfig_1->particle_radius = 0.075;
    cubeConfig_1->particle_mat = IMSCT_NONNEWTON;
    cubeConfig_1->phases.assign({0.5, 0.5});
    cubeConfig_1->shape = "cube";
    cubeConfig_1->lb = {-1, -1.2, -1};
    cubeConfig_1->size = {2, 4, 2};
    cube_1->setName("cube_1");
    cube_1->update();

    auto plane_1 = objectManager.createObject();
    auto planeConfig = plane_1->getParticleObjectConfig();
    planeConfig->particle_radius = 0.075;
    planeConfig->particle_mat = FIXED_BOUND;
    planeConfig->phases.assign({0, 0});
    planeConfig->shape = "plane";
    planeConfig->lb = {-3, -2, -3};
    planeConfig->size = {6, 2, 6};
    planeConfig->layer = 2;
    plane_1->setName("plane_1");
    plane_1->update();

    auto box_1 = objectManager.createObject();
    auto boxConfig = box_1->getParticleObjectConfig();
    boxConfig->particle_radius = 0.075;
    boxConfig->particle_mat = FIXED_BOUND;
    boxConfig->phases.assign({0, 0});
    boxConfig->shape = "box";
    boxConfig->lb = {-2, -2, -2};
    boxConfig->size = {4, 6, 4};
    boxConfig->layer = 2;
    box_1->setName("box_1");
    box_1->update();


    /**  =============================================================
     * create solvers
     */
    SolverManager solverManager;
    auto solver = solverManager.createSolver<IMSCTSolver>();
    auto solver_config = dynamic_cast<IMSCTSolverConfig *>(solver->getConfig().get());
    solver_config->dt = 0.002;
    solver_config->gravity = {0, -9.8, 0};
    solver_config->scene_lb = {-30, -30, -30};
    solver_config->scene_size = {60, 60, 60};
    solver_config->export_data = true;

    solver_config->max_neighborNum = 60;
    solver_config->rest_viscosity = 0.005;
    solver_config->rest_density = {900, 1000};
    solver_config->rest_rigid_density = 1000;
    solver_config->rest_bound_density = 1000;
    solver_config->div_free_threshold = 1e-4;
    solver_config->incompressible_threshold = 1e-4;
    solver_config->Cf = 0.01;
    solver_config->Cd = 1;
    solver_config->phase1_color = {50, 0, 200};
    solver_config->phase2_color = {200, 0, 50};

    solver_config->Cd0 = 0.4;
    solver_config->ct_thinning_exp0 = 0.6;
    solver_config->ct_relaxation_time = 0.01;
    solver_config->solution_vis0 = 2500;

    solver->attachObject(cube_1);
    solver->attachObject(plane_1);
//    solver->attachObject(box_1);

    /**  =============================================================
     * run simulation
     */
    solver->run(2);

}
