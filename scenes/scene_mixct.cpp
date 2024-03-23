//
// Created by ADMIN on 2024/3/22.
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
    cubeConfig_1->particle_mat = JL21CT_NONNEWTON;
    cubeConfig_1->phases.assign({0.6, 0.4});
    cubeConfig_1->shape = "cube";
    cubeConfig_1->lb = {-1, 0.2, -1};
    cubeConfig_1->size = {2, 2, 2};
    cube_1->setName("cube_1");
    cube_1->update();

    auto plane_1 = objectManager.createObject();
    auto planeConfig = plane_1->getParticleObjectConfig();
    planeConfig->particle_radius = 0.05;
    planeConfig->particle_mat = FIXED_BOUND;
    planeConfig->phases.assign({0, 0});
    planeConfig->shape = "plane";
    planeConfig->lb = {-2, -0.4, -2};
    planeConfig->size = {4, 2, 4};
    planeConfig->layer = 2;
    plane_1->setName("plane_1");
    plane_1->update();


    /**  =============================================================
     * create solvers
     */
    SolverManager solverManager;
    auto solver = solverManager.createSolver<JL21CTSolver>();
    auto solver_config = dynamic_cast<JL21CTSolverConfig *>(solver->getConfig().get());
    solver_config->dt = 0.0015;
    solver_config->gravity = {0, -9.8, 0};
    solver_config->rest_density = {995, 1000};
    solver_config->rest_vis = 0.001;
    solver_config->wc_stiff = 1000;
    solver_config->max_neighborNum = 40;

    solver_config->phase1_color = {254, 254, 254};
    solver_config->phase2_color = {255, 2, 2};
    solver_config->rest_density = {995, 1000};
    solver_config->kd = 0.01;
    solver_config->Cd = 0.0005;
    solver_config->relation_time = 0.005;
    solver_config->compliance_0 = 2.5;
    solver_config->compliance_inf = 0.001;
    solver_config->alpha_up_threshold = 0.6;
    solver_config->alpha_low_threshold = 0.001;

    solver_config->scene_lb = {-15, -15, -15};
    solver_config->scene_size = {30, 30, 30};
    solver_config->export_data = true;
    solver->attachObject(cube_1);
//    solver->attachObject(plane_1);

    /**  =============================================================
     * run simulation
     */
    solver->run(2);

}
