//
// Created by ADMIN on 2024/4/10.
//

#include "framework/object_manager.hpp"
#include "framework/solver_manager.hpp"

using namespace SoSim;

int main() {

    /**  =============================================================
     * create objects
     */
    ObjectManager objectManager;
    auto diamond = objectManager.createObject();
    auto di_config = diamond->getParticleObjectConfig();
    di_config->particle_radius = 0.075;
    di_config->particle_mat = IMSCT_NONNEWTON;
    di_config->phases.assign({0.5, 0.5});
    di_config->vel_start = {0, 0, 0};
    di_config->model_file = "F:\\DataSet.Research\\SCA2024\\models\\scene_diamonds\\diamond.ply";
    diamond->setName("diamond");
    diamond->update();

    auto ball = objectManager.createObject();
    auto ball_config = ball->getParticleObjectConfig();
    ball_config->particle_radius = 0.075;
    ball_config->particle_mat = IMSCT_NONNEWTON;
    ball_config->phases.assign({0.3, 0.7});
    ball_config->shape = ObjectShape::Sphere;
    ball_config->center = {0, 5, 0};
    ball_config->volume_radius = 4;
    ball->setName("ball");
    ball->update();

    auto ground = objectManager.createObject();
    auto ground_config = ground->getParticleObjectConfig();
    ground_config->particle_radius = 0.075;
    ground_config->particle_mat = FIXED_BOUND;
    ground_config->phases.assign({0, 0});
    ground_config->shape = ObjectShape::Plane;
    ground_config->lb = {-20, -2, -20};
    ground_config->size = {40, 0, 40};
    ground_config->layer = 2;
    ground->setName("ground");
    ground->update();
//    ground->exportAsPly("F:\\DataSet.Research\\SCA2024\\ply\\diamonds\\ground",
//                        "1");

    /**  =============================================================
     * create solvers
     */
    SolverManager solverManager;
    auto solver = solverManager.createSolver<IMSCTSolver>();
    auto solver_config = dynamic_cast<IMSCTSolverConfig *>(solver->getConfig().get());
    solver_config->dt = 0.001;
    solver_config->gravity = {0, -9.8, 0};
    solver_config->scene_lb = {-30, -30, -30};
    solver_config->scene_size = {60, 60, 60};
    solver_config->export_data = true;
    solver_config->export_path = "F:\\DataSet.Research\\SCA2024\\ply\\diamonds\\f3_new";
    solver_config->export_gap = 10;
    solver_config->export_partial = "fluid";
    solver_config->export_phase = false;

    solver_config->max_neighborNum = 60;
    solver_config->rest_viscosity = 0.01;
    solver_config->rest_density = {980, 1000};
    solver_config->rest_rigid_density = 1000;
    solver_config->rest_bound_density = 1000;
    solver_config->div_free_threshold = 1e-4;
    solver_config->incompressible_threshold = 1e-4;
    solver_config->Cf = 0.01;
    solver_config->phase1_color = {50, 0, 200};
    solver_config->phase2_color = {200, 0, 50};

    solver_config->Cd0 = 1;
    solver_config->ct_thinning_exp0 = 0.01;
    solver_config->ct_relaxation_time = 0.005;
    solver_config->solution_vis_base = 10000;
    solver_config->polymer_vol_frac0 = 0; // threshold

//    solver->attachObject(ball);
    solver->attachObject(diamond);
    solver->attachObject(ground);

    /** world boundary **/
//    auto world = objectManager.createObject();
//    auto world_config = world->getParticleObjectConfig();
//    world_config->particle_radius = 0.075;
//    world_config->particle_mat = IMSCT_NONNEWTON;
//    world_config->phases.assign({0, 0});
//    world_config->shape = ObjectShape::Box;
//    world_config->lb = solver_config->scene_lb;
//    world_config->size = solver_config->scene_size;
//    world_config->layer = 1;
//    world->setName("world_boundary");
//    world->update();
//    world->exportAsPly("F:\\DataSet.Research\\SCA2024\\models\\world_boundary",
//                       "world");

    /**  =============================================================
     * run simulation
     */
    solver->run(2.1);

}
