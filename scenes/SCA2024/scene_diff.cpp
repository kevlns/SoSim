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

    auto cube = objectManager.createObject();
    auto cube_config = cube->getParticleObjectConfig();
    cube_config->particle_radius = 0.075;
    cube_config->particle_mat = IMSCT_NONNEWTON;
    cube_config->phases.assign({0.3, 0.7});
    cube_config->shape = ObjectShape::Cube;
    cube_config->lb = {-2, -4, -2};
    cube_config->size = {4, 3, 4};
    cube->setName("cube");
    cube->update();

    auto ball = objectManager.createObject();
    auto ball_config = ball->getParticleObjectConfig();
    ball_config->particle_radius = 0.075;
    ball_config->particle_mat = IMSCT_NONNEWTON;
    ball_config->phases.assign({0, 1});
    ball_config->shape = ObjectShape::Sphere;
    ball_config->center = {0, -1.25, 0};
    ball_config->volume_radius = 1.5;
    ball_config->vel_start = {0, -1.7, 0};
    ball->setName("ball");
    ball->update();

    auto box = objectManager.createObject();
    auto box_config = box->getParticleObjectConfig();
    box_config->particle_radius = 0.075;
    box_config->particle_mat = FIXED_BOUND;
    box_config->phases.assign({0, 0});
    box_config->shape = ObjectShape::Box;
    box_config->lb = {-2.4, -4.3, -2.4};
    box_config->size = {4.8, 10, 4.8};
    box_config->layer = 2;
    box->setName("box");
    box->update();
    box->exportAsPly("F:\\DataSet.Research\\SCA2024\\ply\\diff\\box",
                     "box");

    auto world = objectManager.createObject();
    auto world_config = world->getParticleObjectConfig();
    world_config->particle_radius = 0.075;
    world_config->particle_mat = FIXED_BOUND;
    world_config->phases.assign({0, 0});
    world_config->model_file = "F:\\DataSet.Research\\SCA2024\\models\\world_boundary\\world.ply";
    world->setName("world");
    world->update();

//    auto cube = objectManager.createObject();
//    auto cube_config = cube->getParticleObjectConfig();
//    cube_config->particle_radius = 0.075;
//    cube_config->particle_mat = IMSCT_NONNEWTON;
//    cube_config->phases.assign({1, 0});
//    cube_config->vel_start = {0, 0, 0};
//    cube_config->model_file = "F:\\DataSet.Research\\SCA2024\\models\\scene_diff\\cube_fluid.ply";
//    cube->setName("cube");
//    cube->update();

//    auto ball = objectManager.createObject();
//    auto ball_config = ball->getParticleObjectConfig();
//    ball_config->particle_radius = 0.075;
//    ball_config->particle_mat = IMSCT_NONNEWTON;
//    ball_config->phases.assign({0, 1});
//    ball_config->model_file = "F:\\DataSet.Research\\SCA2024\\models\\scene_diff\\sphere_fluid.ply";
//    ball_config->vel_start = {-3, 0, 0};
//    ball->setName("ball");
//    ball->update();

//    auto container = objectManager.createObject();
//    auto ground_config = container->getParticleObjectConfig();
//    ground_config->particle_radius = 0.075;
//    ground_config->particle_mat = FIXED_BOUND;
//    ground_config->phases.assign({0, 0});
//    ground_config->model_file = "F:\\DataSet.Research\\SCA2024\\models\\scene_diff\\cube_container.ply";
//    container->setName("container");
//    container->update();

    /**  =============================================================
     * create solvers
     */
    SolverManager solverManager;
    auto solver = solverManager.createSolver<IMSCTSolver>();
    auto solver_config = dynamic_cast<IMSCTSolverConfig *>(solver->getConfig().get());
    solver_config->dt = 0.0015;
    solver_config->gravity = {0, -9.8, 0};
    solver_config->scene_lb = {-35, -35, -35};
    solver_config->scene_size = {70, 70, 70};
    solver_config->export_data = true;
    solver_config->export_path = "F:\\DataSet.Research\\SCA2024\\ply\\diff\\fluid_ben_separate_v2";
    solver_config->export_gap = 150;
    solver_config->export_partial = "fluid";
    solver_config->export_phase = false;

    solver_config->max_neighborNum = 60;
    solver_config->rest_viscosity = 0.008;
    solver_config->rest_density = {800, 1000};
    solver_config->rest_rigid_density = 1000;
    solver_config->rest_bound_density = 1000;
    solver_config->div_free_threshold = 1e-4;
    solver_config->incompressible_threshold = 1e-4;
    solver_config->Cf = 0.0005;
    solver_config->phase1_color = {50, 0, 200};
    solver_config->phase2_color = {200, 0, 50};

    solver_config->Cd0 = 0.3;
    solver_config->ct_thinning_exp0 = 0.2;
    solver_config->ct_relaxation_time = 0.004;
    solver_config->solution_vis_base = 100;
    solver_config->solution_vis_max = 100;
    solver_config->polymer_vol_frac0 = 0; // threshold

    solver->attachObject(cube);
//    solver->attachObject(ball);
    solver->attachObject(box);
    solver->attachObject(world);
//    solver->attachObject(container);

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
    solver->run(20);

}
