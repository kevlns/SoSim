//
// Created by ADMIN on 2024/4/9.
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
    ball_config->center = {0, 3, 0};
    ball_config->volume_radius = 3;
    ball->setName("ball");
    ball->update();

    auto ground = objectManager.createObject();
    auto ground_config = ground->getParticleObjectConfig();
    ground_config->particle_radius = 0.075;
    ground_config->particle_mat = FIXED_BOUND;
    ground_config->phases.assign({0, 0});
    ground_config->model_file = "F:\\DataSet.Research\\SCA\\model\\diamonds\\rescale\\ground_plane.ply";
    ground->setName("ground");
    ground->update();

    auto bound = objectManager.createObject();
    auto bound_config = bound->getParticleObjectConfig();
    bound_config->particle_radius = 0.075;
    bound_config->particle_mat = FIXED_BOUND;
    bound_config->phases.assign({0, 0});
    bound_config->model_file = "F:\\DataSet.Research\\SCA\\ply\\diamonds\\bound\\bound.ply";
    bound->setName("bound");
    bound->update();

    /**  =============================================================
     * create solver
     */
    SolverManager solverManager;
    auto solver = solverManager.createSolver<IMMCTSolver>();
    auto solver_config = dynamic_cast<IMMCTSolverConfig *>(solver->getConfig().get());
    solver_config->dt = 0.0002;
    solver_config->gravity = {0, -9.8, 0};
    solver_config->scene_lb = {-30, -15, -30};
    solver_config->scene_size = {60, 30, 60};

    // export config
    solver_config->export_data = true;
    solver_config->export_path = "F:\\DataSet.Research\\SCA\\ply\\diamonds\\e1";
    solver_config->export_fps = 33;
    solver_config->export_partial = "fluid";
    solver_config->export_phase = false;

    // common
    solver_config->rest_density = {980, 1000};
    solver_config->rest_rigid_density = 1000;
    solver_config->rest_bound_density = 1000;
    solver_config->div_free_threshold = 1e-4;
    solver_config->incompressible_threshold = 1e-4;
    solver_config->phase1_color = {5, 0, 250};
    solver_config->phase2_color = {250, 0, 5};
    solver_config->Cf = 0;
    solver_config->Cd0 = 1;

    // if enable_CT == true
    solver_config->ct_thinning_exp0 = 0;
    solver_config->ct_relaxation_time = 0.01;
    solver_config->solution_vis_base = 10000;
    solver_config->solution_vis_max = 10000;
    solver_config->polymer_vol_frac0 = 0; // threshold
    solver_config->vis_bound_damp_factor = 0.01;

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
//    world->exportAsPly("F:\\DataSet.Research\\SCA\\ply\\diamonds\\bound",
//                       "bound");

    /**  =============================================================
     * attach objects to solver
     */
    solver->attachObject(diamond);
//    solver->attachObject(ball);
    solver->attachObject(ground);
    solver->attachObject(bound);

    /**  =============================================================
     * run simulation
     */
    solver->run(2.5);

}