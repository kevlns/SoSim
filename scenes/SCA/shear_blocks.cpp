//
// Created by ADMIN on 2024/6/29.
//

#include "framework/object_manager.hpp"
#include "framework/solver_manager.hpp"

using namespace SoSim;

int main() {

    /**  =============================================================
     * create objects
     */
    ObjectManager objectManager;
    auto newton = objectManager.createObject();
    auto newton_config = newton->getParticleObjectConfig();
    newton_config->particle_radius = 0.075;
    newton_config->particle_mat = IMSCT_NONNEWTON;
    newton_config->phases.assign({1, 0});
    newton_config->vel_start = {0, 0, 0};
    newton_config->model_file = "F:\\DataSet.Research\\SCA2024\\models\\scene_shear\\ply\\newton.ply";
    newton->setName("newton");
    newton->update();

    auto thinning = objectManager.createObject();
    auto thinning_config = thinning->getParticleObjectConfig();
    thinning_config->particle_radius = 0.075;
    thinning_config->particle_mat = IMSCT_NONNEWTON;
    thinning_config->phases.assign({0.6, 0.4});
    thinning_config->model_file = "F:\\DataSet.Research\\SCA2024\\models\\scene_shear\\ply\\thinning.ply";
    thinning->setName("thinning");
    thinning->update();

    auto thickening = objectManager.createObject();
    auto thickening_config = thickening->getParticleObjectConfig();
    thickening_config->particle_radius = 0.075;
    thickening_config->particle_mat = IMSCT_NONNEWTON;
    thickening_config->phases.assign({0.2, 0.8});
    thickening_config->model_file = "F:\\DataSet.Research\\SCA2024\\models\\scene_shear\\ply\\thickening.ply";
    thickening->setName("thickening");
    thickening->update();

    auto plane = objectManager.createObject();
    auto plane_config = plane->getParticleObjectConfig();
    plane_config->particle_radius = 0.075;
    plane_config->particle_mat = FIXED_BOUND;
    plane_config->phases.assign({0, 0});
    plane_config->model_file = "F:\\DataSet.Research\\SCA2024\\models\\scene_shear\\ply\\plane.ply";
    plane->setName("plane");
    plane->update();

    auto world = objectManager.createObject();
    auto world_config = world->getParticleObjectConfig();
    world_config->particle_radius = 0.075;
    world_config->particle_mat = FIXED_BOUND;
    world_config->phases.assign({0, 0});
    world_config->model_file = "F:\\DataSet.Research\\SCA2024\\models\\world_boundary\\world.ply";
    world->setName("world");
    world->update();


    /**  =============================================================
     * create solvers
     */
    SolverManager solverManager;
    auto solver = solverManager.createSolver<IMMCTSolver>();
    auto solver_config = dynamic_cast<IMMCTSolverConfig *>(solver->getConfig().get());
    solver_config->dt = 0.0001;
    solver_config->gravity = {0, -9.8, 0};
    solver_config->scene_lb = {-15, -15, -15};
    solver_config->scene_size = {30, 30, 30};

    // export config
    solver_config->export_data = true;
    solver_config->export_path = "F:\\DataSet.Research\\SCA\\ply\\shear\\fluid_v1";
    solver_config->export_fps = 33;
    solver_config->export_partial = "fluid";
    solver_config->export_phase = true;

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
    solver_config->ct_thinning_exp0 = 0.15;
    solver_config->ct_relaxation_time = 0.007;
    solver_config->solution_vis_base = 10000;
    solver_config->solution_vis_max = 18000;
    solver_config->polymer_vol_frac0 = 0.6; // threshold
    solver_config->vis_bound_damp_factor = 0.05;

//    auto world = objectManager.createObject();
//    auto world_config = world->getParticleObjectConfig();
//    world_config->particle_radius = 0.075;
//    world_config->particle_mat = FIXED_BOUND;
//    world_config->phases.assign({0, 0});
//    world_config->shape = Box;
//    world_config->lb = solver_config->scene_lb;
//    world_config->size = solver_config->scene_size;
//    world_config->layer = 1;
//    world->setName("world");
//    world->update();
//    world->exportAsPly("F:\\DataSet.Research\\SCA2024\\models\\world_boundary", "world");

    solver->attachObject(newton);
    solver->attachObject(thinning);
    solver->attachObject(thickening);
    solver->attachObject(plane);
    solver->attachObject(world);

    /**  =============================================================
     * run simulation
     */
    solver->run(5);

}