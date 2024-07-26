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
    auto fluid = objectManager.createObject();
    auto fluid_config = fluid->getParticleObjectConfig();
    fluid_config->particle_radius = 0.075;
    fluid_config->particle_mat = IMSCT_NONNEWTON;
    fluid_config->phases.assign({0.1, 0.9});
    fluid_config->vel_start = {0, 0, 0};
    fluid_config->model_file = "F:\\DataSet.Research\\SCA\\model\\buckling\\rescale\\fluid.ply";
    fluid->setName("fluid");
    fluid->update();

    auto ground = objectManager.createObject();
    auto ground_config = ground->getParticleObjectConfig();
    ground_config->particle_radius = 0.075;
    ground_config->particle_mat = FIXED_BOUND;
    ground_config->phases.assign({0, 0});
    ground_config->model_file = "F:\\DataSet.Research\\SCA\\model\\buckling\\rescale\\ground.ply";
    ground->setName("ground");
    ground->update();

    /**  =============================================================
     * create solver
     */
    SolverManager solverManager;
    auto solver = solverManager.createSolver<IMMCTSolver>();
    auto solver_config = dynamic_cast<IMMCTSolverConfig *>(solver->getConfig().get());
    solver_config->dt = 0.0005;
    solver_config->gravity = {0, -9.8, 0};
    solver_config->scene_lb = {-30, -15, -30};
    solver_config->scene_size = {60, 30, 60};

    // export config
    solver_config->export_data = true;
    solver_config->export_path = "F:\\DataSet.Research\\SCA\\ply\\buckling";
    solver_config->export_fps = 33;
    solver_config->export_partial = "fluid";
    solver_config->export_phase = false;

    // common
    solver_config->rest_density = {1000, 1000};
    solver_config->rest_rigid_density = 1000;
    solver_config->rest_bound_density = 1000;
    solver_config->div_free_threshold = 1e-4;
    solver_config->incompressible_threshold = 1e-4;
    solver_config->phase1_color = {254, 254, 254};
    solver_config->phase2_color = {255, 2, 2};
    solver_config->Cf = 0.0;
    solver_config->Cd0 = 1.0;

    // if enable_CT == true
    solver_config->ct_thinning_exp0 = 0;
    solver_config->ct_relaxation_time = 0.009; // 0.012
    solver_config->solution_vis_base = 18000;
    solver_config->solution_vis_max = 18000;
    solver_config->polymer_vol_frac0 = 0; // threshold
    solver_config->vis_bound_damp_factor = 0.9;

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
    solver->attachObject(fluid);
    solver->attachObject(ground);

    /**  =============================================================
     * run simulation
     */
    solver->run(5);

}