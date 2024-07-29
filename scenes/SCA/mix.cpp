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
    auto water_cylinder = objectManager.createObject();
    auto water_cylinder_config = water_cylinder->getParticleObjectConfig();
    water_cylinder_config->particle_radius = 0.075;
    water_cylinder_config->particle_mat = IMSCT_NONNEWTON;
    water_cylinder_config->phases.assign({0.95, 0.05});
    water_cylinder_config->model_file = "F:\\DataSet.Research\\SCA2024\\models\\scene_stir\\ply\\water_cylinder.ply";
    water_cylinder->setName("water_cylinder");
    water_cylinder->update();

    auto polymer_cylinder = objectManager.createObject();
    auto polymer_cylinder_config = polymer_cylinder->getParticleObjectConfig();
    polymer_cylinder_config->particle_radius = 0.075;
    polymer_cylinder_config->particle_mat = IMSCT_NONNEWTON;
    polymer_cylinder_config->phases.assign({0.05, 0.95});
    polymer_cylinder_config->model_file = "F:\\DataSet.Research\\SCA2024\\models\\scene_stir\\ply\\polymer_cylinder.ply";
    polymer_cylinder->setName("polymer_cylinder");
    polymer_cylinder->update();

    auto tube = objectManager.createObject();
    auto tube_config = tube->getParticleObjectConfig();
    tube_config->particle_radius = 0.075;
    tube_config->particle_mat = FIXED_BOUND;
    tube_config->phases.assign({0, 0});
    tube_config->model_file = "F:\\DataSet.Research\\SCA2024\\models\\scene_stir\\ply\\tube.ply";
    tube->setName("tube");
    tube->update();

    auto rocks = objectManager.createObject();
    auto rocks_config = rocks->getParticleObjectConfig();
    rocks_config->particle_radius = 0.075;
    rocks_config->particle_mat = FIXED_BOUND;
    rocks_config->phases.assign({0, 0});
    rocks_config->model_file = "F:\\DataSet.Research\\SCA2024\\models\\scene_stir\\ply\\rocks.ply";
    rocks->setName("rocks");
    rocks->update();

//    auto world = objectManager.createObject();
//    auto world_config = world->getParticleObjectConfig();
//    world_config->particle_radius = 0.075;
//    world_config->particle_mat = FIXED_BOUND;
//    world_config->phases.assign({0, 0});
//    world_config->model_file = "F:\\DataSet.Research\\SCA2024\\models\\world_boundary\\world.ply";
//    world->setName("world");
//    world->update();

    auto ground = objectManager.createObject();
    auto ground_config = ground->getParticleObjectConfig();
    ground_config->particle_radius = 0.075;
    ground_config->particle_mat = FIXED_BOUND;
    ground_config->phases.assign({0, 0});
    ground_config->shape = ObjectShape::Plane;
    ground_config->lb = {-30, -20.5, -30};
    ground_config->size = {60, 0, 60};
    ground_config->layer = 2;
    ground->setName("ground");
    ground->update();
//    ground->exportAsPly("F:\\DataSet.Research\\SCA2024\\ply\\stir\\ground",
//                        "ground");

    auto mass_cylinder = objectManager.createObject();
    auto mass_cylinder_config = mass_cylinder->getParticleObjectConfig();
    mass_cylinder_config->particle_radius = 0.075;
    mass_cylinder_config->particle_mat = FIXED_BOUND;
    mass_cylinder_config->phases.assign({0, 0});
    mass_cylinder_config->model_file = "F:\\DataSet.Research\\SCA2024\\models\\scene_stir\\ply\\mass_cylinder.ply";
    mass_cylinder->setName("mass_cylinder");
    mass_cylinder->update();

    auto bowl = objectManager.createObject();
    auto bowl_config = bowl->getParticleObjectConfig();
    bowl_config->particle_radius = 0.075;
    bowl_config->particle_mat = MOVING_BOWL;
    bowl_config->phases.assign({0, 0});
    bowl_config->model_file = "F:\\DataSet.Research\\SCA2024\\models\\scene_stir\\ply\\bowl.ply";
    bowl->setName("bowl");
    bowl->update();

    auto stir_fan = objectManager.createObject();
    auto stir_fan_config = stir_fan->getParticleObjectConfig();
    stir_fan_config->particle_radius = 0.075;
    stir_fan_config->particle_mat = STIR_FAN;
    stir_fan_config->phases.assign({0, 0});
    stir_fan_config->model_file = "F:\\DataSet.Research\\SCA2024\\models\\scene_stir\\ply\\fan_simple.ply";
    stir_fan->setName("stir_fan");
    stir_fan->update();

    /**  =============================================================
     * create solver
     */
    SolverManager solverManager;
    auto solver = solverManager.createSolver<IMMCTSolver>();
    auto solver_config = dynamic_cast<IMMCTSolverConfig *>(solver->getConfig().get());
    solver_config->dt = 0.001;
    solver_config->gravity = {0, -9.8, 0};
    solver_config->scene_lb = {-30, -15, -30};
    solver_config->scene_size = {60, 30, 60};

    // export config
    solver_config->export_data = true;
    solver_config->export_path = "F:\\DataSet.Research\\SCA\\ply\\diamonds\\f4";
    solver_config->export_fps = 33;
    solver_config->export_partial = "fluid";
    solver_config->export_phase = false;

    // common
    solver_config->rest_density = {1000, 1000};
    solver_config->rest_rigid_density = 1000;
    solver_config->rest_bound_density = 1000;
    solver_config->div_free_threshold = 1e-4;
    solver_config->incompressible_threshold = 1e-4;
    solver_config->phase1_color = {5, 0, 250};
    solver_config->phase2_color = {250, 0, 5};
    solver_config->Cf = 0.2;
    solver_config->Cd0 = 0.6;

    // if enable_CT == true
    solver_config->ct_thinning_exp0 = 0.35;
    solver_config->ct_relaxation_time = 0.005;
    solver_config->solution_vis_base = 5000;
    solver_config->solution_vis_max = 5000;
    solver_config->polymer_vol_frac0 = 0.7; // threshold
    solver_config->vis_bound_damp_factor = 0.01;

    /**  =============================================================
     * attach objects to solver
     */
//    solver->attachObject(water_cylinder);
//    solver->attachObject(polymer_cylinder);
//    solver->attachObject(tube);
//    solver->attachObject(rocks);
//    solver->attachObject(ground);
//    solver->attachObject(bowl);
//    solver->attachObject(world);
//    solver->attachObject(stir_fan);
//    solver->attachObject(mass_cylinder);

        auto world = objectManager.createObject();
    auto world_config = world->getParticleObjectConfig();
    world_config->particle_radius = 0.075;
    world_config->particle_mat = IMSCT_NONNEWTON;
    world_config->phases.assign({0, 0});
    world_config->shape = ObjectShape::Box;
    world_config->lb = solver_config->scene_lb;
    world_config->size = solver_config->scene_size;
    world_config->layer = 1;
    world->setName("world_boundary");
    world->update();
    world->exportAsPly("F:\\DataSet.Research\\SCA2024\\models\\world_boundary",
                       "world");

    /**  =============================================================
     * run simulation
     */
//    solver->run(5);

}