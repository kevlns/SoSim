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
    auto buckling_plane = objectManager.createObject();
    auto buckling_plane_config = buckling_plane->getParticleObjectConfig();
    buckling_plane_config->particle_radius = 0.075;
    buckling_plane_config->particle_mat = IMSCT_NONNEWTON;
    buckling_plane_config->phases.assign({0.2, 0.8});
    buckling_plane_config->vel_start = {0, 0, 0};
    buckling_plane_config->model_file = "F:\\DataSet.Research\\SCA2024\\models\\scene_buckling\\ply\\fluid.ply";
    buckling_plane->setName("buckling_plane");
    buckling_plane->update();

    auto donut = objectManager.createObject();
    auto donut_config = donut->getParticleObjectConfig();
    donut_config->particle_radius = 0.075;
    donut_config->particle_mat = FIXED_BOUND;
    donut_config->phases.assign({0, 0});
    donut_config->model_file = "F:\\DataSet.Research\\SCA2024\\models\\scene_buckling\\ply\\donut.ply";
    donut->setName("donut");
    donut->update();

    auto plate = objectManager.createObject();
    auto plate_config = plate->getParticleObjectConfig();
    plate_config->particle_radius = 0.075;
    plate_config->particle_mat = FIXED_BOUND;
    plate_config->phases.assign({0, 0});
    plate_config->model_file = "F:\\DataSet.Research\\SCA2024\\models\\scene_buckling\\ply\\plate.ply";
    plate->setName("plate");
    plate->update();

    auto moving_tube = objectManager.createObject();
    auto moving_tube_config = moving_tube->getParticleObjectConfig();
    moving_tube_config->particle_radius = 0.075;
    moving_tube_config->particle_mat = MOVING_TUBE;
    moving_tube_config->phases.assign({0, 0});
    moving_tube_config->vel_start = {0, 0, 2};
    moving_tube_config->model_file = "F:\\DataSet.Research\\SCA2024\\models\\scene_buckling\\ply\\moving_tube.ply";
    moving_tube->setName("moving_tube");
    moving_tube->update();

    auto moving_cover = objectManager.createObject();
    auto moving_cover_config = moving_cover->getParticleObjectConfig();
    moving_cover_config->particle_radius = 0.075;
    moving_cover_config->particle_mat = MOVING_COVER;
    moving_cover_config->phases.assign({0, 0});
    moving_cover_config->vel_start = {0, -5, 2};
    moving_cover_config->model_file = "F:\\DataSet.Research\\SCA2024\\models\\scene_buckling\\ply\\moving_cover.ply";
    moving_cover->setName("moving_cover");
    moving_cover->update();


    /**  =============================================================
     * create solvers
     */
    SolverManager solverManager;
    auto solver = solverManager.createSolver<IMSCTSolver>();
    auto solver_config = dynamic_cast<IMSCTSolverConfig *>(solver->getConfig().get());
    solver_config->dt = 0.0001;
    solver_config->gravity = {0, -9.8, 0};
    solver_config->scene_lb = {-30, -30, -30};
    solver_config->scene_size = {60, 60, 60};
    solver_config->export_data = true;
    solver_config->export_path = "F:\\DataSet.Research\\SCA2024\\ply\\buckling\\fluid_v2";
    solver_config->export_gap = 250;
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
    solver_config->ct_thinning_exp0 = 0.2;
    solver_config->ct_relaxation_time = 0.0075;
    solver_config->solution_vis_base = 13000;
    solver_config->solution_vis_max = 13000;
    solver_config->polymer_vol_frac0 = 0; // threshold

    solver->attachObject(buckling_plane);
    solver->attachObject(donut);
    solver->attachObject(plate);
    solver->attachObject(moving_tube);
    solver->attachObject(moving_cover);

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
    solver->run(4);

}
