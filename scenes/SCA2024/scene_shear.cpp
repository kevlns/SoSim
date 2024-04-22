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
    auto solver = solverManager.createSolver<IMSCTSolver>();
    auto solver_config = dynamic_cast<IMSCTSolverConfig *>(solver->getConfig().get());
    solver_config->dt = 0.0001;
    solver_config->gravity = {0, -9.8, 0};
    solver_config->scene_lb = {-35, -35, -35};
    solver_config->scene_size = {70, 70, 70};
    solver_config->export_data = false;
    solver_config->export_path = "F:\\DataSet.Research\\SCA2024\\ply\\shear\\fluid_v1";
    solver_config->export_gap = 500;
    solver_config->export_partial = "fluid";
    solver_config->export_phase = true;

    solver_config->max_neighborNum = 60;
    solver_config->rest_viscosity = 0.008;
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
    solver_config->ct_relaxation_time = 0.006;
    solver_config->solution_vis_base = 4000;
    solver_config->solution_vis_max = 12000;
    solver_config->polymer_vol_frac0 = 0.6; // threshold

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
