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
    fluid_config->particle_radius = 0.025;
    fluid_config->particle_mat = IMSCT_NONNEWTON;
    fluid_config->phases.assign({1, 0});
    fluid_config->model_file = "C:\\Users\\ADMIN\\Downloads\\simulation\\armadillo_ply\\fluid_armadillo.ply";
    fluid->setName("fluid");
    fluid->update();

    auto bound = objectManager.createObject();
    auto bound_config = bound->getParticleObjectConfig();
    bound_config->particle_radius = 0.025;
    bound_config->particle_mat = FIXED_BOUND;
    bound_config->phases.assign({0, 0});
    bound_config->model_file = "C:\\Users\\ADMIN\\Downloads\\simulation\\armadillo_ply\\box_armadillo.ply";
    bound->setName("bound");
    bound->update();

    /**  =============================================================
     * create solvers
     */
    SolverManager solverManager;
    auto solver = solverManager.createSolver<IMSCTSolver>();
    auto solver_config = dynamic_cast<IMSCTSolverConfig *>(solver->getConfig().get());
    solver_config->dt = 0.001;
    solver_config->gravity = {0, -9.8, 0};
    solver_config->scene_lb = {-30, -20, -30};
    solver_config->scene_size = {60, 60, 60};
    solver_config->export_data = true;
    solver_config->export_path = "F:\\DataSet.Research\\XYH\\armadillo";
    solver_config->export_fps = 25;
    solver_config->export_partial = "fluid";
    solver_config->export_phase = false;

    solver_config->max_neighborNum = 60;
    solver_config->rest_viscosity = 0.001;
    solver_config->phase1_vis = 0.001;
    solver_config->phase2_vis = 0.001;
    solver_config->rest_density = {980, 1000};
    solver_config->rest_rigid_density = 1000;
    solver_config->rest_bound_density = 1000;
    solver_config->div_free_threshold = 1e-4;
    solver_config->incompressible_threshold = 1e-4;
    solver_config->Cf = 0.25;
    solver_config->phase1_color = {50, 0, 200};
    solver_config->phase2_color = {200, 0, 50};

    solver_config->Cd0 = 0.4;

    solver->attachObject(fluid);
    solver->attachObject(bound);

    /**  =============================================================
     * run simulation
     */
    solver->run(40);

}
