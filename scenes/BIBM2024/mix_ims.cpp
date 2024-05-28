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
    auto water_col = objectManager.createObject();
    auto water_col_config = water_col->getParticleObjectConfig();
    water_col_config->particle_radius = 0.075;
    water_col_config->particle_mat = IMSCT_NONNEWTON;
    water_col_config->phases.assign({0.9, 0.1});
    water_col_config->model_file = "F:\\DataSet.Research\\BIBM2024\\model\\compare\\rescale\\water_col.ply";
    water_col->setName("water_col");
    water_col->update();

    auto polymer_col = objectManager.createObject();
    auto polymer_col_config = polymer_col->getParticleObjectConfig();
    polymer_col_config->particle_radius = 0.075;
    polymer_col_config->particle_mat = IMSCT_NONNEWTON;
    polymer_col_config->phases.assign({0.1, 0.9});
    polymer_col_config->model_file = "F:\\DataSet.Research\\BIBM2024\\model\\compare\\rescale\\polymer_col.ply";
    polymer_col->setName("polymer_col");
    polymer_col->update();

    auto stir_fan = objectManager.createObject();
    auto stir_fan_config = stir_fan->getParticleObjectConfig();
    stir_fan_config->particle_radius = 0.075;
    stir_fan_config->particle_mat = STIR_FAN;
    stir_fan_config->phases.assign({0, 0});
    stir_fan_config->model_file = "F:\\DataSet.Research\\BIBM2024\\model\\compare\\rescale\\fan.ply";
    stir_fan->setName("stir_fan");
    stir_fan->update();

    auto bowl = objectManager.createObject();
    auto bowl_config = bowl->getParticleObjectConfig();
    bowl_config->particle_radius = 0.075;
    bowl_config->particle_mat = FIXED_BOUND;
    bowl_config->phases.assign({0, 0});
    bowl_config->model_file = "F:\\DataSet.Research\\BIBM2024\\model\\compare\\rescale\\bowl.ply";
    bowl->setName("bowl");
    bowl->update();

    auto bowl_up = objectManager.createObject();
    auto bowl_up_config = bowl_up->getParticleObjectConfig();
    bowl_up_config->particle_radius = 0.075;
    bowl_up_config->particle_mat = FIXED_BOUND;
    bowl_up_config->phases.assign({0, 0});
    bowl_up_config->model_file = "F:\\DataSet.Research\\BIBM2024\\model\\compare\\rescale\\bowl_up.ply";
    bowl_up->setName("bowl_up");
    bowl_up->update();

    /**  =============================================================
     * create solvers
     */

    /*
     * scene 1_10: dt=0.001; export_gap=300; vis_base=vis_max=100
     * scene 1_100: dt=0.001; export_gap=300; vis_base=1000; vis_max=2000
     * scene 1_1000: ...
     */

    SolverManager solverManager;
    auto solver = solverManager.createSolver<IMSCTSolver>();
    auto solver_config = dynamic_cast<IMSCTSolverConfig *>(solver->getConfig().get());
    solver_config->dt = 0.001;
    solver_config->gravity = {0, -9.8, 0};
    solver_config->scene_lb = {-30, -20, -30};
    solver_config->scene_size = {60, 60, 60};
    solver_config->export_data = true;
    solver_config->export_path = "F:\\DataSet.Research\\BIBM2024\\ply\\compare\\ims\\1_10";
    solver_config->export_gap = 300;
    solver_config->export_partial = "fluid";
    solver_config->export_phase = false;

    solver_config->max_neighborNum = 60;
    solver_config->rest_viscosity = 0.01;
    solver_config->phase1_vis = 0.01;
    solver_config->phase2_vis = 0.1;
    solver_config->rest_density = {980, 1000};
    solver_config->rest_rigid_density = 1000;
    solver_config->rest_bound_density = 1000;
    solver_config->div_free_threshold = 1e-4;
    solver_config->incompressible_threshold = 1e-4;
    solver_config->Cf = 0.25;
    solver_config->phase1_color = {50, 0, 200};
    solver_config->phase2_color = {200, 0, 50};

    solver_config->Cd0 = 0.4;

    solver->attachObject(water_col);
    solver->attachObject(polymer_col);
    solver->attachObject(stir_fan);
    solver->attachObject(bowl);
    solver->attachObject(bowl_up);

    /**  =============================================================
     * run simulation
     */
    solver->run(25);

}
