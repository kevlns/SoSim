//
// Created by ADMIN on 2024/6/27.
//

#include "framework/object_manager.hpp"
#include "framework/solver_manager.hpp"

using namespace SoSim;

int main() {

    /**  =============================================================
     * create objects
     */
    ObjectManager objectManager;
    auto fluid_1 = objectManager.createObject();
    auto fluid_1_config = fluid_1->getParticleObjectConfig();
    fluid_1_config->particle_radius = 0.075;
    fluid_1_config->particle_mat = IMSCT_NONNEWTON;
    fluid_1_config->phases.assign({0.99, 0.01});
    fluid_1_config->model_file = "F:\\DataSet.Research\\SCA\\model\\stir_compare\\rescale\\ply\\fluid_1.ply";
    fluid_1->setName("fluid_1");
    fluid_1->update();

    auto fluid_2 = objectManager.createObject();
    auto fluid_2_config = fluid_2->getParticleObjectConfig();
    fluid_2_config->particle_radius = 0.075;
    fluid_2_config->particle_mat = IMSCT_NONNEWTON;
    fluid_2_config->phases.assign({0.1, 0.9});
    fluid_2_config->model_file = "F:\\DataSet.Research\\SCA\\model\\stir\\rescale\\ply\\fluid_2.ply";
    fluid_2->setName("fluid_2");
    fluid_2->update();

    auto stir_fan = objectManager.createObject();
    auto stir_fan_config = stir_fan->getParticleObjectConfig();
    stir_fan_config->particle_radius = 0.075;
    stir_fan_config->particle_mat = STIR_FAN;
    stir_fan_config->phases.assign({0, 0});
    stir_fan_config->model_file = "F:\\DataSet.Research\\SCA\\model\\stir\\rescale\\ply\\stir_fan.ply";
    stir_fan->setName("stir_fan");
    stir_fan->update();

    auto stir_bowl = objectManager.createObject();
    auto stir_bowl_config = stir_bowl->getParticleObjectConfig();
    stir_bowl_config->particle_radius = 0.075;
    stir_bowl_config->particle_mat = DYNAMIC_RIGID;
    stir_bowl_config->phases.assign({0, 0});
    stir_bowl_config->model_file = "F:\\DataSet.Research\\SCA\\model\\stir\\rescale\\ply\\stir_bowl.ply";
    stir_bowl->setName("stir_bowl");
    stir_bowl->update();

    auto scaffold = objectManager.createObject();
    auto scaffold_config = scaffold->getParticleObjectConfig();
    scaffold_config->particle_radius = 0.075;
    scaffold_config->particle_mat = FIXED_BOUND;
    scaffold_config->phases.assign({0, 0});
    scaffold_config->model_file = "F:\\DataSet.Research\\SCA\\model\\stir\\rescale\\ply\\scaffold.ply";
    scaffold->setName("scaffold");
    scaffold->update();

    auto wooden_bowl = objectManager.createObject();
    auto wooden_bowl_config = wooden_bowl->getParticleObjectConfig();
    wooden_bowl_config->particle_radius = 0.075;
    wooden_bowl_config->particle_mat = FIXED_BOUND;
    wooden_bowl_config->phases.assign({0, 0});
    wooden_bowl_config->model_file = "F:\\DataSet.Research\\SCA\\model\\stir\\rescale\\ply\\wooden_bowl.ply";
    wooden_bowl->setName("wooden_bowl");
    wooden_bowl->update();

    auto bound_box = objectManager.createObject();
    auto bound_box_config = bound_box->getParticleObjectConfig();
    bound_box_config->particle_radius = 0.075;
    bound_box_config->particle_mat = FIXED_BOUND;
    bound_box_config->phases.assign({0, 0});
    bound_box_config->model_file = "F:\\DataSet.Research\\SCA\\model\\stir\\rescale\\ply\\bound_box.ply";
    bound_box->setName("bound_box");
    bound_box->update();

    /**  =============================================================
     * create solver
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
    solver_config->export_path = "F:\\DataSet.Research\\SCA\\ply\\stir_v2";
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
    solver_config->Cf = 0.2;
    solver_config->Cd0 = 0.6;

    // if enable_CT == true
    solver_config->ct_thinning_exp0 = 0.1;
    solver_config->ct_relaxation_time = 0.018;
    solver_config->solution_vis_base = 10000;
    solver_config->solution_vis_max = 15000;
    solver_config->polymer_vol_frac0 = 0.65; // threshold
    solver_config->vis_bound_damp_factor = 0.5;

    /**  =============================================================
     * attach objects to solver
     */
    solver->attachObject(fluid_1);
    solver->attachObject(fluid_2);
    solver->attachObject(stir_fan);
    solver->attachObject(stir_bowl);
    solver->attachObject(scaffold);
    solver->attachObject(wooden_bowl);
    solver->attachObject(bound_box);

    /**  =============================================================
     * run simulation
     */
    solver->run(20);

}