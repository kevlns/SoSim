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
    auto fluid_obj = objectManager.createObject();
    auto fluid_obj_config = fluid_obj->getParticleObjectConfig();
    fluid_obj_config->particle_radius = 0.075;
    fluid_obj_config->particle_mat = IMSCT_NONNEWTON;
    fluid_obj_config->shape = Cube;
    fluid_obj_config->lb = {-2,-2,-2,};
    fluid_obj_config->size = {4,4,4};
    fluid_obj_config->phases.assign({0.2, 0.8});
//    fluid_obj_config->model_file = "C:\\Users\\ADMIN\\Downloads\\simulation\\armadillo_ply\\fluid_armadillo.ply";
    fluid_obj->setName("fluid_obj");
    fluid_obj->update();

    auto bound_box = objectManager.createObject();
    auto bound_box_config = bound_box->getParticleObjectConfig();
    bound_box_config->particle_radius = 0.075;
    bound_box_config->particle_mat = FIXED_BOUND;
    bound_box_config->shape = Plane;
    bound_box_config->lb = {-5, -2.5,-5};
    bound_box_config->size = {10, 0,10};
    bound_box_config->layer = 2;
    bound_box_config->phases.assign({0, 0});
//    bound_box_config->model_file = "C:\\Users\\ADMIN\\Downloads\\simulation\\armadillo_ply\\box_armadillo.ply";
    bound_box->setName("bound_box");
    bound_box->update();

    /**  =============================================================
     * create solver
     */
    SolverManager solverManager;
    auto solver = solverManager.createSolver<IMMCTSolver>();
    auto solver_config = dynamic_cast<IMMCTSolverConfig *>(solver->getConfig().get());
    solver_config->dt = 0.0005;
    solver_config->gravity = {0, -9.8, 0};
    solver_config->scene_lb = {-30, -20, -30};
    solver_config->scene_size = {60, 60, 60};

    // export config
    solver_config->export_data = true;
    solver_config->export_path = "F:\\DataSet.Research\\SCA\\ply\\test2";
    solver_config->export_fps = 25;
    solver_config->export_partial = "fluid";
    solver_config->export_phase = false;

    // common
    solver_config->rest_density = {900, 1000};
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
    solver_config->ct_relaxation_time = 0.0005;
    solver_config->solution_vis_base = 15000;
    solver_config->solution_vis_max = 15000;
    solver_config->polymer_vol_frac0 = 0; // threshold
    solver_config->vis_bound_damp_factor = 0.5;

    /**  =============================================================
     * attach objects to solver
     */
    solver->attachObject(fluid_obj);
    solver->attachObject(bound_box);

    /**  =============================================================
     * run simulation
     */
    solver->run(5);

}