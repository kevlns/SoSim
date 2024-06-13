//
// Created by ADMIN on 2024/6/13.
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
    fluid_obj_config->particle_radius = 0.05;
    fluid_obj_config->particle_mat = COMMON_NEWTON;
    fluid_obj_config->shape = ObjectShape::Cube;
    fluid_obj_config->lb = {-1,-1,-1};
    fluid_obj_config->size = {2,2,2};
    fluid_obj->setName("fluid_obj");
    fluid_obj->update();

    auto box = objectManager.createObject();
    auto box_config = box->getParticleObjectConfig();
    box_config->shape = ObjectShape::Box;
    box_config->particle_radius = 0.05;
    box_config->particle_mat = FIXED_BOUND;
    box_config->lb = {-2, -1.5, -2};
    box_config->size = {4, 4, 4};
    box_config->layer = 1;
    box->setName("box");
    box->update();
    box->exportAsPly("F:\\DataSet.Research\\SoSimExamples\\PBF\\bound",
                       "box");

    /**  =============================================================
     * create solver
     */
    SolverManager solverManager;
    auto solver = solverManager.createSolver<PBFSolver>();
    auto solver_config = dynamic_cast<PBFSolverConfig *>(solver->getConfig().get());
    solver_config->dt = 0.001;
    solver_config->gravity = {0, -9.8, 0};
    solver_config->scene_lb = {-30, -20, -30};
    solver_config->scene_size = {60, 60, 60};

    // export config
    solver_config->export_data = true;
    solver_config->export_path = "F:\\DataSet.Research\\SoSimExamples\\PBF";
    solver_config->export_fps = 30;
    solver_config->export_partial = "fluid";

    // common
    solver_config->rest_density = 1000;
    solver_config->rest_rigid_density = 2000;
    solver_config->rest_bound_density = 2000;
    solver_config->pbf_iter_num = 5;

    /**  =============================================================
     * attach objects to solver
     */
    solver->attachObject(fluid_obj);
    solver->attachObject(box);

    /**  =============================================================
     * run simulation
     */
    solver->run(5);

}

