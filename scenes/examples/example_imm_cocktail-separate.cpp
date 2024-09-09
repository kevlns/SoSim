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
    fluid_config->particle_radius = 0.05;
    fluid_config->particle_mat = COMMON_NEWTON;
    fluid_config->model_file = R"(D:\c-disk-data\Desktop\coffee\simulate\cocktail\300.ply)";
    fluid_config->phases = {0.1, 0.2, 0.3, 0.4};
    fluid_config->vel_start = {0, 0, 0};
    fluid->setName("fluid");
    fluid->update();

//    auto fluid2 = objectManager.createObject();
//    auto fluid2_config = fluid2->getParticleObjectConfig();
//    fluid2_config->particle_radius = 0.05;
//    fluid2_config->particle_mat = COMMON_NEWTON;
//    fluid2_config->shape = ObjectShape::Cube;
//    fluid2_config->lb = {-1, -1, -1};
//    fluid2_config->size = {2, 4, 2};
//    fluid2_config->phases = {0.25, 0.25, 0.25, 0.25};
//    fluid2_config->vel_start = {0, 0, 0};
//    fluid2->setName("fluid2");
//    fluid2->update();

    auto bound = objectManager.createObject();
    auto bound_config = bound->getParticleObjectConfig();
    bound_config->particle_radius = 0.05;
    bound_config->particle_mat = FIXED_BOUND;
    bound_config->shape = ObjectShape::Box;
    bound_config->lb = {-1.2, -1.2, -1.2};
    bound_config->size = {2.4, 8, 2.4};
    bound_config->layer = 1;
    bound_config->phases = {0, 0, 0 ,0};
    bound->setName("bound");
    bound->update();

    /**  =============================================================
     * create solver
     */
    SolverManager solverManager;
    auto solver = solverManager.createSolver<IMMSolver_v2>();
    auto solver_config = dynamic_cast<IMMSolverConfig_v2 *>(solver->getConfig().get());
    solver_config->dt = 0.002;
    solver_config->gravity = {0, -9.8, 0};
    solver_config->scene_lb = {-15, -15, -15};
    solver_config->scene_size = {30, 30, 30};

    // export config
    solver_config->export_data = true;
    solver_config->export_phase = true;
    solver_config->export_fps = 25;
    solver_config->export_partial = "fluid"; // fluid or all
    solver_config->export_path = "D:\\c-disk-data\\Desktop\\coffee\\simulate\\cocktail-separate";

    // common
    solver_config->phase_rest_density = {500, 650, 800, 1000};
    solver_config->phase_color = {
            {255, 0, 0},
            {0, 255, 0},
            {0, 0, 255},
            {255, 255, 0},
    };
    solver_config->phase_vis = {0.0, 0.0, 0.0, 0.0};
    solver_config->rest_rigid_density = 1000;
    solver_config->rest_bound_density = 1000;
    solver_config->rest_viscosity = 0.005;
    solver_config->div_free_threshold = 1e-4;
    solver_config->incompressible_threshold = 1e-4;
    solver_config->Cf = 0.05;
    solver_config->Cd = 0.6;

    /**  =============================================================
     * attach objects to solver
     */
    solver->attachObject(fluid);
//    solver->attachObject(fluid2);
    solver->attachObject(bound);

    /**  =============================================================
     * run simulation
     */
    solver->run(20);
}
