//
// Created by ADMIN on 2024/6/27.
//

#include "framework/object_manager.hpp"
#include "framework/solver_manager.hpp"
#include "core/math/matrix.hpp"

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
    fluid_config->shape = ObjectShape::Cube;
    fluid_config->lb = {-1,-1,-1};
    fluid_config->size = {2,2,2};
    fluid_config->phases = {1};
    fluid->setName("fluid");
    fluid->update();

    auto bound = objectManager.createObject();
    auto bound_config = bound->getParticleObjectConfig();
    bound_config->particle_radius = 0.05;
    bound_config->particle_mat = FIXED_BOUND;
    bound_config->shape = ObjectShape::Plane;
    bound_config->lb = {-5,-3,-5};
    bound_config->size = {10,0,10};
    bound_config->layer = 2;
    bound_config->phases = {0};
    bound->setName("bound");
    bound->update();

    /**  =============================================================
     * create solver
     */
    SolverManager solverManager;
    auto solver = solverManager.createSolver<IMMSolver_v2>();
    auto solver_config = dynamic_cast<IMMSolverConfig_v2 *>(solver->getConfig().get());
    solver_config->dt = 0.0015;
    solver_config->gravity = {0, -9.8, 0};
    solver_config->scene_lb = {-15, -15, -15};
    solver_config->scene_size = {30, 30, 30};

    // export config
    solver_config->export_data = true;
    solver_config->export_path = "E:\\Research_DATA_TEMP\\imm_v2_test";
    solver_config->export_fps = 500;
    solver_config->export_partial = "fluid";
    solver_config->export_phase = false;

    // common
    solver_config->phase_rest_densities = {1000.f};
    solver_config->phase_colors = {{0.f, 0.f, 1.f}};
    solver_config->phase_vis = {0.01f};
    solver_config->div_free_threshold = 1e-4;
    solver_config->incompressible_threshold = 1e-4;
    solver_config->Cf = 0.15;
    solver_config->Cd = 0.7;

    /**  =============================================================
     * attach objects to solver
     */
    solver->attachObject(fluid);
    solver->attachObject(bound);

    /**  =============================================================
     * run simulation
     */
    solver->run(0.1);

    return 0;
}