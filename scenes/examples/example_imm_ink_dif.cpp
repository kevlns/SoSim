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
    fluid_config->particle_mat = COMMON_NEWTON;
    fluid_config->model_file = R"(D:\xuyuhang\simulation\multi-fluid\ink_dif\cylinder\ply_phase\1590.ply)";
//    fluid_config->model_file = R"(D:\xuyuhang\sample_ply\multiple-cylinder-water-v2.ply)";
    fluid_config->phases = {1.0, 0.0};
    fluid_config->vel_start = {0, 0, 0};
    fluid->setName("fluid");
    fluid->update();

    auto fluid2 = objectManager.createObject();
    auto fluid2_config = fluid2->getParticleObjectConfig();
    fluid2_config->particle_radius = 0.025;
    fluid2_config->particle_mat = COMMON_NEWTON;
    fluid2_config->model_file = R"(D:\xuyuhang\sample_ply\multiple-waterdrops.ply)";
    fluid2_config->phases = {0.3, 0.7};
    fluid2_config->vel_start = {0.2, 1.0, 0.0};
    fluid2->setName("fluid2");
    fluid2->update();

    auto bound = objectManager.createObject();
    auto bound_config = bound->getParticleObjectConfig();
    bound_config->particle_radius = 0.025;
    bound_config->particle_mat = FIXED_BOUND;
    bound_config->model_file = R"(D:\xuyuhang\sample_ply\multiple-cylinder-box.ply)";
//    bound_config->model_file = R"(D:\xuyuhang\sample_ply\cylinder-box\multiple-cylinder-box-ratio_10.ply)";
    bound_config->phases = {0, 0};
    bound->setName("bound");
    bound->update();

    /**  =============================================================
     * create solver
     */
    SolverManager solverManager;
    auto solver = solverManager.createSolver<IMMSolver_v2>();
    auto solver_config = dynamic_cast<IMMSolverConfig_v2 *>(solver->getConfig().get());
    solver_config->dt = 0.001;
    solver_config->gravity = {0, -9.8, 0};
    solver_config->scene_lb = {-15, -15, -15};
    solver_config->scene_size = {30, 30, 30};

    // export config
    solver_config->export_data = true;
    solver_config->export_phase = true;
    solver_config->export_fps = 50;
    solver_config->export_partial = "fluid"; // fluid or all
//    solver_config->export_path = "D:\\xuyuhang\\simulation\\multi-fluid\\render\\ink-dif\\cylinder\\cf0.01-cd0.8-ink-vel_0.6_3.0_0.0-phase_0.7_0.3\\ply";
//    solver_config->export_path = "D:\\xuyuhang\\simulation\\multi-fluid\\render\\ink-dif\\cylinder\\standing-water";

    // common
    solver_config->phase_rest_density = {1000, 2000};
    solver_config->phase_color = {
            {255, 0, 0},
            {0, 255, 0},
    };
    solver_config->phase_vis = {0.0, 0.0};
    solver_config->rest_rigid_density = 1000;
    solver_config->rest_bound_density = 2000;
    solver_config->rest_viscosity = 0.005;
    solver_config->div_free_threshold = 1e-4;
    solver_config->incompressible_threshold = 1e-4;
    solver_config->Cf = 0.05;
    solver_config->Cd = 0.6;

    // 根据cf cd vel 等设置export_path
    solver_config->export_path = "D:\\xuyuhang\\simulation\\multi-fluid\\render\\ink-dif\\cylinder\\cf"
            + std::to_string(solver_config->Cf) + "-cd" + std::to_string(solver_config->Cd)
            + "-ink-vel_" + std::to_string(fluid2_config->vel_start.x) + "_" + std::to_string(fluid2_config->vel_start.y) + "_" + std::to_string(fluid2_config->vel_start.z)
            + "-phase_" + std::to_string(fluid2_config->phases[0]) + "_" + std::to_string(fluid2_config->phases[1])
            + "-density" + std::to_string(solver_config->phase_rest_density[1])
            + "-bound-density" + std::to_string(solver_config->rest_bound_density)
            + "\\ply";

//    solver_config->export_path = "D:\\xuyuhang\\simulation\\multi-fluid\\render\\ink-dif\\cylinder\\test";
    /**  =============================================================
     * attach objects to solver
     */
    solver->attachObject(fluid);
    solver->attachObject(fluid2);
    solver->attachObject(bound);

    /**  =============================================================
     * run simulation
     */
//    solver->run(5);
    solver->run(20);
}
