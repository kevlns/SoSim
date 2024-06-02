//
// Created by ADMIN on 2024/4/9.
//
#include "framework/object_manager.hpp"
#include "framework/solver_manager.hpp"

using namespace SoSim;

int main() {

    /**  =============================================================
     * create emitters
     */
    std::shared_ptr<ParticleEmitter> emitter = std::make_shared<ParticleEmitter>();
    auto emitter_config = emitter->getConfig();
    emitter_config->particle_radius = 0.075;
    emitter_config->max_particle_num = 100000;
    emitter_config->use_unified_buffer = true;
    emitter_config->use_emit_agent = true;
    emitter_config->emit_mat = IMSCT_NONNEWTON;
    emitter_config->agent_file = "F:\\DataSet.Research\\BIBM2024\\model\\inject\\rescale\\agent.ply";
    emitter_config->agent_normal_file = "F:\\DataSet.Research\\BIBM2024\\model\\inject\\raw\\agent_normal.ply";
    emitter_config->emit_vel = 8;
    emitter_config->phases.assign({0.5, 0.5});

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
    solver_config->export_path = "F:\\DataSet.Research\\BIBM2024\\ply\\inject\\emitter";
    solver_config->export_fps = 35;
    solver_config->export_partial = "fluid";
    solver_config->export_phase = false;

    solver_config->max_neighborNum = 60;
    solver_config->rest_viscosity = 0.01;
    solver_config->phase1_vis = 0.01;
    solver_config->phase2_vis = 0.05;
    solver_config->rest_density = {980, 1000};
    solver_config->rest_rigid_density = 1000;
    solver_config->rest_bound_density = 1000;
    solver_config->div_free_threshold = 1e-4;
    solver_config->incompressible_threshold = 1e-4;
    solver_config->Cf = 0.25;
    solver_config->phase1_color = {50, 0, 200};
    solver_config->phase2_color = {200, 0, 50};

    solver_config->Cd0 = 0.4;

    solver->attachParticleEmitter(emitter);

    /**  =============================================================
     * run simulation
     */
    solver->run(5);

}
