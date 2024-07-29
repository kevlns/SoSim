//
// Created by ADMIN on 2024/6/15.
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
    emitter_config->emit_vel = 0.8;
    emitter_config->phases.assign({0.25, 0.75});

    ObjectManager objectManager;
    auto fluid = objectManager.createObject();
    auto fluid_config = fluid->getParticleObjectConfig();
    fluid_config->particle_radius = 0.075;
    fluid_config->particle_mat = IMSCT_NONNEWTON;
    fluid_config->phases.assign({0.15, 0.85});
    fluid_config->model_file = "F:\\DataSet.Research\\BIBM2024\\model\\inject\\rescale\\fluid_cylinder.ply";
    fluid->setName("fluid");
    fluid->update();

    auto vertebra = objectManager.createObject();
    auto vertebra_config = vertebra->getParticleObjectConfig();
    vertebra_config->particle_radius = 0.075;
    vertebra_config->particle_mat = FIXED_BOUND;
    vertebra_config->phases.assign({0, 0});
    vertebra_config->model_file = "F:\\DataSet.Research\\BIBM2024\\model\\inject\\rescale\\vertebra_tube.ply";
    vertebra->setName("vertebra");
    vertebra->update();

    auto vertebra_cover = objectManager.createObject();
    auto vertebra_cover_config = vertebra_cover->getParticleObjectConfig();
    vertebra_cover_config->particle_radius = 0.075;
    vertebra_cover_config->particle_mat = FIXED_BOUND;
    vertebra_cover_config->phases.assign({0, 0});
    vertebra_cover_config->model_file = "F:\\DataSet.Research\\BIBM2024\\model\\inject\\rescale\\box.ply";
    vertebra_cover->setName("vertebra_cover");
    vertebra_cover->update();
//
//    auto network = objectManager.createObject();
//    auto network_config = network->getParticleObjectConfig();
//    network_config->particle_radius = 0.075;
//    network_config->particle_mat = FIXED_BOUND;
//    network_config->phases.assign({0, 0});
//    network_config->model_file = "F:\\DataSet.Research\\BIBM2024\\model\\inject\\rescale\\network.ply";
//    network->setName("network");
//    network->update();

    /**  =============================================================
     * create solvers
     */
    SolverManager solverManager;
    auto solver = solverManager.createSolver<IMMCTSolver>();
    auto solver_config = dynamic_cast<IMMCTSolverConfig *>(solver->getConfig().get());
    solver_config->dt = 0.0005;
    solver_config->gravity = {0, -3, 0};
    solver_config->scene_lb = {-25, -40, -25};
    solver_config->scene_size = {50, 160, 50};
    solver_config->export_data = true;
    solver_config->export_path = "F:\\DataSet.Research\\BIBM2024\\ply\\inject\\new\\0.15_0.85";
    solver_config->export_fps = 35;
    solver_config->export_partial = "fluid";
    solver_config->export_phase = false;

    solver_config->max_neighborNum = 60;
    solver_config->rest_viscosity = 0.01;
    solver_config->rest_density = {980, 1000};
    solver_config->rest_rigid_density = 1000;
    solver_config->rest_bound_density = 1000;
    solver_config->div_free_threshold = 1e-4;
    solver_config->incompressible_threshold = 1e-4;
    solver_config->Cf = 0.1;
    solver_config->phase1_color = {50, 0, 200};
    solver_config->phase2_color = {200, 0, 50};

    solver_config->Cd0 = 0.5;
    solver_config->ct_thinning_exp0 = 0.65;
    solver_config->ct_relaxation_time = 0.005;
    solver_config->solution_vis_base = 7000;
    solver_config->solution_vis_max = 10000;
    solver_config->polymer_vol_frac0 = 0.65; // threshold
    solver_config->vis_bound_damp_factor = 0.9;

    solver->attachObject(fluid);
    solver->attachObject(vertebra);
    solver->attachObject(vertebra_cover);
//    solver->attachObject(network);

//    auto world = objectManager.createObject();
//    auto world_config = world->getParticleObjectConfig();
//    world_config->particle_radius = 0.075;
//    world_config->particle_mat = FIXED_BOUND;
//    world_config->shape = Box;
//    world_config->lb = solver_config->scene_lb;
//    world_config->size = solver_config->scene_size;
//    world_config->layer = 1;
//    world->update();
//    world->exportAsPly("F:\\DataSet.Research\\BIBM2024\\model\\bound",
//                       "world_bound");

    /**  =============================================================
     * run simulation
     */
    solver->run(10);

}