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
    std::shared_ptr<ParticleEmitter> emitter_1 = std::make_shared<ParticleEmitter>();
    auto emitter1_config = emitter_1->getConfig();
    emitter1_config->particle_radius = 0.075;
    emitter1_config->max_particle_num = 100000;
    emitter1_config->use_unified_buffer = true;
    emitter1_config->use_emit_agent = true;
    emitter1_config->emit_mat = IMSCT_NONNEWTON;
    emitter1_config->agent_file = "F:\\DataSet.Research\\CMPB\\models\\agent1_parts.ply";
    emitter1_config->agent_normal_file = "F:\\DataSet.Research\\CMPB\\models\\agent1_normal.ply";
    emitter1_config->emit_vel = 0.65;
    emitter1_config->phases.assign({0.45, 0.55});

    std::shared_ptr<ParticleEmitter> emitter_2 = std::make_shared<ParticleEmitter>();
    auto emitter2_config = emitter_2->getConfig();
    emitter2_config->particle_radius = 0.075;
    emitter2_config->max_particle_num = 100000;
    emitter2_config->use_unified_buffer = true;
    emitter2_config->use_emit_agent = true;
    emitter2_config->emit_mat = IMSCT_NONNEWTON;
    emitter2_config->agent_file = "F:\\DataSet.Research\\CMPB\\models\\agent2_parts.ply";
    emitter2_config->agent_normal_file = "F:\\DataSet.Research\\CMPB\\models\\agent2_normal.ply";
    emitter2_config->emit_vel = 0.65;
    emitter2_config->phases.assign({0.45, 0.55});

    std::shared_ptr<ParticleEmitter> emitter_3 = std::make_shared<ParticleEmitter>();
    auto emitter3_config = emitter_3->getConfig();
    emitter3_config->particle_radius = 0.075;
    emitter3_config->max_particle_num = 100000;
    emitter3_config->use_unified_buffer = true;
    emitter3_config->use_emit_agent = true;
    emitter3_config->emit_mat = IMSCT_NONNEWTON;
    emitter3_config->agent_file = "F:\\DataSet.Research\\CMPB\\models\\agent3_parts.ply";
    emitter3_config->agent_normal_file = "F:\\DataSet.Research\\CMPB\\models\\agent3_normal.ply";
    emitter3_config->emit_vel = 0.65;
    emitter3_config->phases.assign({0.45, 0.55});


    ObjectManager objectManager;
    auto vertebra = objectManager.createObject();
    auto vertebra_config = vertebra->getParticleObjectConfig();
    vertebra_config->particle_radius = 0.075;
    vertebra_config->particle_mat = FIXED_BOUND;
    vertebra_config->phases.assign({0, 0});
    vertebra_config->model_file = "F:\\DataSet.Research\\CMPB\\models\\inner-surface.ply";
    vertebra->setName("inner");
    vertebra->update();

    auto vertebra_cover = objectManager.createObject();
    auto vertebra_cover_config = vertebra_cover->getParticleObjectConfig();
    vertebra_cover_config->particle_radius = 0.075;
    vertebra_cover_config->particle_mat = FIXED_BOUND;
    vertebra_cover_config->phases.assign({0, 0});
    vertebra_cover_config->model_file = "F:\\DataSet.Research\\CMPB\\models\\bound.ply";
    vertebra_cover->setName("bound");
    vertebra_cover->update();

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
    solver_config->dt = 0.001;
    solver_config->gravity = {0, -2, 0};
    solver_config->scene_lb = {-12, -5, -15};
    solver_config->scene_size = {24, 25, 30};
    solver_config->export_data = true;
    solver_config->export_path = "F:\\DataSet.Research\\CMPB\\res_ply\\emitter_1";
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

    solver_config->Cd0 = 0.6;
    solver_config->ct_thinning_exp0 = 0.6;
    solver_config->ct_relaxation_time = 0.0025;
    solver_config->solution_vis_base = 8000;
    solver_config->solution_vis_max = 9000;
    solver_config->polymer_vol_frac0 = 0.85; // threshold

//    auto world = objectManager.createObject();
//    auto world_config = world->getParticleObjectConfig();
//    world_config->particle_radius = 0.075;
//    world_config->particle_mat = FIXED_BOUND;
//    world_config->shape = Box;
//    world_config->lb = solver_config->scene_lb;
//    world_config->size = solver_config->scene_size;
//    world_config->layer = 1;
//    world->update();
//    world->exportAsPly("F:\\DataSet.Research\\CMPB\\res_ply\\bound",
//                       "world_bound");

    solver->attachParticleEmitter(emitter_1);
    solver->attachParticleEmitter(emitter_2);
    solver->attachParticleEmitter(emitter_3);
    solver->attachObject(vertebra);
    solver->attachObject(vertebra_cover);
//    solver->attachObject(network);

    /**  =============================================================
     * run simulation
     */
    solver->run(15);

}