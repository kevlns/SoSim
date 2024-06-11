//
// Created by ADMIN on 2024/3/26.
//

#include "solvers/IMM-CT/immct_solver.hpp"

#include "immct_cuda_api.cuh"
#include "libs/ModelL/model_helper.hpp"
#include "libs/AnalysisL/statistic_util.hpp"


namespace SoSim {
    IMMCTSolver::IMMCTSolver() {
        m_config = std::make_shared<IMMCTSolverConfig>();
        std::cout << "Create IMMCTSolver.\n";
    }

    IMMCTSolver::~IMMCTSolver() {
        destroy();
    }

    std::shared_ptr<SolverConfig> IMMCTSolver::getConfig() {
        return m_config;
    }

    void IMMCTSolver::attachObject(std::shared_ptr<Object> object) {
        if (m_objects.count(object) == 0)
            m_objects.insert(object);

        m_change_occur = true;
        std::cout << "IMMCTSolver attach object: " << object->getName() << ".\n";
    }

    void IMMCTSolver::attachParticleEmitter(std::shared_ptr<ParticleEmitter> emitter) {
        if (m_emitters.count(emitter) == 0)
            m_emitters.insert(emitter);

        m_change_occur = true;
        std::cout << "IMMCTSolver attach a ParticleEmitter.\n";
    }

    void IMMCTSolver::detachObject(std::shared_ptr<Object> object) {
        if (m_objects.count(object) > 0)
            m_objects.erase(object);

        m_change_occur = true;
        std::cout << "IMMCTSolver detach object: " << object->getName() << ".\n";
    }

    void IMMCTSolver::detachParticleEmitter(std::shared_ptr<ParticleEmitter> emitter) {
        if (m_emitters.count(emitter) > 0)
            m_emitters.erase(emitter);

        m_change_occur = true;
        std::cout << "IMMCTSolver detach a ParticleEmitter.\n";
    }

    void IMMCTSolver::mergeObjects() {
        pos_all.clear();
        vel_all.clear();
        mat_all.clear();
        vol_frac_all.clear();

        // object push order: IMSCT_NONNEWTON, DYNAMIC_RIGID, FIXED_BOUND
        std::set<std::shared_ptr<Object>> obj_offline;
        m_host_const.particle_num = 0;

        m_unified_part_type_start_index.x = 0;
        for (const auto &obj: m_objects) {
            if (obj->getParticleObjectConfig()->particle_mat.value() == IMSCT_NONNEWTON &&
                obj_offline.count(obj) < 1) {

                // TODO: if solver attach order doesn't follow object push order above, this push policy if wrong
                // TODO: also syncObjectDeviceJitData() is wrong
                m_obj_start_index.emplace_back(m_host_const.particle_num);

                std::vector<Vec3f> pos_tmp = obj->getParticles();
                pos_all.insert(pos_all.end(), pos_tmp.begin(), pos_tmp.end());

                std::vector<Vec3f> vel_tmp(pos_tmp.size(), obj->getParticleObjectConfig()->vel_start);
                vel_all.insert(vel_all.end(), vel_tmp.begin(), vel_tmp.end());

                std::vector<Material> mat_tmp(pos_tmp.size(), obj->getParticleObjectConfig()->particle_mat.value());
                mat_all.insert(mat_all.end(), mat_tmp.begin(), mat_tmp.end());

                if (obj->getParticleObjectConfig()->phases.size() != 2)
                    throw std::runtime_error("IMMSolver only solve two-phase fluid now.\n");
                std::vector<Vec2f> alpha_tmp(pos_tmp.size(), {obj->getParticleObjectConfig()->phases[0],
                                                              obj->getParticleObjectConfig()->phases[1]});
                vol_frac_all.insert(vol_frac_all.end(), alpha_tmp.begin(), alpha_tmp.end());

                m_host_const.particle_num += static_cast<int>(pos_tmp.size());
                obj_offline.insert(obj);
            }
        }

        // add emitter particles
        for (const auto &emitter: m_emitters) {
            m_obj_start_index.emplace_back(m_host_const.particle_num);

            auto part_num = emitter->getConfig()->max_particle_num;
            emitter->getConfig()->insert_index = m_host_const.particle_num;

            std::vector<Vec3f> pos_tmp(part_num);
            pos_all.insert(pos_all.end(), pos_tmp.begin(), pos_tmp.end());

            std::vector<Vec3f> vel_tmp(part_num);
            vel_all.insert(vel_all.end(), vel_tmp.begin(), vel_tmp.end());

            std::vector<Material> mat_tmp(part_num, Emitter_Particle);
            mat_all.insert(mat_all.end(), mat_tmp.begin(), mat_tmp.end());

            if (emitter->getConfig()->phases.size() != 2)
                throw std::runtime_error("IMMSolver only solve two-phase fluid now.\n");
            std::vector<Vec2f> alpha_tmp(part_num, {emitter->getConfig()->phases[0],
                                                    emitter->getConfig()->phases[1]});
            vol_frac_all.insert(vol_frac_all.end(), alpha_tmp.begin(), alpha_tmp.end());

            m_host_const.particle_num += static_cast<int>(part_num);
        }

        m_unified_part_type_start_index.y = m_host_const.particle_num;
        for (const auto &obj: m_objects) {
            if (obj->getParticleObjectConfig()->particle_mat.value() == DYNAMIC_RIGID &&
                obj_offline.count(obj) < 1) {
                m_obj_start_index.emplace_back(m_host_const.particle_num);

                std::vector<Vec3f> pos_tmp = obj->getParticles();
                pos_all.insert(pos_all.end(), pos_tmp.begin(), pos_tmp.end());

                std::vector<Vec3f> vel_tmp(pos_tmp.size(), obj->getParticleObjectConfig()->vel_start);
                vel_all.insert(vel_all.end(), vel_tmp.begin(), vel_tmp.end());

                std::vector<Material> mat_tmp(pos_tmp.size(), obj->getParticleObjectConfig()->particle_mat.value());
                mat_all.insert(mat_all.end(), mat_tmp.begin(), mat_tmp.end());

                if (obj->getParticleObjectConfig()->phases.size() != 2)
                    throw std::runtime_error("IMMSolver only solve two-phase fluid now.\n");
                std::vector<Vec2f> alpha_tmp(pos_tmp.size(), {obj->getParticleObjectConfig()->phases[0],
                                                              obj->getParticleObjectConfig()->phases[1]});
                vol_frac_all.insert(vol_frac_all.end(), alpha_tmp.begin(), alpha_tmp.end());

                m_host_const.particle_num += static_cast<int>(pos_tmp.size());
                obj_offline.insert(obj);
            }
        }

        m_unified_part_type_start_index.z = m_host_const.particle_num;
        for (const auto &obj: m_objects) {
            if (obj->getParticleObjectConfig()->particle_mat.value() != IMSCT_NONNEWTON &&
                obj->getParticleObjectConfig()->particle_mat.value() != DYNAMIC_RIGID &&
                obj_offline.count(obj) < 1) {
                m_obj_start_index.emplace_back(m_host_const.particle_num);

                std::vector<Vec3f> pos_tmp = obj->getParticles();
                pos_all.insert(pos_all.end(), pos_tmp.begin(), pos_tmp.end());

                std::vector<Vec3f> vel_tmp(pos_tmp.size(), obj->getParticleObjectConfig()->vel_start);
                vel_all.insert(vel_all.end(), vel_tmp.begin(), vel_tmp.end());

                std::vector<Material> mat_tmp(pos_tmp.size(), obj->getParticleObjectConfig()->particle_mat.value());
                mat_all.insert(mat_all.end(), mat_tmp.begin(), mat_tmp.end());

                if (obj->getParticleObjectConfig()->phases.size() != 2)
                    throw std::runtime_error("IMMSolver only solve two-phase fluid now.\n");
                std::vector<Vec2f> alpha_tmp(pos_tmp.size(), {obj->getParticleObjectConfig()->phases[0],
                                                              obj->getParticleObjectConfig()->phases[1]});
                vol_frac_all.insert(vol_frac_all.end(), alpha_tmp.begin(), alpha_tmp.end());

                m_host_const.particle_num += static_cast<int>(pos_tmp.size());
                obj_offline.insert(obj);
            }
        }
    }

    bool IMMCTSolver::initialize() {
        if (!m_config) {
            std::cout << "ERROR:: solver config empty.\n";
            return false;
        }

        if (m_objects.empty()) {
            std::cout << "ERROR:: solver attach no object.\n";
            return false;
        }

        mergeObjects();

        auto solver_config = dynamic_cast<IMMCTSolverConfig *>(m_config.get());
        int device;
        cudaGetDevice(&device);
        cudaDeviceProp prop{};
        cudaGetDeviceProperties(&prop, device);
        solver_config->kernel_threads = prop.maxThreadsPerBlock;
        solver_config->kernel_blocks = std::ceil(
                (m_host_const.particle_num + solver_config->kernel_threads - 1) / solver_config->kernel_threads);

        auto particle_radius = m_objects.begin()->get()->getParticleObjectConfig()->particle_radius.value();
        auto particle_num = m_host_const.particle_num;

        // TODO setup m_host_const
        m_host_const.dt = solver_config->dt;
        m_host_const.inv_dt = 1 / solver_config->dt;
        m_host_const.inv_dt2 = std::powf(m_host_const.inv_dt, 2);
        m_host_const.cur_sim_time = solver_config->cur_sim_time;
        m_host_const.gravity = solver_config->gravity;
        m_host_const.particle_num = particle_num;
        m_host_const.particle_radius = particle_radius;
        m_host_const.phase1_color = solver_config->phase1_color;
        m_host_const.phase2_color = solver_config->phase2_color;
        m_host_const.rest_density = solver_config->rest_density;
        m_host_const.rest_volume = std::powf(2 * particle_radius, 3);
        m_host_const.rest_rigid_density = solver_config->rest_rigid_density;
        m_host_const.rest_bound_density = solver_config->rest_bound_density;
        m_host_const.sph_h = 4 * particle_radius;
        m_host_const.rest_viscosity = solver_config->rest_viscosity;
        m_host_const.Cf = solver_config->Cf;
        m_host_const.Cd = solver_config->Cd;
        m_host_const.div_free_threshold = solver_config->div_free_threshold;
        m_host_const.incompressible_threshold = solver_config->incompressible_threshold;
        m_host_const.block_num = solver_config->kernel_blocks;
        m_host_const.thread_num = solver_config->kernel_threads;

        m_host_const.Cd0 = solver_config->Cd0;
        m_host_const.ct_thinning_exp0 = solver_config->ct_thinning_exp0;
        m_host_const.ct_relaxation_time = solver_config->ct_relaxation_time;
        m_host_const.solution_vis_base = solver_config->solution_vis_base;
        m_host_const.solution_vis_max = solver_config->solution_vis_max;
        m_host_const.polymer_vol_frac0 = solver_config->polymer_vol_frac0;

        m_host_const.phase1_vis = solver_config->phase1_vis;
        m_host_const.phase2_vis = solver_config->phase2_vis;

        // setup neighbor search
        NeighborSearchUGConfig ns_config;
        ns_config.sceneLB = solver_config->scene_lb;
        ns_config.sceneSize = solver_config->scene_size;
        ns_config.cellLength = m_host_const.sph_h;
        ns_config.particle_num = particle_num;
        ns_config.maxNeighborNum = solver_config->max_neighborNum;
        ns_config.kernel_threads = solver_config->kernel_threads;
        ns_config.kernel_blocks = solver_config->kernel_blocks;
        m_neighborSearch.setConfig(ns_config);

        // malloc
        cudaMalloc((void **) &m_device_const, sizeof(IMMCTConstantParams));
        cudaMalloc((void **) &m_device_data, sizeof(IMMCTDynamicParams));
        m_host_data.malloc(particle_num);
        m_neighborSearch.malloc();

        // TODO data copy
        cudaMemcpy(m_host_data.mat, mat_all.data(), particle_num * sizeof(Material), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.pos, pos_all.data(), particle_num * sizeof(Vec3f), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.pos_adv, pos_all.data(), particle_num * sizeof(Vec3f), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.vel, vel_all.data(), particle_num * sizeof(Vec3f), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.vel_adv, vel_all.data(), particle_num * sizeof(Vec3f), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.vel_phase_1, vel_all.data(), particle_num * sizeof(Vec3f), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.vel_phase_2, vel_all.data(), particle_num * sizeof(Vec3f), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.vol_frac, vol_frac_all.data(), particle_num * sizeof(Vec2f), cudaMemcpyHostToDevice);

        cudaMemcpy(m_device_const, &m_host_const, sizeof(IMMCTConstantParams), cudaMemcpyHostToDevice);
        cudaMemcpy(m_device_data, &m_host_data, sizeof(IMMCTDynamicParams), cudaMemcpyHostToDevice);

        if (cudaGetLastError() == cudaSuccess) {
            std::cout << "IMMCTSolver initialized.\n";
            m_is_init = true;
            return true;
        }
        return false;
    }

    void IMMCTSolver::destroy() {
        m_objects.clear();

        if (m_is_init) {
            // delete neighbor_search
            cudaFree(m_device_const);

            // delete data
            m_host_data.freeMemory();
            cudaFree(m_device_data);

            // delete const
            m_neighborSearch.freeMemory();

            if (cudaGetLastError() == cudaSuccess)
                std::cout << "IMMCTSolver destroyed.\n";
        }
    }

    void IMMCTSolver::exportAsPly() {
        static int counter = 0;
        static int frame = 1;
        auto config = dynamic_cast<IMMCTSolverConfig *>(m_config.get());

        static float gap = 1.f / config->export_fps;

        auto part_num = m_host_const.particle_num;
        if (config->export_partial == "fluid")
            part_num = m_unified_part_type_start_index.y;

        if (config->dt * counter >= gap) {
            std::cout << "export index: " << frame << "\n";

            std::vector<Vec3f> pos(part_num);
            std::vector<Vec3f> color(part_num);
            std::vector<Vec2f> phase(part_num);
            if (config->export_phase)
                cudaMemcpy(phase.data(),
                           m_host_data.vol_frac,
                           part_num * sizeof(Vec2f),
                           cudaMemcpyDeviceToHost);
            cudaMemcpy(pos.data(),
                       m_host_data.pos,
                       part_num * sizeof(Vec3f),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(color.data(),
                       m_host_data.color,
                       part_num * sizeof(Vec3f),
                       cudaMemcpyDeviceToHost);

            if (config->export_phase)
                ModelHelper::export3DModelAsPly(pos,
                                                phase,
                                                config->export_path.value() + "_phase",
                                                std::to_string(frame));
            ModelHelper::export3DModelAsPly(pos,
                                            color,
                                            config->export_path.value(),
                                            std::to_string(frame));

            frame++;
            counter = 0;
        }
        counter++;
    }

    void IMMCTSolver::syncObjectDeviceJitData() {
        int cnt = 0;
        for (auto &obj: m_objects) {
            auto offset = m_obj_start_index[cnt++];

            cudaMemcpy(obj->m_device_cuda_jit_particles,
                       m_host_data.pos + offset,
                       obj->getParticleNum() * sizeof(Vec3f),
                       cudaMemcpyDeviceToDevice);
        }
    }

    void IMMCTSolver::run(float total_time) {
        if (!m_is_init)
            initialize();

        std::cout << "IMMCTSolver run.\n";

        auto solver_config = dynamic_cast<IMMCTSolverConfig *>(m_config.get());
        if (m_is_start) {

            m_neighborSearch.update(m_host_data.pos);

            init_data(m_host_const,
                      m_device_const,
                      m_device_data,
                      m_neighborSearch.d_params);

            prepare_ims(m_host_const,
                        m_device_const,
                        m_device_data,
                        m_neighborSearch.d_config,
                        m_neighborSearch.d_params);

            m_is_start = false;
        }


        int frame = 1;
        while (solver_config->cur_sim_time < total_time) {

            std::cout << "*========== Frame: " << frame << " ==========* \n";

            for (const auto &emitter: m_emitters) {
                auto config = emitter->getConfig();
                if (config->start_time.has_value() && config->start_time.value() > solver_config->cur_sim_time)
                    continue;
                if (config->end_time.has_value() && config->end_time.value() < solver_config->cur_sim_time)
                    continue;

                if (emitter->emit(solver_config->cur_sim_time)) {
                    auto size_1 = emitter->getTemplatePartNum() * sizeof(Vec3f);
                    // copy vel to phase vel
                    cudaMemcpy(m_host_data.vel_adv + config->insert_index,
                               emitter->getCudaTemplateVelBuffer(), size_1, cudaMemcpyDeviceToDevice);
                    cudaMemcpy(m_host_data.pos_adv + config->insert_index,
                               emitter->getCudaTemplatePosBuffer(), size_1, cudaMemcpyHostToDevice);

                    config->insert_index += emitter->getTemplatePartNum();
                }
            }

            step();

            if (m_is_crash)
                break;

            if (dynamic_cast<IMMCTSolverConfig *>(m_config.get())->export_data &&
                dynamic_cast<IMMCTSolverConfig *>(m_config.get())->export_path.has_value())
                exportAsPly();

            frame++;
        }
    }

    void IMMCTSolver::step() {
        auto solver_config = dynamic_cast<IMMCTSolverConfig *>(m_config.get());

        auto d_nsConfig = m_neighborSearch.d_config;
        auto d_nsParams = m_neighborSearch.d_params;

        // neighbor search
        m_neighborSearch.update(m_host_data.pos);

        sph_precompute(m_host_const,
                       m_device_const,
                       m_device_data,
                       d_nsConfig,
                       d_nsParams);

        update_CT_parameters(m_host_const,
                             m_device_const,
                             m_device_data,
                             d_nsConfig,
                             d_nsParams);

        vfsph_div(m_host_const,
                  m_host_data,
                  m_unified_part_type_start_index,
                  m_device_const,
                  m_device_data,
                  d_nsConfig,
                  d_nsParams,
                  m_is_crash);

        apply_pressure_acc(m_host_const,
                           m_device_const,
                           m_device_data,
                           d_nsParams);

        ism_gravity_vis_surface(m_host_const,
                                m_device_const,
                                m_device_data,
                                d_nsConfig,
                                d_nsParams);

        vfsph_incomp(m_host_const,
                     m_host_data,
                     m_unified_part_type_start_index,
                     m_device_const,
                     m_device_data,
                     d_nsConfig,
                     d_nsParams,
                     m_is_crash);

        apply_pressure_acc(m_host_const,
                           m_device_const,
                           m_device_data,
                           d_nsParams);

        ism_viscoelastic(m_host_const,
                         m_device_const,
                         m_device_data,
                         d_nsConfig,
                         d_nsParams);

        artificial_vis_bound(m_host_const,
                             m_device_const,
                             m_device_data,
                             d_nsConfig,
                             d_nsParams);

        update_pos(m_host_const,
                   m_device_const,
                   m_device_data,
                   d_nsParams);

        phase_transport_ism(m_host_const,
                            m_device_const,
                            m_device_data,
                            d_nsConfig,
                            d_nsParams,
                            m_is_crash);

        update_mass_and_vel(m_host_const,
                            m_device_const,
                            m_device_data,
                            d_nsParams);

        update_color(m_host_const,
                     m_device_const,
                     m_device_data,
                     d_nsParams);

        syncObjectDeviceJitData();

        cudaGetLastError();

        solver_config->cur_sim_time += solver_config->dt;
    }
}