//
// Created by ADMIN on 2024/3/26.
//

#include "solvers/IMM-v2/imm_v2_solver.hpp"
#include "imm_v2_cuda_api.cuh"

#include <chrono>

#include "libs/ModelL/model_helper.hpp"
#include "libs/AnalysisL/statistic_util.hpp"


namespace SoSim {
    IMMSolver_v2::IMMSolver_v2() {
        m_config = std::make_shared<IMMSolverConfig_v2>();
        std::cout << "Create IMMSolver_v2.\n";
    }

    IMMSolver_v2::~IMMSolver_v2() {
        destroy();
    }

    std::shared_ptr<SolverConfig> IMMSolver_v2::getConfig() {
        return m_config;
    }

    void IMMSolver_v2::attachObject(std::shared_ptr<Object> object) {
        if (m_objects.count(object) == 0)
            m_objects.insert(object);

        m_change_occur = true;
        std::cout << "IMMSolver_v2 attach object: " << object->getName() << ".\n";
    }

    void IMMSolver_v2::attachParticleEmitter(std::shared_ptr<ParticleEmitter> emitter) {
        if (m_emitters.count(emitter) == 0)
            m_emitters.insert(emitter);

        m_change_occur = true;
        std::cout << "IMMSolver_v2 attach a ParticleEmitter.\n";
    }

    void IMMSolver_v2::detachObject(std::shared_ptr<Object> object) {
        if (m_objects.count(object) > 0)
            m_objects.erase(object);

        m_change_occur = true;
        std::cout << "IMMSolver_v2 detach object: " << object->getName() << ".\n";
    }

    void IMMSolver_v2::detachParticleEmitter(std::shared_ptr<ParticleEmitter> emitter) {
        if (m_emitters.count(emitter) > 0)
            m_emitters.erase(emitter);

        m_change_occur = true;
        std::cout << "IMMSolver_v2 detach a ParticleEmitter.\n";
    }

    void IMMSolver_v2::mergeObjects() {
        pos_all.clear();
        vel_all.clear();
        mat_all.clear();
        vol_frac_all.clear();

        // object push order: COMMON_NEWTON, DYNAMIC_RIGID, FIXED_BOUND
        std::set<std::shared_ptr<Object>> obj_offline;
        m_host_const.particle_num = 0;

        m_unified_part_type_start_index.x = 0;
        for (const auto &obj: m_objects) {
            if (obj->getParticleObjectConfig()->particle_mat.value() == COMMON_NEWTON &&
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

                if (obj->getParticleObjectConfig()->phases.empty())
                    throw std::runtime_error(obj->getName() + " doesn't have phase info, please check.\n");

                auto s = obj->getParticleObjectConfig()->phases;
                std::vector<float> alpha_tmp(pos_tmp.size() * obj->getParticleObjectConfig()->phases.size());
                for (int i = 0; i < pos_tmp.size(); i++) {
                    for (int j = 0; j < obj->getParticleObjectConfig()->phases.size(); j++) {
                        alpha_tmp[i * obj->getParticleObjectConfig()->phases.size() + j] = s[j];
                    }
                }
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

            if (emitter->getConfig()->phases.empty())
                throw std::runtime_error(" This emitter doesn't have phase info, please check.\n");

            auto s = emitter->getConfig()->phases;
            std::vector<float> alpha_tmp(pos_tmp.size() * emitter->getConfig()->phases.size());
            for (int i = 0; i < pos_tmp.size(); i++) {
                for (int j = 0; j < emitter->getConfig()->phases.size(); j++) {
                    alpha_tmp[i * emitter->getConfig()->phases.size() + j] = s[j];
                }
            }
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

                if (obj->getParticleObjectConfig()->phases.empty())
                    throw std::runtime_error(obj->getName() + " doesn't have phase info, please check.\n");

                auto s = obj->getParticleObjectConfig()->phases;
                std::vector<float> alpha_tmp(pos_tmp.size() * obj->getParticleObjectConfig()->phases.size());
                for (int i = 0; i < pos_tmp.size(); i++) {
                    for (int j = 0; j < obj->getParticleObjectConfig()->phases.size(); j++) {
                        alpha_tmp[i * obj->getParticleObjectConfig()->phases.size() + j] = s[j];
                    }
                }
                vol_frac_all.insert(vol_frac_all.end(), alpha_tmp.begin(), alpha_tmp.end());

                m_host_const.particle_num += static_cast<int>(pos_tmp.size());
                obj_offline.insert(obj);
            }
        }

        m_unified_part_type_start_index.z = m_host_const.particle_num;
        for (const auto &obj: m_objects) {
            if (obj->getParticleObjectConfig()->particle_mat.value() != COMMON_NEWTON &&
                obj->getParticleObjectConfig()->particle_mat.value() != DYNAMIC_RIGID &&
                obj_offline.count(obj) < 1) {
                m_obj_start_index.emplace_back(m_host_const.particle_num);

                std::vector<Vec3f> pos_tmp = obj->getParticles();
                pos_all.insert(pos_all.end(), pos_tmp.begin(), pos_tmp.end());

                std::vector<Vec3f> vel_tmp(pos_tmp.size(), obj->getParticleObjectConfig()->vel_start);
                vel_all.insert(vel_all.end(), vel_tmp.begin(), vel_tmp.end());

                std::vector<Material> mat_tmp(pos_tmp.size(), obj->getParticleObjectConfig()->particle_mat.value());
                mat_all.insert(mat_all.end(), mat_tmp.begin(), mat_tmp.end());

                if (obj->getParticleObjectConfig()->phases.empty())
                    throw std::runtime_error(obj->getName() + " doesn't have phase info, please check.\n");

                auto s = obj->getParticleObjectConfig()->phases;
                std::vector<float> alpha_tmp(pos_tmp.size() * obj->getParticleObjectConfig()->phases.size());
                for (int i = 0; i < pos_tmp.size(); i++) {
                    for (int j = 0; j < obj->getParticleObjectConfig()->phases.size(); j++) {
                        alpha_tmp[i * obj->getParticleObjectConfig()->phases.size() + j] = s[j];
                    }
                }
                vol_frac_all.insert(vol_frac_all.end(), alpha_tmp.begin(), alpha_tmp.end());

                m_host_const.particle_num += static_cast<int>(pos_tmp.size());
                obj_offline.insert(obj);
            }
        }

        std::vector<int> phase_nums;
        for (const auto &obj: m_objects){
            std::cout << obj->getName() << " has " << obj->getParticleObjectConfig()->phases.size() << " phases.\n";
            phase_nums.push_back(obj->getParticleObjectConfig()->phases.size());
        }

        for(int i = 1; i < phase_nums.size(); ++i){
            if(phase_nums[i] != phase_nums[i-1]){
                std::cout << "ERROR:: phase number of different objects are not equal.\n";
                m_is_crash = true;
                return;
            }
        }
    }

    bool IMMSolver_v2::initialize() {
        if (!m_config) {
            std::cout << "ERROR:: solver config empty.\n";
            return false;
        }

        if (m_objects.empty() && m_emitters.empty()) {
            std::cout << "ERROR:: solver attach no object.\n";
            return false;
        }

        mergeObjects();

        auto solver_config = dynamic_cast<IMMSolverConfig_v2 *>(m_config.get());
        int device;
        cudaGetDevice(&device);
        cudaDeviceProp prop{};
        cudaGetDeviceProperties(&prop, device);
        solver_config->kernel_threads = prop.maxThreadsPerBlock;
        solver_config->kernel_blocks = std::ceil(
                (m_host_const.particle_num + solver_config->kernel_threads - 1) / solver_config->kernel_threads);

        float particle_radius;
        if (!m_objects.empty())
            particle_radius = m_objects.begin()->get()->getParticleObjectConfig()->particle_radius.value();
        if (!m_emitters.empty())
            particle_radius = m_emitters.begin()->get()->getConfig()->particle_radius.value();
        auto particle_num = m_host_const.particle_num;

        // TODO setup m_host_const
        m_host_const.dt = solver_config->dt;
        m_host_const.inv_dt = 1 / solver_config->dt;
        m_host_const.inv_dt2 = std::powf(m_host_const.inv_dt, 2);
        m_host_const.cur_sim_time = solver_config->cur_sim_time;
        m_host_const.gravity = solver_config->gravity;
        m_host_const.particle_num = particle_num;
        m_host_const.particle_radius = particle_radius;
        m_host_const.rest_volume = std::powf(2 * particle_radius, 3);
        m_host_const.sph_h = 4 * particle_radius;
        m_host_const.Cf = solver_config->Cf;
        m_host_const.Cd = solver_config->Cd;
        m_host_const.div_free_threshold = solver_config->div_free_threshold;
        m_host_const.incompressible_threshold = solver_config->incompressible_threshold;
        m_host_const.block_num = solver_config->kernel_blocks;
        m_host_const.thread_num = solver_config->kernel_threads;

        m_host_const.phase_rest_densities = solver_config->phase_rest_densities;
        m_host_const.phase_colors = solver_config->phase_colors;
        m_host_const.phase_vis = solver_config->phase_vis;
        m_host_const.phase_num = solver_config->phase_colors.size();

        m_host_const.rest_rigid_density = solver_config->rigid_rest_density;
        m_host_const.rest_bound_density = solver_config->bound_rest_density;

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
        cudaMalloc((void **) &m_device_const, sizeof(m_host_const));
        cudaMalloc((void **) &m_device_data, sizeof(m_host_data));
        m_host_data.malloc(particle_num, m_host_const.phase_num);
        m_neighborSearch.malloc();

        // data copy
        cudaMemcpy(m_host_data.mat, mat_all.data(), particle_num * sizeof(Material), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.pos, pos_all.data(), particle_num * sizeof(Vec3f), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.pos_adv, pos_all.data(), particle_num * sizeof(Vec3f), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.vel, vel_all.data(), particle_num * sizeof(Vec3f), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.vel_adv, vel_all.data(), particle_num * sizeof(Vec3f), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.vol_frac, vol_frac_all.data(), particle_num * sizeof(float) * m_host_const.phase_num,
                   cudaMemcpyHostToDevice);

        cudaMemcpy(m_host_data.const_phase_rest_densities, m_host_const.phase_rest_densities.data(),
                   sizeof(float) * m_host_const.phase_num,
                   cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.const_phase_colors, m_host_const.phase_colors.data(),
                   sizeof(Vec3f) * m_host_const.phase_num,
                   cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.const_phase_vis, m_host_const.phase_vis.data(), sizeof(float) * m_host_const.phase_num,
                   cudaMemcpyHostToDevice);

        cudaMemcpy(m_device_const, &m_host_const, sizeof(m_host_const), cudaMemcpyHostToDevice);
        cudaMemcpy(m_device_data, &m_host_data, sizeof(m_host_data), cudaMemcpyHostToDevice);

        // post-init emitter
        for (auto &emitter: m_emitters) {
            auto config = emitter->getConfig();
            if (config->use_unified_buffer) {
                config->unified_pos_buffer = m_host_data.pos;
                config->unified_vel_buffer = m_host_data.vel;
                config->unified_mat_buffer = m_host_data.mat;
            }

            emitter->update();
        }

        if (cudaGetLastError() == cudaSuccess) {
            std::cout << "IMMSolver_v2 initialized.\n";
            m_is_init = true;
            return true;
        }
        return false;
    }

    void IMMSolver_v2::destroy() {
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
                std::cout << "IMMSolver_v2 destroyed.\n";
        }
    }

    void IMMSolver_v2::exportAsPly() {
        static int counter = 0;
        static int frame = 1;
        auto config = dynamic_cast<IMMSolverConfig_v2 *>(m_config.get());

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

    void IMMSolver_v2::syncObjectDeviceJitData() {
        int cnt = 0;
        for (auto &obj: m_objects) {
            auto offset = m_obj_start_index[cnt++];

            cudaMemcpy(obj->m_device_cuda_jit_particles,
                       m_host_data.pos + offset,
                       obj->getParticleNum() * sizeof(Vec3f),
                       cudaMemcpyDeviceToDevice);
        }
    }

    void IMMSolver_v2::run(float total_time) {
        if (!m_is_init)
            initialize();

        std::cout << "IMMSolver_v2 run.\n";

        auto solver_config = dynamic_cast<IMMSolverConfig_v2 *>(m_config.get());
        if (m_is_start) {
            m_neighborSearch.update(m_host_data.pos);

            m_neighborSearch.dump();

            init_data(m_host_const,
                      m_device_const,
                      m_device_data,
                      m_neighborSearch.d_params);

            prepare_ims(m_host_const,
                        m_device_const,
                        m_device_data,
                        m_neighborSearch.d_config,
                        m_neighborSearch.d_params);

            std::vector<Vec3f> pos(m_host_const.particle_num);
            cudaMemcpy(pos.data(), m_host_data.pos,
                       sizeof(Vec3f) * m_host_const.particle_num, cudaMemcpyDeviceToHost);

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

            if (dynamic_cast<IMMSolverConfig_v2 *>(m_config.get())->export_data &&
                dynamic_cast<IMMSolverConfig_v2 *>(m_config.get())->export_path.has_value())
                exportAsPly();

            frame++;
        }
    }

    void IMMSolver_v2::step() {
        auto solver_config = dynamic_cast<IMMSolverConfig_v2 *>(m_config.get());

        auto d_nsConfig = m_neighborSearch.d_config;
        auto d_nsParams = m_neighborSearch.d_params;

        // neighbor search
        m_neighborSearch.update(m_host_data.pos);

        std::vector<Vec3f> pos(m_host_const.particle_num);
        cudaMemcpy(pos.data(), m_host_data.pos,
                   sizeof(Vec3f) * m_host_const.particle_num, cudaMemcpyDeviceToHost);

        sph_precompute(m_host_const,
                       m_device_const,
                       m_device_data,
                       d_nsConfig,
                       d_nsParams);

        std::vector<float> mass(m_host_const.particle_num);
        cudaMemcpy(mass.data(), m_host_data.mass,
                   sizeof(float) * m_host_const.particle_num, cudaMemcpyDeviceToHost);


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

        std::vector<Vec3f> acc(m_host_const.particle_num);
        cudaMemcpy(acc.data(), m_host_data.acc,
                   sizeof(Vec3f) * m_host_const.particle_num, cudaMemcpyDeviceToHost);


        ism_gravity_vis_surface(m_host_const,
                                m_device_const,
                                m_device_data,
                                d_nsConfig,
                                d_nsParams);

        std::vector<Vec3f> vel_phase(m_host_const.particle_num * m_host_const.phase_num);
        cudaMemcpy(vel_phase.data(), m_host_data.vel_phase,
                   m_host_const.phase_num * sizeof(Vec3f) * m_host_const.particle_num, cudaMemcpyDeviceToHost);

        std::vector<float> vol_frac(m_host_const.particle_num * m_host_const.phase_num);
        cudaMemcpy(vol_frac.data(), m_host_data.vol_frac,
                   m_host_const.phase_num * sizeof(float) * m_host_const.particle_num, cudaMemcpyDeviceToHost);

        std::vector<float> cd(m_host_const.particle_num * m_host_const.phase_num);
        cudaMemcpy(cd.data(), m_host_data.Cd,
                   sizeof(float) * m_host_const.particle_num, cudaMemcpyDeviceToHost);

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

        update_pos(m_host_const,
                   m_device_const,
                   m_device_data,
                   d_nsParams);

//        phase_transport_ism(m_host_const,
//                            m_device_const,
//                            m_device_data,
//                            d_nsConfig,
//                            d_nsParams,
//                            m_is_crash);

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