//
// Created by ADMIN on 2024/6/13.
//

#include "solvers/PBF/pbf_solver.hpp"


#include "pbf_cuda_api.cuh"
#include "libs/ModelL/model_helper.hpp"
#include "libs/AnalysisL/statistic_util.hpp"

namespace SoSim {
    PBFSolver::PBFSolver() {
        m_config = std::make_shared<PBFSolverConfig>();
        std::cout << "Create PBFSolver.\n";
    }

    PBFSolver::~PBFSolver() {
        destroy();
    }

    std::shared_ptr<SolverConfig> PBFSolver::getConfig() {
        return m_config;
    }

    void PBFSolver::attachObject(std::shared_ptr<Object> object) {
        if (m_objects.count(object) == 0)
            m_objects.insert(object);

        m_change_occur = true;
        std::cout << "PBFSolver attach object: " << object->getName() << ".\n";
    }

    void PBFSolver::attachParticleEmitter(std::shared_ptr<ParticleEmitter> emitter) {
        if (m_emitters.count(emitter) == 0)
            m_emitters.insert(emitter);

        m_change_occur = true;
        std::cout << "PBFSolver attach a ParticleEmitter.\n";
    }

    void PBFSolver::detachObject(std::shared_ptr<Object> object) {
        if (m_objects.count(object) > 0)
            m_objects.erase(object);

        m_change_occur = true;
        std::cout << "PBFSolver detach object: " << object->getName() << ".\n";
    }

    void PBFSolver::detachParticleEmitter(std::shared_ptr<ParticleEmitter> emitter) {
        if (m_emitters.count(emitter) > 0)
            m_emitters.erase(emitter);

        m_change_occur = true;
        std::cout << "PBFSolver detach a ParticleEmitter.\n";
    }

    void PBFSolver::mergeObjects() {
        pos_all.clear();
        vel_all.clear();
        mat_all.clear();

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

                m_host_const.particle_num += static_cast<int>(pos_tmp.size());
                obj_offline.insert(obj);
            }
        }
    }

    bool PBFSolver::initialize() {
        if (!m_config) {
            std::cout << "ERROR:: solver config empty.\n";
            return false;
        }

        if (m_objects.empty() && m_emitters.empty()) {
            std::cout << "ERROR:: solver attach no object.\n";
            return false;
        }

        mergeObjects();

        auto solver_config = dynamic_cast<PBFSolverConfig *>(m_config.get());
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
        m_host_const.rest_density = solver_config->rest_density;
        m_host_const.rest_volume = std::powf(2 * particle_radius, 3);
        m_host_const.rest_rigid_density = solver_config->rest_rigid_density;
        m_host_const.rest_bound_density = solver_config->rest_bound_density;
        m_host_const.sph_h = 4 * particle_radius;
        m_host_const.block_num = solver_config->kernel_blocks;
        m_host_const.thread_num = solver_config->kernel_threads;

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
        cudaMalloc((void **) &m_device_const, sizeof(PBFConstantParams));
        cudaMalloc((void **) &m_device_data, sizeof(PBFDynamicParams));
        m_host_data.malloc(particle_num);
        m_neighborSearch.malloc();

        // TODO data copy
        cudaMemcpy(m_host_data.mat, mat_all.data(), particle_num * sizeof(Material), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.pos, pos_all.data(), particle_num * sizeof(Vec3f), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.vel, vel_all.data(), particle_num * sizeof(Vec3f), cudaMemcpyHostToDevice);

        cudaMemcpy(m_device_const, &m_host_const, sizeof(PBFConstantParams), cudaMemcpyHostToDevice);
        cudaMemcpy(m_device_data, &m_host_data, sizeof(PBFDynamicParams), cudaMemcpyHostToDevice);

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
            std::cout << "PBFSolver initialized.\n";
            m_is_init = true;
            return true;
        }
        return false;
    }

    void PBFSolver::destroy() {
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
                std::cout << "PBFSolver destroyed.\n";
        }
    }

    void PBFSolver::exportAsPly() {
        static int counter = 0;
        static int frame = 1;
        auto config = dynamic_cast<PBFSolverConfig *>(m_config.get());

        static float gap = 1.f / config->export_fps;

        auto part_num = m_host_const.particle_num;
        if (config->export_partial == "fluid")
            part_num = m_unified_part_type_start_index.y;

        if (config->dt * counter >= gap) {
            std::cout << "export index: " << frame << "\n";

            std::vector<Vec3f> pos(part_num);
            std::vector<Vec3f> color(part_num);

            cudaMemcpy(pos.data(),
                       m_host_data.pos,
                       part_num * sizeof(Vec3f),
                       cudaMemcpyDeviceToHost);

            ModelHelper::export3DModelAsPly(pos,
                                            config->export_path.value(),
                                            std::to_string(frame));

            frame++;
            counter = 0;
        }
        counter++;
    }

    void PBFSolver::run(float total_time) {
        if (!m_is_init)
            initialize();

        std::cout << "DFSPH run.\n";

        auto solver_config = dynamic_cast<PBFSolverConfig *>(m_config.get());
        if (m_is_start) {

            m_neighborSearch.update(m_host_data.pos);

            init_data(m_host_const,
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

                if (emitter->emit(solver_config->cur_sim_time))
                    config->insert_index += emitter->getTemplatePartNum();
            }

            step();

            if (m_is_crash)
                break;

            if (dynamic_cast<PBFSolverConfig *>(m_config.get())->export_data &&
                dynamic_cast<PBFSolverConfig *>(m_config.get())->export_path.has_value())
                exportAsPly();

            frame++;
        }
    }

    void PBFSolver::step() {
        auto solver_config = dynamic_cast<PBFSolverConfig *>(m_config.get());

        PBFConstantParams cp;
        cudaMemcpy(&cp,
                   m_device_const,
                   sizeof(PBFConstantParams),
                   cudaMemcpyDeviceToHost);

        apply_ext_force(m_host_const,
                        m_device_const,
                        m_device_data,
                        m_neighborSearch.d_config,
                        m_neighborSearch.d_params);

        m_neighborSearch.update(m_host_data.pos);

        for (int i = 0; i < solver_config->pbf_iter_num; ++i) {
            compute_sph_density_and_error(m_host_const,
                                          m_device_const,
                                          m_device_data,
                                          m_neighborSearch.d_config,
                                          m_neighborSearch.d_params);

            update_lamb(m_host_const,
                        m_device_const,
                        m_device_data,
                        m_neighborSearch.d_config,
                        m_neighborSearch.d_params);

            compute_dx(m_host_const,
                       m_device_const,
                       m_device_data,
                       m_neighborSearch.d_config,
                       m_neighborSearch.d_params);

            apply_dx(m_host_const,
                     m_device_const,
                     m_device_data,
                     m_neighborSearch.d_config,
                     m_neighborSearch.d_params);
        }

        post_correct(m_host_const,
                     m_device_const,
                     m_device_data,
                     m_neighborSearch.d_config,
                     m_neighborSearch.d_params);

        cudaGetLastError();

        solver_config->cur_sim_time += solver_config->dt;
    }

}