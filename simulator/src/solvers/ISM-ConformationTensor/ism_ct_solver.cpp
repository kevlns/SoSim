//
// Created by ADMIN on 2024/3/26.
//

#include "solvers/ISM-ConformationTensor/ism_ct_solver.hpp"
#include "ism_ct_cuda_api.cuh"

#include "libs/ModelL/model_helper.hpp"
#include "libs/AnalysisL/statistic_util.hpp"


namespace SoSim {
    IMSCTSolver::IMSCTSolver() {
        m_config = std::make_shared<IMSCTSolverConfig>();
        std::cout << "Create ISMCTSolver.\n";
    }

    IMSCTSolver::~IMSCTSolver() {
        destroy();
    }

    std::shared_ptr<SolverConfig> IMSCTSolver::getConfig() {
        return m_config;
    }

    void IMSCTSolver::attachObject(std::shared_ptr<Object> object) {
        if (m_objects.count(object) == 0)
            m_objects.insert(object);

        m_change_occur = true;
        std::cout << "ISMCTSolver attach object: " << object->getName() << ".\n";
    }

    void IMSCTSolver::detachObject(std::shared_ptr<Object> object) {
        if (m_objects.count(object) > 0)
            m_objects.erase(object);

        m_change_occur = true;
        std::cout << "ISMCTSolver detach object: " << object->getName() << ".\n";
    }

    void IMSCTSolver::mergeObjects() {
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
                    throw std::runtime_error("ISMCTSolver only solve two-phase fluid now.\n");
                std::vector<Vec2f> alpha_tmp(pos_tmp.size(), {obj->getParticleObjectConfig()->phases[0],
                                                              obj->getParticleObjectConfig()->phases[1]});
                vol_frac_all.insert(vol_frac_all.end(), alpha_tmp.begin(), alpha_tmp.end());

                m_host_const.particle_num += static_cast<int>(pos_tmp.size());
                obj_offline.insert(obj);
            }
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
                    throw std::runtime_error("ISMCTSolver only solve two-phase fluid now.\n");
                std::vector<Vec2f> alpha_tmp(pos_tmp.size(), {obj->getParticleObjectConfig()->phases[0],
                                                              obj->getParticleObjectConfig()->phases[1]});
                vol_frac_all.insert(vol_frac_all.end(), alpha_tmp.begin(), alpha_tmp.end());

                m_host_const.particle_num += static_cast<int>(pos_tmp.size());
                obj_offline.insert(obj);
            }
        }

        m_unified_part_type_start_index.z = m_host_const.particle_num;
        for (const auto &obj: m_objects) {
            if (obj->getParticleObjectConfig()->particle_mat.value() == FIXED_BOUND &&
                obj_offline.count(obj) < 1) {
                m_obj_start_index.emplace_back(m_host_const.particle_num);

                std::vector<Vec3f> pos_tmp = obj->getParticles();
                pos_all.insert(pos_all.end(), pos_tmp.begin(), pos_tmp.end());

                std::vector<Vec3f> vel_tmp(pos_tmp.size(), obj->getParticleObjectConfig()->vel_start);
                vel_all.insert(vel_all.end(), vel_tmp.begin(), vel_tmp.end());

                std::vector<Material> mat_tmp(pos_tmp.size(), obj->getParticleObjectConfig()->particle_mat.value());
                mat_all.insert(mat_all.end(), mat_tmp.begin(), mat_tmp.end());

                if (obj->getParticleObjectConfig()->phases.size() != 2)
                    throw std::runtime_error("ISMCTSolver only solve two-phase fluid now.\n");
                std::vector<Vec2f> alpha_tmp(pos_tmp.size(), {obj->getParticleObjectConfig()->phases[0],
                                                              obj->getParticleObjectConfig()->phases[1]});
                vol_frac_all.insert(vol_frac_all.end(), alpha_tmp.begin(), alpha_tmp.end());

                m_host_const.particle_num += static_cast<int>(pos_tmp.size());
                obj_offline.insert(obj);
            }
        }
    }

    bool IMSCTSolver::initialize() {
        if (!m_config) {
            std::cout << "ERROR:: solver config empty.\n";
            return false;
        }

        if (m_objects.empty()) {
            std::cout << "ERROR:: solver attach no object.\n";
            return false;
        }

        mergeObjects();

        auto solver_config = dynamic_cast<IMSCTSolverConfig *>(m_config.get());
        int device;
        cudaGetDevice(&device);
        cudaDeviceProp prop{};
        cudaGetDeviceProperties(&prop, device);
        solver_config->kernel_threads = prop.maxThreadsPerBlock;
        solver_config->kernel_blocks = std::ceil(
                (m_host_const.particle_num + solver_config->kernel_threads - 1) / solver_config->kernel_threads);

        auto particle_radius = m_objects.begin()->get()->getParticleObjectConfig()->particle_radius;
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
        cudaMalloc((void **) &m_device_const, sizeof(IMSCTConstantParams));
        cudaMalloc((void **) &m_device_data, sizeof(IMSCTDynamicParams));
        m_host_data.malloc(particle_num);
        m_neighborSearch.malloc();

        // TODO data copy
        cudaMemcpy(m_host_data.mat, mat_all.data(), particle_num * sizeof(Material), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.pos, pos_all.data(), particle_num * sizeof(Vec3f), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.pos_adv, pos_all.data(), particle_num * sizeof(Vec3f), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.vel, vel_all.data(), particle_num * sizeof(Vec3f), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.vel_adv, vel_all.data(), particle_num * sizeof(Vec3f), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.vol_frac, vol_frac_all.data(), particle_num * sizeof(Vec2f), cudaMemcpyHostToDevice);

        cudaMemcpy(m_device_const, &m_host_const, sizeof(IMSCTConstantParams), cudaMemcpyHostToDevice);
        cudaMemcpy(m_device_data, &m_host_data, sizeof(IMSCTDynamicParams), cudaMemcpyHostToDevice);

        if (cudaGetLastError() == cudaSuccess) {
            std::cout << "ISMCTSolver initialized.\n";
            m_is_init = true;
            return true;
        }
        return false;
    }

    void IMSCTSolver::destroy() {
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
                std::cout << "ISMCTSolver destroyed.\n";
        }
    }

    void IMSCTSolver::exportAsPly() {
        static int counter = 1;
        std::vector<Vec3f> pos(m_unified_part_type_start_index.y);
        std::vector<Vec3f> color(m_unified_part_type_start_index.y);
        cudaMemcpy(pos.data(),
                   m_host_data.pos,
                   m_unified_part_type_start_index.y * sizeof(Vec3f),
                   cudaMemcpyDeviceToHost);
        cudaMemcpy(color.data(),
                   m_host_data.color,
                   m_unified_part_type_start_index.y * sizeof(Vec3f),
                   cudaMemcpyDeviceToHost);
        ModelHelper::export3DModelAsPly(pos,
                                        color,
                                        "F:\\DataSet.Research\\ITEM.NN_NewModel\\ply\\imsct_test\\phase_trans",
                                        std::to_string(counter++));
    }

    void IMSCTSolver::syncObjectDeviceJitData() {
        int cnt = 0;
        for (auto &obj: m_objects) {
            auto offset = m_obj_start_index[cnt++];

            cudaMemcpy(obj->m_device_cuda_jit_particles,
                       m_host_data.pos + offset,
                       obj->getParticleNum() * sizeof(Vec3f),
                       cudaMemcpyDeviceToDevice);
        }
    }

    void IMSCTSolver::run(float total_time) {
        if (!m_is_init)
            initialize();

        std::cout << "IMSCT run.\n";

        auto solver_config = dynamic_cast<IMSCTSolverConfig *>(m_config.get());
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

            step();

            if (dynamic_cast<IMSCTSolverConfig *>(m_config.get())->export_data)
                exportAsPly();

            frame++;
        }
    }

    void IMSCTSolver::step() {
        auto solver_config = dynamic_cast<IMSCTSolverConfig *>(m_config.get());

        auto d_nsConfig = m_neighborSearch.d_config;
        auto d_nsParams = m_neighborSearch.d_params;

        // update color

        // neighbor search
        m_neighborSearch.update(m_host_data.pos);

        sph_precompute(m_host_const,
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
                  d_nsParams);

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
                     d_nsParams);

        apply_pressure_acc(m_host_const,
                           m_device_const,
                           m_device_data,
                           d_nsParams);
//
        update_pos(m_host_const,
                   m_device_const,
                   m_device_data,
                   d_nsParams);

        phase_transport_ism(m_host_const,
                            m_device_const,
                            m_device_data,
                            d_nsConfig,
                            d_nsParams);

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