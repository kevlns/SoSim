//
// Created by ADMIN on 2024/3/22.
//

#include "solvers/JL21-ConformationTensor/jl21_ct_solver.hpp"
#include "jl21_ct_cuda_api.cuh"
//#include "jl21_solver_host_invoke_api.cuh"

#include "libs/ModelL/model_exporter.hpp"
#include "libs/AnalysisL/statistic_util.hpp"

namespace SoSim {

    JL21CTSolver::JL21CTSolver() {
        m_config = std::make_shared<JL21CTSolverConfig>();
        std::cout << "Create JL21CTSolver.\n";
    }

    JL21CTSolver::~JL21CTSolver() {
        destroy();
    }

    std::shared_ptr<SolverConfig> JL21CTSolver::getConfig() {
        return m_config;
    }

    void JL21CTSolver::attachObject(std::shared_ptr<Object> object) {
        if (m_objects.count(object) == 0)
            m_objects.insert(object);

        m_change_occur = true;
        std::cout << "JL21CTSolver attach object: " << object->getName() << ".\n";
    }

    void JL21CTSolver::detachObject(std::shared_ptr<Object> object) {
        if (m_objects.count(object) > 0)
            m_objects.erase(object);

        m_change_occur = true;
        std::cout << "JL21CTSolver detach object: " << object->getName() << ".\n";
    }

    void JL21CTSolver::mergeObjects() {
        pos_all.clear();
        vel_all.clear();
        mat_all.clear();
        alpha_all.clear();

        m_host_const.total_particle_num = 0;
        for (const auto &obj: m_objects) {
            std::vector<Vec3f> pos_tmp = obj->getParticles();
            pos_all.insert(pos_all.end(), pos_tmp.begin(), pos_tmp.end());

            std::vector<Vec3f> vel_tmp(pos_tmp.size(), obj->getParticleObjectConfig()->vel_start);
            vel_all.insert(vel_all.end(), vel_tmp.begin(), vel_tmp.end());

            std::vector<Material> mat_tmp(pos_tmp.size(), obj->getParticleObjectConfig()->particle_mat.value());
            mat_all.insert(mat_all.end(), mat_tmp.begin(), mat_tmp.end());

            std::vector<Vec2f> alpha_tmp(pos_tmp.size(), {obj->getParticleObjectConfig()->phases[0],
                                                          obj->getParticleObjectConfig()->phases[1]});
            alpha_all.insert(alpha_all.end(), alpha_tmp.begin(), alpha_tmp.end());

            m_host_const.total_particle_num += pos_tmp.size();
        }
    }

    bool JL21CTSolver::initialize() {
        if (!m_config) {
            std::cout << "ERROR:: solver config empty.\n";
            return false;
        }

        if (m_objects.empty()) {
            std::cout << "ERROR:: solver attach no object.\n";
            return false;
        }

        mergeObjects();

        auto solver_config = dynamic_cast<JL21CTSolverConfig *>(m_config.get());
        int device;
        cudaGetDevice(&device);
        cudaDeviceProp prop{};
        cudaGetDeviceProperties(&prop, device);
        solver_config->kernel_threads = prop.maxThreadsPerBlock;
        solver_config->kernel_blocks = std::ceil(
                (m_host_const.total_particle_num + solver_config->kernel_threads - 1) / solver_config->kernel_threads);

        auto particle_radius = (*m_objects.begin())->getParticleObjectConfig()->particle_radius.value();
        auto particle_num = m_host_const.total_particle_num;

        // setup m_host_const
        m_host_const.dt = solver_config->dt;
        m_host_const.cur_sim_time = solver_config->cur_sim_time;
        m_host_const.gravity = solver_config->gravity;
        m_host_const.block_num = solver_config->kernel_blocks;
        m_host_const.thread_num = solver_config->kernel_threads;
        m_host_const.sph_h = 4 * particle_radius;
        m_host_const.rest_volume = powf(particle_radius * 2, 3);
        m_host_const.total_particle_num = pos_all.size();

        m_host_const.phase1_color = solver_config->phase1_color;
        m_host_const.phase2_color = solver_config->phase2_color;
        m_host_const.rest_density = solver_config->rest_density;
        m_host_const.kd = solver_config->kd;
        m_host_const.Cd = solver_config->Cd;
        m_host_const.relation_time = solver_config->relation_time;
        m_host_const.alpha_up_threshold = solver_config->alpha_up_threshold;
        m_host_const.alpha_low_threshold = solver_config->alpha_low_threshold;
        m_host_const.compliance_0 = solver_config->compliance_0;
        m_host_const.compliance_inf = solver_config->compliance_inf;

        m_host_const.shear_modulus = 10000;
        m_host_const.Y = 100;

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
        cudaMalloc((void **) &m_device_const, sizeof(JL21CTConstantParams));
        cudaMalloc((void **) &m_device_data, sizeof(JL21CTDynamicParams));
        m_host_data.malloc(particle_num);
        m_neighborSearch.malloc();

        // data copy
        cudaMemcpy(m_device_const, &m_host_const, sizeof(JL21CTConstantParams), cudaMemcpyHostToDevice);
        cudaMemcpy(m_device_data, &m_host_data, sizeof(JL21CTDynamicParams), cudaMemcpyHostToDevice);

        cudaMemcpy(m_host_data.pos, pos_all.data(), particle_num * sizeof(Vec3f), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.predictPos, pos_all.data(), particle_num * sizeof(Vec3f), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.mat, mat_all.data(), particle_num * sizeof(Material), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.alpha, alpha_all.data(), particle_num * sizeof(Vec2f), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.alpha_last, alpha_all.data(), particle_num * sizeof(Vec2f), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.vel_mix, vel_all.data(), particle_num * sizeof(Vec3f), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.vel_k1, vel_all.data(), particle_num * sizeof(Vec3f), cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_data.vel_k2, vel_all.data(), particle_num * sizeof(Vec3f), cudaMemcpyHostToDevice);

        if (cudaGetLastError() == cudaSuccess) {
            std::cout << "JL21CTSolver initialized.\n";
            m_is_init = true;
            return true;
        }
        return false;
    }

    void JL21CTSolver::deviceMalloc() {
        // const
        cudaMalloc((void **) &m_device_const, sizeof(m_host_const));

        // data
        m_host_data.malloc(m_host_const.total_particle_num);
        cudaMalloc((void **) &m_device_data, sizeof(m_host_data));

        // neighbor_search
        m_neighborSearch.malloc();
    }

    void JL21CTSolver::destroy() {
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
                std::cout << "JL21CTSolver destroyed.\n";
        }
    }

    void JL21CTSolver::exportAsPly() {
        static int counter = 1;
        std::vector<Vec3f> pos(m_host_const.total_particle_num);
        cudaMemcpy(pos.data(),
                   m_host_data.pos,
                   m_host_const.total_particle_num * sizeof(Vec3f),
                   cudaMemcpyDeviceToHost);
        ModelExporter::exportVecSetAsPly("F:\\DataSet.Research\\ITEM.NN_NewModel\\ply\\jl21ct_test",
                                         std::to_string(counter++),
                                         pos);
    }

    void JL21CTSolver::run(float total_time) {
        if (!m_is_init)
            initialize();

        std::cout << "JL21CT run.\n";

        auto solver_config = dynamic_cast<JL21CTSolverConfig *>(m_config.get());
        if (m_is_start) {

            m_neighborSearch.update(m_host_data.pos);
            InitOtherData(m_host_const, m_device_const, m_device_data);

            m_is_start = false;
        }

        while (solver_config->cur_sim_time <= total_time) {

            if (dynamic_cast<JL21CTSolverConfig *>(m_config.get())->export_data)
                exportAsPly();

            step();


            std::cout << "Frame out: " << unsigned(solver_config->cur_sim_time / solver_config->dt) << '\n';
        }
    }

    void JL21CTSolver::step() {
        auto solver_config = dynamic_cast<JL21CTSolverConfig *>(m_config.get());

        auto d_nsConst = m_neighborSearch.d_config;
        auto d_nsParams = m_neighborSearch.d_params;

        EstimateDensity(m_host_const, m_device_const, m_device_data, d_nsConst, d_nsParams);

        ComputeVelGrad(m_host_const, m_device_const, m_device_data, d_nsConst, d_nsParams);

        ComputeCompliance(m_host_const, m_device_const, m_device_data, d_nsParams);

        ComputeConformationTensorForce(m_host_const, m_device_const, m_device_data, d_nsConst, d_nsParams);

        ComputeFluidVisForce(m_host_const, m_device_const, m_device_data, d_nsConst, d_nsParams);

        ComputeCohesionForce(m_host_const, m_device_const, m_device_data, d_nsConst, d_nsParams);

        ComputeFluidPressureForce(m_host_const, m_device_const, m_device_data, d_nsConst, d_nsParams);

        ComputePhaseVel(m_host_const, m_device_const, m_device_data, d_nsParams);

        ComputeMixVel(m_host_const, m_device_const, m_device_data, d_nsParams);

        ComputeDragForce(m_host_const, m_device_const, m_device_data, d_nsParams);

        ComputeMixVel(m_host_const, m_device_const, m_device_data, d_nsParams);

        Advect(m_host_const, m_device_const, m_device_data, d_nsParams);

        Reset(m_host_const, m_device_const, m_device_data, d_nsParams);

        m_neighborSearch.update(m_host_data.pos);

        PhaseTransport(m_host_const, m_device_const, m_device_data, d_nsConst, d_nsParams);

        UpdateColor(m_host_const, m_device_const, m_device_data, d_nsParams);

        solver_config->cur_sim_time += solver_config->dt;
    }

}