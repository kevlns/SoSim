//
// Created by ADMIN on 2024/3/14.
//
#include <memory>
#include <cuda_runtime.h>

#include "solvers/WCSPH/wcsph_solver.hpp"
#include "wcsph_cuda_api.cuh"
#include "libs/ModelL/model_exporter.hpp"
#include "libs/AnalysisL/statistic_util.hpp"

namespace SoSim {

    WCSPHSolver::WCSPHSolver() {
        m_config = std::make_shared<WCSPHSolverConfig>();
        std::cout << "Create WCSPHSolver.\n";
    }

    WCSPHSolver::~WCSPHSolver() {
        destroy();
    }

    std::shared_ptr<SolverConfig> WCSPHSolver::getConfig() {
        return m_config;
    }

    void WCSPHSolver::attachObject(std::shared_ptr<Object> object) {
        if (m_objects.count(object) == 0)
            m_objects.insert(object);

        m_change_occur = true;
        std::cout << "WCSPHSolver attach object: " << object->getName() << ".\n";
    }

    void WCSPHSolver::detachObject(std::shared_ptr<Object> object) {
        if (m_objects.count(object) > 0)
            m_objects.erase(object);

        m_change_occur = true;
        std::cout << "WCSPHSolver detach object: " << object->getName() << ".\n";
    }

    void WCSPHSolver::mergeObjects() {
        m_host_const.particle_num = 0;
        for (const auto &obj: m_objects)
            m_host_const.particle_num += obj->getParticleNum();
    }

    bool WCSPHSolver::initialize() {
        if (!m_config) {
            std::cout << "ERROR:: solver config empty.\n";
            return false;
        }

        if (m_objects.empty()) {
            std::cout << "ERROR:: solver attach no object.\n";
            return false;
        }

        mergeObjects();

        auto solver_config = dynamic_cast<WCSPHSolverConfig *>(m_config.get());
        auto particle_radius = (*m_objects.begin())->getParticleObjectConfig()->particle_radius;

        int device;
        cudaGetDevice(&device);
        cudaDeviceProp prop{};
        cudaGetDeviceProperties(&prop, device);
        solver_config->kernel_threads = prop.maxThreadsPerBlock;
        solver_config->kernel_blocks = std::ceil(
                (m_host_const.particle_num + solver_config->kernel_threads - 1) / solver_config->kernel_threads);

        m_host_const.dt = solver_config->dt;
        m_host_const.h = 4 * particle_radius;
        m_host_const.rest_vis = solver_config->rest_vis;
        m_host_const.rest_density = solver_config->rest_density;
        m_host_const.rest_volume = std::pow(2 * particle_radius, 3);
        m_host_const.rest_rigid_density = solver_config->rest_rigid_density;
        m_host_const.stiff = m_host_const.rest_density * solver_config->cs * solver_config->cs / 7;
        m_host_const.stiff = 1000000;
//        m_host_const.stiff = 1000000;
        m_host_const.gravity = solver_config->gravity;

        NeighborSearchUGConfig ns_config;
        ns_config.sceneLB = solver_config->scene_lb;
        ns_config.sceneSize = solver_config->scene_size;
        ns_config.cellLength = m_host_const.h;
        ns_config.gridSize = {
                static_cast<uint32_t>(std::ceil(ns_config.sceneSize.x / ns_config.cellLength)),
                static_cast<uint32_t>(std::ceil(ns_config.sceneSize.y / ns_config.cellLength)),
                static_cast<uint32_t>(std::ceil(ns_config.sceneSize.z / ns_config.cellLength))};
        ns_config.cellNum = ns_config.gridSize.x * ns_config.gridSize.y * ns_config.gridSize.z;
        ns_config.particle_num = m_host_const.particle_num;
        ns_config.maxNeighborNum = solver_config->max_neighborNum;
        ns_config.kernel_threads = solver_config->kernel_threads;
        ns_config.kernel_blocks = solver_config->kernel_blocks;
        m_neighborSearch.setConfig(ns_config);

        deviceMalloc();

        // data cpy
        unsigned offset = 0;
        for (const auto &obj: m_objects) {
            std::vector<Vec3f> vel_start(obj->getParticleNum(), obj->getParticleObjectConfig()->vel_start);
            std::vector<Material> mat(obj->getParticleNum(), obj->getParticleObjectConfig()->particle_mat.value());
            cudaMemcpy(m_host_data.m_device_pos + offset,
                       obj->getParticles().data(),
                       obj->getParticleNum() * sizeof(Vec3f),
                       cudaMemcpyHostToDevice);
            cudaMemcpy(m_host_data.m_device_vel + offset,
                       vel_start.data(),
                       obj->getParticleNum() * sizeof(Vec3f),
                       cudaMemcpyHostToDevice);
            cudaMemcpy(m_host_data.m_device_mat + offset,
                       mat.data(),
                       obj->getParticleNum() * sizeof(Material),
                       cudaMemcpyHostToDevice);
            offset += obj->getParticleNum();
        }

        cudaMemcpy(m_device_const, &m_host_const, sizeof(m_host_const), cudaMemcpyHostToDevice);
        cudaMemcpy(m_device_data, &m_host_data, sizeof(m_host_data), cudaMemcpyHostToDevice);

        if (cudaGetLastError() == cudaSuccess) {
            std::cout << "WCSPHSolver initialized.\n";
            m_is_init = true;
            return true;
        }
        return false;
    }

    void WCSPHSolver::deviceMalloc() {
        // const
        cudaMalloc((void **) &m_device_const, sizeof(m_host_const));

        // data
        m_host_data.malloc(m_host_const.particle_num);
        cudaMalloc((void **) &m_device_data, sizeof(m_host_data));

        // neighbor_search
        m_neighborSearch.malloc();
    }

    void WCSPHSolver::destroy() {
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
                std::cout << "WCSPHSolver destroyed.\n";
        }
    }

    void WCSPHSolver::exportAsPly() {
        static int counter = 1;
        std::vector<Vec3f> pos(m_host_const.particle_num);
        cudaMemcpy(pos.data(),
                   m_host_data.m_device_pos,
                   m_host_const.particle_num * sizeof(Vec3f),
                   cudaMemcpyDeviceToHost);
        ModelExporter::exportVecSetAsPly("F:\\DataSet.Research\\ITEM.NN_NewModel\\ply\\wcsph_test",
                                         std::to_string(counter++),
                                         pos);
    }

    void WCSPHSolver::run(float total_time) {
        if (!m_is_init)
            initialize();

        std::cout << "WCSPH run.\n";

        auto solver_config = dynamic_cast<WCSPHSolverConfig *>(m_config.get());
        if (m_is_start) {

            m_neighborSearch.update(m_host_data.m_device_pos);
            init(m_neighborSearch.h_config,
                 m_device_const,
                 m_device_data,
                 m_neighborSearch.d_config,
                 m_neighborSearch.d_params);

            m_is_start = false;
        }

        while (solver_config->cur_sim_time <= total_time) {

            if (dynamic_cast<WCSPHSolverConfig *>(m_config.get())->export_data)
                exportAsPly();

            step();


            std::cout << "Frame out: " << unsigned(solver_config->cur_sim_time / solver_config->dt) << '\n';
        }
    }

    void WCSPHSolver::step() {
        auto solver_config = dynamic_cast<WCSPHSolverConfig *>(m_config.get());

        computeGravityForce(
                m_neighborSearch.h_config,
                m_device_const,
                m_device_data,
                m_neighborSearch.d_config,
                m_neighborSearch.d_params
        );

        advect(
                m_neighborSearch.h_config,
                m_device_const,
                m_device_data,
                m_neighborSearch.d_config,
                m_neighborSearch.d_params
        );

        m_neighborSearch.update(m_host_data.m_device_pos);

        computeDensityAndPressure(
                m_neighborSearch.h_config,
                m_device_const,
                m_device_data,
                m_neighborSearch.d_config,
                m_neighborSearch.d_params
        );

        dump_mean(m_host_data.m_device_pressure,
                  m_host_const.particle_num,
                  m_host_const.particle_num);

        computePressureForce(
                m_neighborSearch.h_config,
                m_device_const,
                m_device_data,
                m_neighborSearch.d_config,
                m_neighborSearch.d_params
        );

        computeViscousForce(
                m_neighborSearch.h_config,
                m_device_const,
                m_device_data,
                m_neighborSearch.d_config,
                m_neighborSearch.d_params
        );

        solver_config->cur_sim_time += solver_config->dt;
    }
}

