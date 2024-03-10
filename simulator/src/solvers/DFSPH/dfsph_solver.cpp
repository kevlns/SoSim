//
// Created by ADMIN on 2024/3/7.
//

#include "solvers/DFSPH/dfsph_solver.hpp"

namespace SoSim {


    void DFSPHSolver::setConfig(DFSPHSolverConfig *config) {
        m_config = config;
    }

    DFSPHSolverConfig *DFSPHSolver::getConfig() {
        return m_config.value();
    }

    void DFSPHSolver::attachObject(Object *object) {
        if (m_objects.count(object) == 0)
            m_objects.insert(object);

        m_changeOccur = true;
        std::cout << "DFSPHSolver attach object: " << object->getName() << ".\n";
    }

    void DFSPHSolver::detachObject(Object *object) {
        if (m_objects.count(object) > 0)
            m_objects.erase(object);

        m_changeOccur = true;
        std::cout << "DFSPHSolver detach object: " << object->getName() << ".\n";
    }

    void DFSPHSolver::mergeObjects() {
        m_host_const.particle_num = 0;
        for (auto obj: m_objects)
            m_host_const.particle_num += obj->getParticleNum();
    }

    bool DFSPHSolver::initialize() {
        if (!m_config.has_value()) {
            std::cout << "ERROR:: solver config empty.\n";
            return false;
        }

        if (m_objects.empty()) {
            std::cout << "ERROR:: solver attach no object.\n";
            return false;
        }

        auto solver_config = reinterpret_cast
                <DFSPHSolverConfig *>(m_config.value());
        auto particle_radius = (*m_objects.begin())->getParticleObjectConfig()->particle_radius;

        int device;
        cudaGetDevice(&device);
        cudaDeviceProp prop{};
        cudaGetDeviceProperties(&prop, device);
        solver_config->kernel_threads = prop.maxThreadsPerBlock;
        solver_config->kernel_blocks = std::ceil(
                (m_host_const.particle_num + solver_config->kernel_threads - 1) / solver_config->kernel_threads);

        mergeObjects();
        m_host_const.dt = solver_config->dt;
        m_host_const.h = 4 * particle_radius;
        m_host_const.rest_vis = solver_config->rest_vis;
        m_host_const.rest_density = solver_config->rest_density;
        m_host_const.rest_volume = std::pow(particle_radius, 3);
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
        for (auto obj: m_objects) {
            std::vector<Vec3f> vel_start(obj->getParticleNum(), obj->getParticleObjectConfig()->vel_start);
            std::vector<Material> mat(obj->getParticleNum(), obj->getParticleObjectConfig()->particle_mat.value());
            cudaMemcpy(m_host_data.pos + offset * sizeof(Vec3f),
                       obj->getParticles().data(),
                       obj->getParticleNum() * sizeof(Vec3f),
                       cudaMemcpyHostToDevice);
            cudaMemcpy(m_host_data.vel + offset * sizeof(Vec3f),
                       vel_start.data(),
                       obj->getParticleNum() * sizeof(Vec3f),
                       cudaMemcpyHostToDevice);
            cudaMemcpy(m_host_data.mat + offset * sizeof(Material),
                       mat.data(),
                       obj->getParticleNum() * sizeof(Material),
                       cudaMemcpyHostToDevice);
            offset += obj->getParticleNum();
        }

        cudaMemcpy(m_device_const, &m_host_const, sizeof(m_host_const), cudaMemcpyHostToDevice);
        cudaMemcpy(m_device_data, &m_host_data, sizeof(m_host_data), cudaMemcpyHostToDevice);

        if (cudaGetLastError() == cudaSuccess) {
            std::cout << "DFSPHSolver initialized.\n";
            m_isInit = true;
            return true;
        }
        return false;
    }

    void DFSPHSolver::deviceMalloc() {
        // const
        cudaMalloc((void **) &m_device_const, sizeof(m_host_const));

        // data
        m_host_data.malloc(m_host_const.particle_num);
        cudaMalloc((void **) &m_device_data, sizeof(m_host_data));

        // neighbor_search
        m_neighborSearch.malloc();
    }

    void DFSPHSolver::destroy() {
        m_objects.clear();

        if (m_config.has_value())
            delete m_config.value();

        if (m_isInit) {
            // delete neighbor_search
            cudaFree(m_device_const);

            // delete data
            m_host_data.freeMemory();
            cudaFree(m_device_data);

            // delete const
            m_neighborSearch.freeMemory();

            if (cudaGetLastError() == cudaSuccess)
                std::cout << "DFSPHSolver destroyed.\n";

        }

    }

    void DFSPHSolver::run(float total_time) {
        m_neighborSearch.update(m_host_data.pos);
        m_neighborSearch.dump();
    }

    void DFSPHSolver::step() {

    }
}
