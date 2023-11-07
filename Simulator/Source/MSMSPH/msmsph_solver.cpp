//@author        : Long Shen
//@date          : 2023/10/24
//@description   :
//@version       : 1.0

#include <fstream>
#include <string>

#include "Public/MSMSPH/msmsph_solver.hpp"
#include "Private/MSMSPH/cuda_api.cuh"
#include "Public/Shared/Analyzer/params_analysis_helper.hpp"
#include "Public/Shared/ModelUtils/model_tool.hpp"

namespace SoSim::MSMSPH {

    void MSMSPHSolver::initialize() {
        if (!(m_rIsInit & m_sIsInit)) {
            std::cout << "ERROR:: Pre-Setup is incomplete.\n";
            return;
        }

        m_isInit = true;

        initCP();

        initDP();

        int device;
        cudaGetDevice(&device);

        cudaDeviceProp prop{};
        cudaGetDeviceProperties(&prop, device);

        m_threadNum = prop.maxThreadsPerBlock;
        m_blockNum = (m_host_cp.total_particle_num + m_threadNum - 1) / m_threadNum;
    }

    void MSMSPHSolver::initCP() {
        m_host_cp.gravity = m_solverConfig.gravity;
        m_host_cp.dt = static_cast<float>(m_solverConfig.dt);
        m_host_cp.rest_volume = static_cast<float>(pow(m_particleRadius, 3));
        m_host_cp.ns_maxNeighborNum = 35;
        m_host_cp.sph_h = 4 * m_particleRadius;

        m_neighborSearcher = new NeighborSearcher;
        m_neighborSearcher->initialize(m_sceneLB, m_sceneSize, m_host_cp.total_particle_num, m_host_cp.sph_h);

        cudaMalloc_t((void **) &m_device_cp, sizeof(ConstParams), m_mem);
        cudaMemcpy(m_device_cp, &m_host_cp, sizeof(ConstParams), cudaMemcpyHostToDevice);

        m_isInit &= cudaGetLastError_t("ERROR::MSMSPHSolver::initCP() failed.");
    }

    void MSMSPHSolver::initDP() {
        size_t size1 = m_host_cp.total_particle_num * sizeof(float3);
        size_t size2 = m_host_cp.total_particle_num * sizeof(float2);
        size_t size3 = m_host_cp.total_particle_num * sizeof(float);
        size_t size4 = m_host_cp.total_particle_num * sizeof(Material);
        size_t size5 = m_host_cp.total_particle_num * sizeof(Phase);

        cudaMalloc_t((void **) &m_host_dp.pos, size1, m_mem);
        cudaMalloc_t((void **) &m_host_dp.predictPos, size1, m_mem);
        cudaMalloc_t((void **) &m_host_dp.v_m, size1, m_mem);
        cudaMalloc_t((void **) &m_host_dp.v_mk1, size1, m_mem);
        cudaMalloc_t((void **) &m_host_dp.v_mk2, size1, m_mem);
        cudaMalloc_t((void **) &m_host_dp.acc, size1, m_mem);
        cudaMalloc_t((void **) &m_host_dp.volF_k, size2, m_mem);
        cudaMalloc_t((void **) &m_host_dp.d_volF, size2, m_mem);
        cudaMalloc_t((void **) &m_host_dp.pressure_k, size2, m_mem);
        cudaMalloc_t((void **) &m_host_dp.pressure, size3, m_mem);
        cudaMalloc_t((void **) &m_host_dp.mass, size3, m_mem);
        cudaMalloc_t((void **) &m_host_dp.density, size3, m_mem);
        cudaMalloc_t((void **) &m_host_dp.density_ba, size3, m_mem);
        cudaMalloc_t((void **) &m_host_dp.mat, size4, m_mem);
        cudaMalloc_t((void **) &m_host_dp.original_phase, size5, m_mem);

        cudaMemcpy(m_host_dp.pos, m_host_pos.data(), size1, cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_dp.predictPos, m_host_pos.data(), size1, cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_dp.v_m, m_host_vel.data(), size1, cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_dp.density, m_host_den.data(), size3, cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_dp.mat, m_host_mat.data(), size4, cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_dp.original_phase, m_host_oPhase.data(), size5, cudaMemcpyHostToDevice);

        cudaMalloc_t((void **) &m_device_dp, sizeof(DynamicParams), m_mem);
        cudaMemcpy(m_device_dp, &m_host_dp, sizeof(DynamicParams), cudaMemcpyHostToDevice);

        m_isInit &= cudaGetLastError_t("ERROR::MSMSPHSolver::initDP() failed.");
    }

    bool MSMSPHSolver::isInitialized() const {
        return m_isInit;
    }

    void MSMSPHSolver::run() {
        for (int i = 1; i < 101; ++i) {
            step();

            size_t size = m_host_cp.total_particle_num * sizeof(float3);
            cudaMemcpy(m_host_pos.data(), m_host_dp.pos, size, cudaMemcpyDeviceToHost);
            write_ply("C:\\Users\\sl936\\Desktop\\experiment\\ply\\sosim-test\\" + std::to_string(i) + ".ply",
                      m_host_pos);
            std::cout << "Frame out: " << i << "\n";

        }
    }

    void MSMSPHSolver::destroy() {

        this->m_host_dp.destroy();

        if (m_neighborSearcher) {
            m_neighborSearcher->destroy();
            delete m_neighborSearcher;
            m_neighborSearcher = nullptr;
        }

        if (m_bNeighborSearcher) {
            m_bNeighborSearcher->destroy();
            delete m_bNeighborSearcher;
            m_bNeighborSearcher = nullptr;
        }

        m_device_cp = nullptr;
        m_device_dp = nullptr;
    }

    void MSMSPHSolver::setConfig(const SolverConfig &config) {
        m_solverConfig = config;
    }

    void MSMSPHSolver::step() {

        if (m_isStart) {
            init_data(m_device_cp, m_device_dp, m_blockNum, m_threadNum);
            m_neighborSearcher->update(m_host_dp.pos);
            update_density_and_pressure(m_device_cp, m_device_dp, m_neighborSearcher->getPartIndexDevicePtr(),
                                        m_neighborSearcher->getNeighborsDevicePtr(), m_blockNum, m_threadNum);
            m_isStart = false;
        }

//        compute_drift_vel(m_device_cp, m_device_dp, m_neighborSearcher->getPartIndexDevicePtr(),
//                          m_neighborSearcher->getNeighborsDevicePtr(), m_blockNum, m_threadNum);

//        advect_volFrac(m_device_cp, m_device_dp, m_neighborSearcher->getPartIndexDevicePtr(),
//                       m_neighborSearcher->getNeighborsDevicePtr(), m_blockNum, m_threadNum);

        update_density_and_pressure(m_device_cp, m_device_dp, m_neighborSearcher->getPartIndexDevicePtr(),
                                    m_neighborSearcher->getNeighborsDevicePtr(), m_blockNum, m_threadNum);

        compute_overall_acc(m_device_cp, m_device_dp, m_neighborSearcher->getPartIndexDevicePtr(),
                            m_neighborSearcher->getNeighborsDevicePtr(), m_blockNum, m_threadNum);

        dump_max(m_host_cp.total_particle_num, m_host_dp.acc);

        advect_particles(m_device_cp, m_device_dp, m_neighborSearcher->getPartIndexDevicePtr(), m_blockNum,
                         m_threadNum);

        m_neighborSearcher->update(m_host_dp.predictPos);
    }

    void MSMSPHSolver::setPhaseDensity(float den1, float den2) {
        m_host_cp.rest_density = {den1, den2};
    }

    void MSMSPHSolver::setParticleRadius(float particle_radius) {
        m_particleRadius = particle_radius;
        m_host_cp.sph_h = 4 * particle_radius;
        m_rIsInit = true;
    }

    void MSMSPHSolver::setSceneInfo(float3 scene_lb, float3 scene_size) {
        m_sceneLB = scene_lb;
        m_sceneSize = scene_size;
        m_sIsInit = true;
    }

    void MSMSPHSolver::setBNeighborSearcherPtr(NeighborSearcher *bNS) {
        m_bNeighborSearcher = bNS;
    }

    NeighborSearcher *MSMSPHSolver::getNeighborSearcherPtr() const {
        return m_neighborSearcher;
    }

    void MSMSPHSolver::addObject() {}

//    void MSMSPHSolver::addParticles(const std::string &obj_json) {
//        genObjFromJson(obj_json, m_host_pos, m_host_vel, m_host_den, m_host_mat, m_particleRadius);
//        m_host_cp.total_particle_num = m_host_pos.size();
//        m_rIsInit = true;
//    }

    void MSMSPHSolver::addParts(const std::string &obj_json) {
        genObjFromJson(obj_json, m_host_pos, m_host_vel, m_host_den, m_host_mat, m_host_oPhase, m_particleRadius);
        m_host_cp.total_particle_num = m_host_pos.size();
        m_rIsInit = true;
        setPhaseDensity(1000, 1000);
        setSceneInfo({-20, -20, -20}, {40, 40, 40});
    }

}
