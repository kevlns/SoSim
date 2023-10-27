//@author        : Long Shen
//@date          : 2023/10/24
//@description   :
//@version       : 1.0

#include "Public/MSMSPH/msmsph_solver.hpp"
#include "Private/MSMSPH/cuda_api.cuh"

namespace SoSim::MSMSPH {

    MSMSPHSolver::~MSMSPHSolver() {}

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

        cudaMalloc_t((void **) &m_device_dp, sizeof(DynamicParams), m_mem);
        cudaMemcpy(m_device_dp, &m_host_dp, sizeof(DynamicParams), cudaMemcpyHostToDevice);

        m_isInit &= cudaGetLastError_t("ERROR::MSMSPHSolver::initDP() failed.");
    }

    bool MSMSPHSolver::isInitialized() const {
        return m_isInit;
    }

    void MSMSPHSolver::run() {

    }

    void MSMSPHSolver::destroy() {
        m_neighborSearcher->destroy();

        m_host_dp.destroy();

        cudaFree(m_device_cp);
        cudaFree(m_device_dp);

        cudaGetLastError_t("MSMSPHSolver::destroy() failed.");

        m_neighborSearcher = nullptr;
        m_bNeighborSearcher = nullptr;
        m_device_cp = nullptr;
        m_device_dp = nullptr;
    }

    void MSMSPHSolver::setConfig(const Solver::SolverConfig &config) {
        m_solverConfig = config;
    }

    void MSMSPHSolver::step() {
        step_cuda(m_device_cp, m_device_dp, m_neighborSearcher->getPartIndexDevicePtr(),
                  m_neighborSearcher->getNeighborsDevicePtr(), nullptr, nullptr, m_blockNum, m_threadNum);
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


}
