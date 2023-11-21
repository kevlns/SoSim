//@author        : Long Shen
//@date          : 2023/11/8
//@description   :
//@version       : 1.0

#include <fstream>

#include "Public/PhysicalSolvers/MixtureModelForNonNewtonFlow/mmsph_solver.hpp"
#include "Private/PhysicalSolvers/MixtureModelForNonNewtonFlow/cuda_api.cuh"
#include "Public/Shared/NeighborSearchUGB/neighbor_search_config.hpp"
#include "Public/ThirdParty/json/json.hpp"
#include "Public/Shared/ModelUtils/model_tool.hpp"
#include "Public/Shared/Math/helper_math.hpp"
#include "Public/Shared/Analyzer/params_analysis_helper.hpp"

namespace SoSim::MMSPH {

    void MMSPHSolver::initialize() {
        if (!m_config) {
            std::cout << "ERROR::MMSPHSolver::initialize(): not set solverConfig yet.\n";
            return;
        }

        if (m_host_pos.empty()) {
            std::cout << "ERROR::MMSPHSolver::initialize(): have no particles.\n";
            return;
        }

        m_isInit = true;
        m_total_particle_num = m_host_pos.size();

        // reset some param in m_config
        m_config->block_num = (m_total_particle_num + m_config->thread_num - 1) / m_config->thread_num;

        // setup neighbor search
        m_neighborSearcher = new NeighborSearcher;
        SoSim::NSUGB::NeighborSearchConfig ns_config;
        ns_config.block_num = m_config->block_num;
        ns_config.thread_num = m_config->thread_num;
        ns_config.scene_lb = m_config->scene_lb;
        ns_config.scene_size = m_config->scene_size;
        ns_config.max_neighbor_num = 35;
        ns_config.total_particle_num = m_total_particle_num;
        ns_config.particle_radius = m_config->unified_particle_radius;
        m_neighborSearcher->setConfig(&ns_config);

        // setup device resource
        m_host_cp.dt = m_config->dt;
        m_host_cp.rest_volume = 8 * powf(m_config->unified_particle_radius, 3);
        m_host_cp.sph_h = m_config->unified_particle_radius * 4;
        m_host_cp.total_particle_num = m_total_particle_num;
        m_host_cp.max_neighbor_num = ns_config.max_neighbor_num;
        m_host_cp.gravity = m_config->gravity;

        cudaMalloc_t((void **) &m_device_cp, sizeof(ConstParams), m_mem);
        cudaMemcpy(m_device_cp, &m_host_cp, sizeof(ConstParams), cudaMemcpyHostToDevice);
        m_isInit &= cudaGetLastError_t("MMSPHSolver::initialize(): init m_device_cp failed.");

        auto partNum = ns_config.total_particle_num;
        size_t size1 = partNum * sizeof(float3);
        size_t size2 = partNum * sizeof(float2);
        size_t size3 = partNum * sizeof(float);
        size_t size4 = partNum * sizeof(Material);
        size_t size5 = partNum * sizeof(uint8_t);

        cudaMalloc_t((void **) &m_host_dp.pos, size1, m_mem);
        cudaMalloc_t((void **) &m_host_dp.predictPos, size1, m_mem);
        cudaMalloc_t((void **) &m_host_dp.acc_mix, size1, m_mem);
        cudaMalloc_t((void **) &m_host_dp.vel_mix, size1, m_mem);
        cudaMalloc_t((void **) &m_host_dp.drift_vel_k1, size1, m_mem);
        cudaMalloc_t((void **) &m_host_dp.drift_vel_k2, size1, m_mem);
        cudaMalloc_t((void **) &m_host_dp.M_m, size1, m_mem);
        cudaMalloc_t((void **) &m_host_dp.alpha, size2, m_mem);
        cudaMalloc_t((void **) &m_host_dp.density_mix, size3, m_mem);
        cudaMalloc_t((void **) &m_host_dp.density_sph, size3, m_mem);
        cudaMalloc_t((void **) &m_host_dp.pressure_mix, size3, m_mem);
        cudaMalloc_t((void **) &m_host_dp.bPart_volume, size3, m_mem);
        cudaMalloc_t((void **) &m_host_dp.mat, size4, m_mem);
        cudaMalloc_t((void **) &m_host_dp.isDynamic, size5, m_mem);

        cudaMemcpy(m_host_dp.pos, m_host_pos.data(), size1, cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_dp.predictPos, m_host_pos.data(), size1, cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_dp.vel_mix, m_host_vel.data(), size1, cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_dp.alpha, m_host_phase.data(), size2, cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_dp.mat, m_host_mat.data(), size4, cudaMemcpyHostToDevice);
        cudaMemcpy(m_host_dp.isDynamic, m_host_isDynamic.data(), size5, cudaMemcpyHostToDevice);
        m_isInit &= cudaGetLastError_t("MMSPHSolver::initialize(): init m_host_dp failed.");

        cudaMalloc_t((void **) &m_device_dp, sizeof(DynamicParams), m_mem);
        cudaMemcpy(m_device_dp, &m_host_dp, sizeof(DynamicParams), cudaMemcpyHostToDevice);
        m_isInit &= cudaGetLastError_t("MMSPHSolver::initialize(): init m_device_dp failed.");

    }

    bool MMSPHSolver::isInitialized() const {
        return m_isInit;
    }

    void MMSPHSolver::run() {
        if (m_isStart) {
            initialize();
            init_data(m_device_cp,
                      m_device_dp,
                      m_config->block_num,
                      m_config->thread_num);
            m_neighborSearcher->update(m_host_dp.pos);
            precompute_bPart_volume(m_device_cp,
                                    m_device_dp,
                                    m_neighborSearcher->getPartIndexDevicePtr(),
                                    m_neighborSearcher->getNeighborsDevicePtr(),
                                    m_config->block_num,
                                    m_config->thread_num);
            m_isStart = false;
        }

        if (!m_isInit) {
            std::cout << "ERROR::MMSPHSolver::run(): initialize failed.\n";
            return;
        }

        std::string restore_path_prefix = "C:\\Users\\sl936\\Desktop\\experiment\\ply\\sosim-test\\";
        for (int i = 1; i < 2001; ++i) {
            step();

            auto restore_path = restore_path_prefix + std::to_string(i) + ".ply";
            write_ply(restore_path, m_host_dp.pos, 0, m_total_particle_num);
            std::cout << "Frame out: " << i << '\n';
        }

    }

    void MMSPHSolver::destroy() {
        bool isOk = true;
        if (m_isInit) {

            m_neighborSearcher->destroy();
            m_neighborSearcher = nullptr;

            m_host_dp.destroy();

            cudaFree(m_device_cp);
            cudaFree(m_device_dp);
            m_device_cp = nullptr;
            m_device_dp = nullptr;

            delete m_config;
            m_config = nullptr;

            isOk = cudaGetLastError_t("MMSPHSolver::destroy(): cudaFree m_device_*p failed.");
        }
        if (isOk)
            std::cout << "MMSPHSolver destructed.\n";
    }

    void MMSPHSolver::setSolverConfig(SolverConfig *solverConfig, const SceneConfig *sceneConfig) {
        solverConfig->block_num = sceneConfig->block_num;
        solverConfig->thread_num = sceneConfig->thread_num;
        solverConfig->scene_lb = sceneConfig->scene_lb;
        solverConfig->scene_size = sceneConfig->scene_size;

        if (!m_config) {
            m_config = new SolverConfig;
            memcpy(m_config, solverConfig, sizeof(SolverConfig));
        } else {
            std::cout << "ERROR:: MMSPHSolverConfig already setup.\n";
        }
    }

    void MMSPHSolver::attachObject(const Object *obj) {
        // TODO
    }

    void MMSPHSolver::addParticles(const std::string &obj_json_path) {
        using json = nlohmann::json;

        std::ifstream f(obj_json_path);
        json config = json::parse(f);

        float r = config["unified_particle_radius"];
        m_config->unified_particle_radius = r;

        for (auto obj: config["objs"]) {
            if (!obj["model_path"].empty())
                // TODO load model
                return;
            else {

                // get default config
                std::vector<float3> newPos = gen_pos(obj["default"], r);

                // get obj phase
                float2 phase{0, 0};
                if (obj["phase"] == 1) {
                    phase.x = 1;
                    m_host_cp.rest_density.x = obj["rest_density"];
                } else if (obj["phase"] == 2) {
                    phase.y = 1;
                    m_host_cp.rest_density.y = obj["rest_density"];
                }
                std::vector<float2> newPhase(newPos.size(), phase);

                // get obj if is dynamic
                uint8_t isDynamic = obj["is_dynamic"];
                std::vector<uint8_t> newIsDynamic(newPos.size(), isDynamic);

                // get obj material
                Material mat;
                if (obj["mat"] == "fluid") {
                    mat = Material::FLUID;
                    m_fluid_particle_num += newPos.size();
                } else if (obj["mat"] == "bound") {
                    mat = Material::BOUND;
                    m_bound_particle_num += newPos.size();
                }
                std::vector<Material> newMat(newPos.size(), mat);

                // get vel_start
                float3 v_start = make_float3(obj["vel_start"].get<std::vector<float>>());
                std::vector<float3> newVelStart(newPos.size(), v_start);

                m_host_pos.insert(m_host_pos.end(), newPos.begin(), newPos.end());
                m_host_vel.insert(m_host_vel.end(), newVelStart.begin(), newVelStart.end());
                m_host_phase.insert(m_host_phase.end(), newPhase.begin(), newPhase.end());
                m_host_mat.insert(m_host_mat.end(), newMat.begin(), newMat.end());
                m_host_isDynamic.insert(m_host_isDynamic.end(), newIsDynamic.begin(), newIsDynamic.end());
            }
        }

    }

    void MMSPHSolver::step() {

        compute_ext_force(m_device_cp,
                          m_device_dp,
                          m_neighborSearcher->getPartIndexDevicePtr(),
                          m_neighborSearcher->getNeighborsDevicePtr(),
                          m_config->block_num,
                          m_config->thread_num);

//        dump_max(m_fluid_particle_num, m_host_dp.acc_mix);

        advect_pos(m_device_cp,
                   m_device_dp,
                   m_neighborSearcher->getPartIndexDevicePtr(),
                   m_neighborSearcher->getNeighborsDevicePtr(),
                   m_config->block_num,
                   m_config->thread_num);

        m_neighborSearcher->update(m_host_dp.pos);

        update_density(m_device_cp,
                       m_device_dp,
                       m_neighborSearcher->getPartIndexDevicePtr(),
                       m_neighborSearcher->getNeighborsDevicePtr(),
                       m_config->block_num,
                       m_config->thread_num);

        dump_max(m_fluid_particle_num, m_host_dp.density_sph);

        compute_WC_pressure_force(m_device_cp,
                                  m_device_dp,
                                  m_neighborSearcher->getPartIndexDevicePtr(),
                                  m_neighborSearcher->getNeighborsDevicePtr(),
                                  m_config->block_num,
                                  m_config->thread_num);

        compute_viscous_force(m_device_cp,
                              m_device_dp,
                              m_neighborSearcher->getPartIndexDevicePtr(),
                              m_neighborSearcher->getNeighborsDevicePtr(),
                              m_config->block_num,
                              m_config->thread_num);

        m_config->cur_sim_time += m_config->dt;
    }

}