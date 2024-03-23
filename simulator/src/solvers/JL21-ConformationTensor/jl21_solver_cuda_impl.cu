//@author        : Long Shen
//@date          : 2024/1/14
//@description   :
//@version       : 1.0

#include "jl21_solver_cuda_api.cuh"
#include "libs/SPHKernelL/kernels.cuh"

namespace SoSim {

    extern __global__ void
    checkParticleActivity(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst) {
        unsigned i = blockDim.x * blockIdx.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        d_data->active[i] = true;
        auto ur = d_nsConst->sceneLB + d_nsConst->sceneSize;
        if (d_data->pos[i].x < d_nsConst->sceneLB.x && d_data->pos[i].y < d_nsConst->sceneLB.y &&
            d_data->pos[i].z < d_nsConst->sceneLB.z && d_data->pos[i].x > ur.x && d_data->pos[i].y > ur.y &&
            d_data->pos[i].z > ur.z)
            d_data->active[i] = false;
    }

    extern __global__ void initOtherData(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data) {
        unsigned i = blockDim.x * blockIdx.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        //    d_data->vel_mix[i] *= 0;
        //    d_data->vel_k1[i] *= 0;
        //    d_data->vel_k2[i] *= 0;
        d_data->vel_drift_k1[i] *= 0;
        d_data->vel_drift_k2[i] *= 0;
        d_data->acc[i] *= 0;
        d_data->d_alpha[i] *= 0;
        d_data->lambda[i] = 1;
        d_data->density_sph[i] *= 0;
        d_data->pressure[i] *= 0;
        d_data->density_mix[i] = dot(d_const->rest_density, d_data->alpha[i]);
        if (d_data->mat[i] != JL21CT_NONNEWTON) {
            d_data->density_mix[i] = 1100;
        }
        d_data->mass[i] = d_data->density_mix[i] * d_const->rest_volume;
        d_data->color[i] = d_data->alpha[i].x * d_const->phase1_color + d_data->alpha[i].y * d_const->phase2_color;
        d_data->Q[i] = Mat33f::eye();
    }

    extern __global__ void
    computeBoundaryVolume(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
                          NeighborSearchUGParams *d_nsParams) {
        unsigned i = blockDim.x * blockIdx.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        auto p_i = d_nsParams->particleIndices_cuData[i];
        auto neib_ind = p_i * d_nsConst->maxNeighborNum;

        if (d_data->mat[p_i] == JL21CT_NONNEWTON)
            return;

        auto pos_i = d_data->predictPos[p_i];
        float delta = 0;

        for (unsigned int p_j = d_nsParams->neighbors_cuData[neib_ind], t = 0;
             p_j != UINT_MAX && t < d_nsConst->maxNeighborNum;
             ++t, p_j = d_nsParams->neighbors_cuData[neib_ind + t]) {

            auto pos_j = d_data->predictPos[p_j];
            if (d_data->mat[p_i] != JL21CT_NONNEWTON)
                delta += cubic_value((pos_i - pos_j).length(), d_const->sph_h);
        }

        d_data->volume[p_i] = 1 / delta;
    }

    extern __global__ void
    estimateDensity(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
                    NeighborSearchUGParams *d_nsParams) {
        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        auto p_i = d_nsParams->particleIndices_cuData[i];
        auto neib_ind = p_i * d_nsConst->maxNeighborNum;

        //    if (d_data->mat[p_i] != JL21CT_NONNEWTON)
        //        return;

        d_data->density_sph[p_i] *= 0.0;
        auto pos_i = d_data->predictPos[p_i];

        for (unsigned int p_j = d_nsParams->neighbors_cuData[neib_ind], t = 0;
             p_j != UINT_MAX && t < d_nsConst->maxNeighborNum;
             ++t, p_j = d_nsParams->neighbors_cuData[neib_ind + t]) {

            auto pos_j = d_data->predictPos[p_j];

            d_data->density_sph[p_i] += d_data->mass[p_j] * cubic_value((pos_i - pos_j).length(), d_const->sph_h);

            //        if (d_data->mat[p_j] == JL21CT_NONNEWTON)
            //            d_data->density_sph[p_i] += d_data->mass[p_j] * cubic_value((pos_i - pos_j).length(), d_const->sph_h);
            //        else if (d_data->mat[p_j] != JL21CT_NONNEWTON)
            //            d_data->density_sph[p_i] +=
            //                    d_data->density_mix[p_i] * d_data->volume[p_j] * cubic_value((pos_i - pos_j).length(), d_const->sph_h);
        }

        d_data->density_sph[p_i] = max(d_data->density_mix[p_i], d_data->density_sph[p_i]);
        d_data->pressure[p_i] = d_const->wc_stiffness * d_data->density_mix[p_i] *
                                (pow(d_data->density_sph[p_i] / d_data->density_mix[p_i], 7) - 1);
    }

    extern __global__ void
    computeFluidVisForce(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
                         NeighborSearchUGParams *d_nsParams) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        auto p_i = d_nsParams->particleIndices_cuData[i];
        auto neib_ind = p_i * d_nsConst->maxNeighborNum;

        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
            return;

        auto pos_i = d_data->predictPos[p_i];
        auto vel_i = d_data->vel_mix[p_i];
        Vec3f acc = {0, 0, 0};
        float vis = 0.01;
        float h2_001 = 0.001f * pow(d_const->sph_h, 2);

        for (unsigned int p_j = d_nsParams->neighbors_cuData[neib_ind], t = 0;
             p_j != UINT_MAX && t < d_nsConst->maxNeighborNum;
             ++t, p_j = d_nsParams->neighbors_cuData[neib_ind + t]) {

            auto pos_j = d_data->predictPos[p_j];
            auto vel_j = d_data->vel_mix[p_j];

            auto pi = -vis * min(0.f, dot(vel_i - vel_j, pos_i - pos_j)) /
                      ((pos_i - pos_j).length() * (pos_i - pos_j).length() + h2_001);
            acc += -d_data->mass[p_i] * d_data->mass[p_j] * pi * cubic_gradient(pos_i - pos_j, d_const->sph_h);
        }

        d_data->acc[p_i] += acc;
    }

    extern __global__ void
    computeFluidPressureForce(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
                              NeighborSearchUGParams *d_nsParams) {
        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        auto p_i = d_nsParams->particleIndices_cuData[i];
        auto neib_ind = p_i * d_nsConst->maxNeighborNum;

        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
            return;

        auto pos_i = d_data->predictPos[p_i];
        auto volume_i = d_data->mass[p_i] / d_data->density_sph[p_i];
        Vec3f acc = {0, 0, 0};

        for (unsigned int p_j = d_nsParams->neighbors_cuData[neib_ind], t = 0;
             p_j != UINT_MAX && t < d_nsConst->maxNeighborNum;
             ++t, p_j = d_nsParams->neighbors_cuData[neib_ind + t]) {

            auto pos_j = d_data->predictPos[p_j];
            auto volume_j = d_data->mass[p_j] / d_data->density_sph[p_j];

            acc -= volume_i * volume_j * (d_data->pressure[p_i] + d_data->pressure[p_j]) *
                   cubic_gradient(pos_i - pos_j, d_const->sph_h);
        }

        d_data->acc[p_i] += acc;
    }

    extern __global__ void
    computePhaseVel(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams) {
        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        auto p_i = d_nsParams->particleIndices_cuData[i];

        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
            return;

        // phase 1
        d_data->vel_k1[p_i] +=
                (d_data->acc[p_i] / (d_const->rest_density.x * d_const->rest_volume) + d_const->gravity) * d_const->dt;

        // phase 2
        d_data->vel_k2[p_i] +=
                (d_data->acc[p_i] / (d_const->rest_density.y * d_const->rest_volume) + d_const->gravity) * d_const->dt;
    }

    extern __global__ void
    computeMixVel(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        auto p_i = d_nsParams->particleIndices_cuData[i];

        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
            return;

        d_data->vel_mix[p_i] = d_data->alpha[p_i].x * d_data->vel_k1[p_i] + d_data->alpha[p_i].y * d_data->vel_k2[p_i];
        d_data->vel_drift_k1[p_i] = d_data->vel_k1[p_i] - d_data->vel_mix[p_i];
        d_data->vel_drift_k2[p_i] = d_data->vel_k2[p_i] - d_data->vel_mix[p_i];
    }

    extern __global__ void
    computeDragForce(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        auto p_i = d_nsParams->particleIndices_cuData[i];

        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
            return;

        // phase 1
        d_data->vel_k1[p_i] -= d_const->kd * d_const->dt * d_data->density_mix[p_i] / d_const->rest_density.x *
                               (d_data->vel_k1[p_i] - d_data->vel_mix[p_i]);

        // phase 2
        d_data->vel_k2[p_i] -= d_const->kd * d_const->dt * d_data->density_mix[p_i] / d_const->rest_density.y *
                               (d_data->vel_k2[p_i] - d_data->vel_mix[p_i]);
    }

    extern __global__ void
    advect(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        if (d_data->mat[i] != JL21CT_NONNEWTON)
            return;

        d_data->predictPos[i] += d_data->vel_mix[i] * d_const->dt;
        d_data->pos[i] = d_data->predictPos[i];
    }

    extern __global__ void
    reset(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        if (d_data->mat[i] != JL21CT_NONNEWTON)
            return;

        d_data->acc[i] *= 0;
        d_data->lambda[i] = 1;
        d_data->stress[i] *= 0;
    }

    extern __global__ void
    computePhaseTransport(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
                          NeighborSearchUGParams *d_nsParams) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        auto p_i = d_nsParams->particleIndices_cuData[i];
        auto neib_ind = p_i * d_nsConst->maxNeighborNum;

        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
            return;

        auto pos_i = d_data->predictPos[p_i];
        float h2_001 = 0.001f * pow(d_const->sph_h, 2);

        float d_alpha_k1 = 0;
        float d_alpha_k2 = 0;

        // phase transport
        for (unsigned int p_j = d_nsParams->neighbors_cuData[neib_ind], t = 0;
             p_j != UINT_MAX && t < d_nsConst->maxNeighborNum;
             ++t, p_j = d_nsParams->neighbors_cuData[neib_ind + t]) {

            if (d_data->mat[p_j] != JL21CT_NONNEWTON)
                continue;

            auto pos_j = d_data->predictPos[p_j];

            auto lam_ij = min(d_data->lambda[p_i], d_data->lambda[p_j]);

            // compute T_m
            auto T_m_k1 = -d_const->rest_volume * dot(d_data->alpha[p_i].x * d_data->vel_drift_k1[p_i] +
                                                      d_data->alpha[p_j].x * d_data->vel_drift_k1[p_j],
                                                      cubic_gradient(pos_i - pos_j, d_const->sph_h));
            auto T_d_k1 = -2 * d_const->Cd * d_const->rest_volume * (d_data->alpha[p_i].x - d_data->alpha[p_j].x) *
                          dot(pos_i - pos_j, cubic_gradient(pos_i - pos_j, d_const->sph_h)) /
                          ((pos_i - pos_j).length() * (pos_i - pos_j).length() + h2_001);
            d_alpha_k1 += d_const->dt * lam_ij * (T_m_k1 + T_d_k1);

            auto T_m_k2 = -d_const->rest_volume * dot(d_data->alpha[p_i].y * d_data->vel_drift_k2[p_i] +
                                                      d_data->alpha[p_j].y * d_data->vel_drift_k2[p_j],
                                                      cubic_gradient(pos_i - pos_j, d_const->sph_h));
            auto T_d_k2 = -2 * d_const->Cd * d_const->rest_volume * (d_data->alpha[p_i].y - d_data->alpha[p_j].y) *
                          dot(pos_i - pos_j, cubic_gradient(pos_i - pos_j, d_const->sph_h)) /
                          ((pos_i - pos_j).length() * (pos_i - pos_j).length() + h2_001);
            d_alpha_k2 += d_const->dt * lam_ij * (T_m_k2 + T_d_k2);
        }

        d_data->d_alpha[p_i] = {d_alpha_k1, d_alpha_k2};
    }

    extern __global__ void
    updateLambda(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams) {
        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        auto p_i = d_nsParams->particleIndices_cuData[i];

        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
            return;

        auto alpha_tmp = d_data->alpha[p_i] + d_data->d_alpha[p_i];
        float min_lam = 1.0;
        if (alpha_tmp.x < 1e-6)
            min_lam = min(min_lam, -d_data->alpha[p_i].x / d_data->d_alpha[p_i].x);
        if (alpha_tmp.y < 1e-6)
            min_lam = min(min_lam, -d_data->alpha[p_i].y / d_data->d_alpha[p_i].y);

        d_data->lambda[p_i] *= fabs(min_lam);
    }

    extern __global__ void
    phaseTransportFinalize(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams) {
        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        if (d_data->mat[i] != JL21CT_NONNEWTON)
            return;

        // update volume fraction
        d_data->alpha[i] += d_data->d_alpha[i];

        // normalize frac
        float frac_sum = 1e-6;
        if (d_data->alpha[i].x < 0)
            d_data->alpha[i].x = 0;
        else
            frac_sum += d_data->alpha[i].x;
        if (d_data->alpha[i].y < 0)
            d_data->alpha[i].y = 0;
        else
            frac_sum += d_data->alpha[i].y;

        d_data->alpha[i] /= frac_sum;

        if (d_data->alpha[i].x <= 0 && d_data->alpha[i].y <= 0)
            d_data->alpha[i] = d_data->alpha_last[i];

        // update density_mix and mass
        d_data->density_mix[i] = dot(d_data->alpha[i], d_const->rest_density);
        d_data->mass[i] = d_data->density_mix[i] * d_const->rest_volume;
        d_data->alpha_last[i] = d_data->alpha[i];
    }

    extern __global__ void
    updateColor(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams) {
        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        if (d_data->mat[i] != JL21CT_NONNEWTON)
            return;

        d_data->color[i] = d_data->alpha[i].x * d_const->phase1_color + d_data->alpha[i].y * d_const->phase2_color;
    }

    extern __global__ void
    computeVelGrad(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
                   NeighborSearchUGParams *d_nsParams) {
        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        auto p_i = d_nsParams->particleIndices_cuData[i];
        auto neib_ind = p_i * d_nsConst->maxNeighborNum;

        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
            return;

        auto pos_i = d_data->predictPos[p_i];
        auto vel_i = d_data->vel_mix[p_i];
        Mat33f vGrad_sum;

        for (unsigned int p_j = d_nsParams->neighbors_cuData[neib_ind], t = 0;
             p_j != UINT_MAX && t < d_nsConst->maxNeighborNum;
             ++t, p_j = d_nsParams->neighbors_cuData[neib_ind + t]) {

            if (d_data->mat[p_j] != JL21CT_NONNEWTON)
                continue;

            auto pos_j = d_data->predictPos[p_j];
            auto vel_j = d_data->vel_mix[p_j];
            auto vel_ji = vel_j - vel_i;
            auto wGrad = cubic_gradient(pos_i - pos_j, d_const->sph_h);
            auto volume_j = d_data->mass[p_j] / d_data->density_sph[p_j];

            vGrad_sum += volume_j * vel_ji * wGrad;
        }

        d_data->vel_grad[p_i] = vGrad_sum;
    }

    extern __global__ void
    computeStress(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams) {
        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        auto p_i = d_nsParams->particleIndices_cuData[i];

        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
            return;

        // elastoplastic stress
        Mat33f tao0 = d_data->stress[p_i] - d_data->stress[p_i].trace() / 3 * Mat33f::eye();
        Mat33f e = (d_data->vel_grad[p_i] + d_data->vel_grad[p_i].transpose()) / 2;
        Mat33f e0 = e - e.trace() / 3 * Mat33f::eye();
        Mat33f w = (d_data->vel_grad[p_i] - d_data->vel_grad[p_i].transpose()) / 2;
        Mat33f stress_rate = w * tao0 - tao0 * w + 2 * d_const->shear_modulus * e0;

        d_data->stress[p_i] += stress_rate * d_const->dt * d_data->alpha[p_i].y;
        tao0 = d_data->stress[p_i] - d_data->stress[p_i].trace() / 3 * Mat33f::eye();

        if (d_const->Y < sqrtf(dot(tao0, tao0)))
            d_data->stress[p_i] /= d_const->Y;

        // fluid stress
        float vis_fluid = 0.005;

        Mat33f D = d_data->vel_grad[p_i] + d_data->vel_grad[p_i].transpose();
        Mat33f D0 = D - D.trace() / 3 * Mat33f::eye();
        d_data->stress[p_i] += 2 * vis_fluid * D0 * d_data->alpha[p_i].x;
    }

    extern __global__ void
    computeFluidStress(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
                       NeighborSearchUGParams *d_nsParams) {
        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        auto p_i = d_nsParams->particleIndices_cuData[i];

        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
            return;

        float vis_fluid = 0.001;

        Mat33f D = d_data->vel_grad[p_i] + d_data->vel_grad[p_i].transpose();
        Mat33f D0 = D - D.trace() / 3 * Mat33f::eye();
        d_data->stress[p_i] += (-d_data->pressure[p_i] * Mat33f::eye() + 2 * vis_fluid * D0) * d_data->alpha[p_i].x;
    }

    extern __global__ void
    computeElastoplasticStress(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
                               NeighborSearchUGParams *d_nsParams) {
        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        auto p_i = d_nsParams->particleIndices_cuData[i];

        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
            return;

        Mat33f tao0 = d_data->stress[p_i] - d_data->stress[p_i].trace() / 3 * Mat33f::eye();
        Mat33f e = (d_data->vel_grad[p_i] + d_data->vel_grad[p_i].transpose()) / 2;
        Mat33f e0 = e - e.trace() / 3 * Mat33f::eye();
        Mat33f w = (d_data->vel_grad[p_i] - d_data->vel_grad[p_i].transpose()) / 2;
        Mat33f stress_rate = w * tao0 - tao0 * w + 2 * d_const->shear_modulus * e0;

        d_data->stress[p_i] +=
                (-d_data->pressure[p_i] * Mat33f::eye() + stress_rate) * d_const->dt * d_data->alpha[p_i].y;
        tao0 = d_data->stress[p_i] - d_data->stress[p_i].trace() / 3 * Mat33f::eye();

        if (d_const->Y < sqrtf(dot(tao0, tao0)))
            d_data->stress[p_i] /= d_const->Y;
    }

    extern __global__ void
    computeStressForce(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
                       NeighborSearchUGParams *d_nsParams) {
        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        auto p_i = d_nsParams->particleIndices_cuData[i];
        auto neib_ind = p_i * d_nsConst->maxNeighborNum;

        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
            return;

        auto pos_i = d_data->predictPos[p_i];
        auto dens_i = d_data->density_sph[p_i];
        Vec3f acc;

        for (unsigned int p_j = d_nsParams->neighbors_cuData[neib_ind], t = 0;
             p_j != UINT_MAX && t < d_nsConst->maxNeighborNum;
             ++t, p_j = d_nsParams->neighbors_cuData[neib_ind + t]) {

            if (d_data->mat[p_j] != JL21CT_NONNEWTON)
                continue;

            auto pos_j = d_data->predictPos[p_j];
            auto wGrad = cubic_gradient(pos_i - pos_j, d_const->sph_h);
            auto dens_j = d_data->density_sph[p_j];

            acc += (d_data->stress[p_i] / powf(dens_i, 2) + d_data->stress[p_j] / powf(dens_j, 2)) * wGrad;
        }

        d_data->acc[p_i] += acc * dens_i;
    }

    extern __global__ void
    updateConformationTensor(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
                             NeighborSearchUGParams *d_nsParams) {
        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        auto p_i = d_nsParams->particleIndices_cuData[i];

        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
            return;

        Mat33f dQ = (d_data->Q[p_i] * d_data->vel_grad[p_i] + d_data->vel_grad[p_i].transpose() * d_data->Q[p_i] -
                     1 / d_const->relation_time * (d_data->Q[p_i] - Mat33f::eye())) * d_const->dt;
        d_data->Q[p_i] += dQ;
    }

    extern __global__ void
    computeElasticStress(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams) {
        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        auto p_i = d_nsParams->particleIndices_cuData[i];

        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
            return;

        //    d_data->stress[p_i] = 1 / d_const->compliance * (d_data->Q[p_i] - Mat33f::eye());
        d_data->stress[p_i] = 1 / d_data->compliance[p_i] * (d_data->Q[p_i] - Mat33f::eye());
    }

    extern __global__ void
    computeCompliance(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams) {
        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        auto p_i = d_nsParams->particleIndices_cuData[i];

        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
            return;

        // compute shear rate
        auto D = 0.5 * (d_data->vel_grad[p_i] + d_data->vel_grad[p_i].transpose());
        auto shearRate = sqrtf(0.5 * D.trace() * D.trace());
        d_data->shear_rate[p_i] = shearRate;

        // compute compliance
        auto scale_k = pow(1 - d_data->alpha[p_i].y, -2);

        if (d_data->alpha[p_i].y < d_const->alpha_low_threshold)
            d_data->compliance[p_i] = d_const->compliance_0;
        else if (d_data->alpha[p_i].y == 1)
            d_data->compliance[p_i] = d_const->compliance_inf;
        else {// non-newton
            if (d_data->alpha[p_i].y <= d_const->alpha_up_threshold)
                // thinning
                d_data->compliance[p_i] = d_const->compliance_inf + (d_const->compliance_0 - d_const->compliance_inf) /
                                                                    (1 + pow(30 * shearRate, -2)) / scale_k;
            else
                // thickening
                d_data->compliance[p_i] = d_const->compliance_inf + (d_const->compliance_0 - d_const->compliance_inf) /
                                                                    (1 + pow(50 * shearRate, 10)) / scale_k;
        }
    }

    extern __global__ void
    computeCohesionForce(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
                         NeighborSearchUGParams *d_nsParams) {
        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        auto p_i = d_nsParams->particleIndices_cuData[i];
        auto neib_ind = p_i * d_nsConst->maxNeighborNum;

        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
            return;

        Vec3f acc = {0, 0, 0};
        auto pos_i = d_data->predictPos[p_i];
        float gamma = 0.001;

        for (unsigned int p_j = d_nsParams->neighbors_cuData[neib_ind], t = 0;
             p_j != UINT_MAX && t < d_nsConst->maxNeighborNum;
             ++t, p_j = d_nsParams->neighbors_cuData[neib_ind + t]) {
            if (d_data->mat[p_j] != JL21CT_NONNEWTON)
                continue;

            auto pos_j = d_data->predictPos[p_j];
            acc += -gamma * d_data->mass[p_i] * d_data->mass[p_j] *
                   surface_tension_C((pos_i - pos_j).length(), d_const->sph_h) * (pos_i - pos_j) /
                   (pos_i - pos_j).length();
        }

        d_data->acc[p_i] += acc;
    }

    extern __global__ void
    cylinderRotate(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data) {
        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        if (d_data->mat[i] != ROTATE_RIGID)
            return;

        const float M_PI = 3.1415926;
        float angleRadians = 0.05 * (M_PI / 180.0f);// 将角度转换为弧度
        float cosAngle = cos(angleRadians);
        float sinAngle = sin(angleRadians);

        auto pos = d_data->predictPos[i];
        d_data->predictPos[i].x = pos.x * cosAngle - pos.z * sinAngle;
        d_data->predictPos[i].z = pos.x * sinAngle + pos.z * cosAngle;
        d_data->pos[i] = d_data->predictPos[i];
    }

    extern __global__ void
    moveStirringRod(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data) {
        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        if (d_data->mat[i] != ROTATE_RIGID)
            return;

        Vec3f offset = {0, 0, 0};

        if (d_const->cur_sim_time > 0.5f && d_const->cur_sim_time < 0.698)
            offset = {0, -8.2, 0};
        else if (d_const->cur_sim_time >= 0.68 && d_const->cur_sim_time < 2.1)
            offset = {0, -0.6, 0};
        else if (d_const->cur_sim_time > 25.5 && d_const->cur_sim_time < 30)
            offset = {0, 1, 0};

        d_data->predictPos[i] += offset * d_const->dt;
        d_data->pos[i] = d_data->predictPos[i];

        if (d_const->cur_sim_time >= 1.38 && d_const->cur_sim_time < 25.5) {// 25.5
            const float M_PI = 3.1415926;
            float angleRadians = 0.08 * (M_PI / 180.0f);// 将角度转换为弧度
            float cosAngle = cos(angleRadians);
            float sinAngle = sin(angleRadians);

            auto pos = d_data->predictPos[i];
            d_data->predictPos[i].x = pos.x * cosAngle - pos.z * sinAngle;
            d_data->predictPos[i].z = pos.x * sinAngle + pos.z * cosAngle;
            d_data->pos[i] = d_data->predictPos[i];
        }
    }

    extern __global__ void
    moveStirringRod_test(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data) {
        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        if (d_data->mat[i] != ROTATE_RIGID)
            return;

        Vec3f offset = {0, 0, 0};

        if (d_const->cur_sim_time > 0.5f && d_const->cur_sim_time < 0.71)
            offset = {0, -8.2, 0};
        else if (d_const->cur_sim_time >= 0.71 && d_const->cur_sim_time < 3.75)
            offset = {0, -0.6, 0};

        d_data->predictPos[i] += offset * d_const->dt;
        d_data->pos[i] = d_data->predictPos[i];

        if (d_const->cur_sim_time >= 3.6 && d_const->cur_sim_time < 30) {// 25.5
            const float M_PI = 3.1415926;
            float angleRadians = 0.075 * (M_PI / 180.0f);// 将角度转换为弧度
            float cosAngle = cos(angleRadians);
            float sinAngle = sin(angleRadians);

            auto pos = d_data->predictPos[i];
            d_data->predictPos[i].x = pos.x * cosAngle - pos.z * sinAngle;
            d_data->predictPos[i].z = pos.x * sinAngle + pos.z * cosAngle;
            d_data->pos[i] = d_data->predictPos[i];
        }
    }

    extern __global__ void
    moveBowl(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data) {
        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        if (d_data->mat[i] != BOUND_BOWL)
            return;

        if (d_const->cur_sim_time < 28)
            return;

        if (d_const->cur_sim_time >= 28 && d_const->cur_sim_time < 29.1) {
            const float M_PI = 3.1415926;
            float angleRadians = 0.075 * (M_PI / 180.0f);// 将角度转换为弧度
            float cosAngle = cos(angleRadians);
            float sinAngle = sin(angleRadians);

            auto pos = d_data->predictPos[i];
            d_data->predictPos[i].y = pos.y * cosAngle - pos.z * sinAngle;
            d_data->predictPos[i].z = pos.y * sinAngle + pos.z * cosAngle;
            d_data->pos[i] = d_data->predictPos[i];
        }
    }

    extern __global__ void
    matTransfer(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data) {
        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        if (d_data->mat[i] != FLUID_PREPARE_1)
            return;

        if (d_const->cur_sim_time > 2.2)
            d_data->mat[i] = JL21CT_NONNEWTON;
    }

    extern __global__ void
    moveRigid(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data) {
        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= d_const->total_particle_num)
            return;

        if (d_data->mat[i] != MOVE_RIGID)
            return;

        Vec3f offset = {1.2, 0, 0};

        d_data->predictPos[i] += offset * d_const->dt;
        d_data->pos[i] = d_data->predictPos[i];
    }

    extern __global__ void
    updateTime(JL21CTConstantParams *d_const) {
        d_const->cur_sim_time += d_const->dt;
    }
}