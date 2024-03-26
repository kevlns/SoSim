////
//// Created by ADMIN on 2024/3/22.
////
//
//#include "jl21_ct_cuda_api.cuh"
//#include "libs/SPHKernelL/kernels.cuh"
//
//namespace SoSim {
//
//    __global__ void initOtherData(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data) {
//        unsigned i = blockDim.x * blockIdx.x + threadIdx.x;
//        if (i >= d_const->total_particle_num)
//            return;
//
//        //    d_data->vel_mix[i] *= 0;
//        //    d_data->vel_k1[i] *= 0;
//        //    d_data->vel_k2[i] *= 0;
//        d_data->vel_drift_k1[i] *= 0;
//        d_data->vel_drift_k2[i] *= 0;
//        d_data->acc[i] *= 0;
//        d_data->d_alpha[i] *= 0;
//        d_data->lambda[i] = 1;
//        d_data->density_sph[i] *= 0;
//        d_data->pressure[i] *= 0;
//        d_data->density_mix[i] = dot(d_const->rest_density, d_data->alpha[i]);
//        if (d_data->mat[i] != JL21CT_NONNEWTON) {
//            d_data->density_mix[i] = 1100;
//        }
//        d_data->mass[i] = d_data->density_mix[i] * d_const->rest_volume;
//        d_data->color[i] = d_data->alpha[i].x * d_const->phase1_color + d_data->alpha[i].y * d_const->phase2_color;
//        d_data->Q[i] = Mat33f::eye();
//    }
//
//    extern __global__ void
//    estimateDensity(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
//                    NeighborSearchUGParams *d_nsParams) {
//        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
//        if (i >= d_const->total_particle_num)
//            return;
//
//        auto p_i = d_nsParams->particleIndices_cuData[i];
//        auto neib_ind = p_i * d_nsConst->maxNeighborNum;
//
//        //    if (d_data->mat[p_i] != JL21CT_NONNEWTON)
//        //        return;
//
//        d_data->density_sph[p_i] *= 0.0;
//        auto pos_i = d_data->predictPos[p_i];
//
//        for (unsigned int p_j = d_nsParams->neighbors_cuData[neib_ind], t = 0;
//             p_j != UINT_MAX && t < d_nsConst->maxNeighborNum;
//             ++t, p_j = d_nsParams->neighbors_cuData[neib_ind + t]) {
//
//            auto pos_j = d_data->predictPos[p_j];
//
//            d_data->density_sph[p_i] += d_data->mass[p_j] * cubic_value((pos_i - pos_j).length(), d_const->sph_h);
//
//            //        if (d_data->mat[p_j] == JL21CT_NONNEWTON)
//            //            d_data->density_sph[p_i] += d_data->mass[p_j] * cubic_value((pos_i - pos_j).length(), d_const->sph_h);
//            //        else if (d_data->mat[p_j] != JL21CT_NONNEWTON)
//            //            d_data->density_sph[p_i] +=
//            //                    d_data->density_mix[p_i] * d_data->volume[p_j] * cubic_value((pos_i - pos_j).length(), d_const->sph_h);
//        }
//
//        d_data->density_sph[p_i] = max(d_data->density_mix[p_i], d_data->density_sph[p_i]);
//        d_data->pressure[p_i] = d_const->wc_stiffness * d_data->density_mix[p_i] *
//                                (pow(d_data->density_sph[p_i] / d_data->density_mix[p_i], 7) - 1);
//    }
//
//    __global__ void
//    computeVelGrad(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
//                   NeighborSearchUGParams *d_nsParams) {
//        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
//        if (i >= d_const->total_particle_num)
//            return;
//
//        auto p_i = d_nsParams->particleIndices_cuData[i];
//        auto neib_ind = p_i * d_nsConst->maxNeighborNum;
//
//        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
//            return;
//
//        auto pos_i = d_data->predictPos[p_i];
//        auto vel_i = d_data->vel_mix[p_i];
//        Mat33f vGrad_sum;
//
//        for (unsigned int p_j = d_nsParams->neighbors_cuData[neib_ind], t = 0;
//             p_j != UINT_MAX && t < d_nsConst->maxNeighborNum;
//             ++t, p_j = d_nsParams->neighbors_cuData[neib_ind + t]) {
//
//            if (d_data->mat[p_j] != JL21CT_NONNEWTON)
//                continue;
//
//            auto pos_j = d_data->predictPos[p_j];
//            auto vel_j = d_data->vel_mix[p_j];
//            auto vel_ji = vel_j - vel_i;
//            auto wGrad = cubic_gradient(pos_i - pos_j, d_const->sph_h);
//            auto volume_j = d_data->mass[p_j] / d_data->density_sph[p_j];
//
//            vGrad_sum += volume_j * vel_ji * wGrad;
//        }
//
//        d_data->vel_grad[p_i] = vGrad_sum;
//    }
//
//    __global__ void
//    computeCompliance(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams) {
//        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
//        if (i >= d_const->total_particle_num)
//            return;
//
//        auto p_i = d_nsParams->particleIndices_cuData[i];
//
//        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
//            return;
//
//        // compute shear rate
//        auto D = 0.5 * (d_data->vel_grad[p_i] + d_data->vel_grad[p_i].transpose());
//        auto shearRate = sqrtf(0.5 * D.trace() * D.trace());
//        d_data->shear_rate[p_i] = shearRate;
//
//        // compute compliance
//        auto scale_k = pow(1 - d_data->alpha[p_i].y, -2);
//
//        if (d_data->alpha[p_i].y < d_const->alpha_low_threshold)
//            d_data->compliance[p_i] = d_const->compliance_0;
//        else if (d_data->alpha[p_i].y == 1)
//            d_data->compliance[p_i] = d_const->compliance_inf;
//        else {// non-newton
//            if (d_data->alpha[p_i].y <= d_const->alpha_up_threshold)
//                // thinning
//                d_data->compliance[p_i] = d_const->compliance_inf + (d_const->compliance_0 - d_const->compliance_inf) /
//                                                                    (1 + pow(30 * shearRate, -2)) / scale_k;
//            else
//                // thickening
//                d_data->compliance[p_i] = d_const->compliance_inf + (d_const->compliance_0 - d_const->compliance_inf) /
//                                                                    (1 + pow(50 * shearRate, 10)) / scale_k;
//        }
//    }
//
//    __global__ void
//    updateConformationTensor(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
//                             NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams) {
//        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
//        if (i >= d_const->total_particle_num)
//            return;
//
//        auto p_i = d_nsParams->particleIndices_cuData[i];
//
//        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
//            return;
//
//        Mat33f dQ = (d_data->Q[p_i] * d_data->vel_grad[p_i] + d_data->vel_grad[p_i].transpose() * d_data->Q[p_i] -
//                     1 / d_const->relation_time * (d_data->Q[p_i] - Mat33f::eye())) * d_const->dt;
//        d_data->Q[p_i] += dQ;
//    }
//
//    __global__ void
//    computeElasticStress(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
//                         NeighborSearchUGParams *d_nsParams) {
//        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
//        if (i >= d_const->total_particle_num)
//            return;
//
//        auto p_i = d_nsParams->particleIndices_cuData[i];
//
//        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
//            return;
//
//        //    d_data->stress[p_i] = 1 / d_const->compliance * (d_data->Q[p_i] - Mat33f::eye());
//        d_data->stress[p_i] = 1 / d_data->compliance[p_i] * (d_data->Q[p_i] - Mat33f::eye());
//    }
//
//    __global__ void
//    computeStressForce(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
//                       NeighborSearchUGParams *d_nsParams) {
//        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
//        if (i >= d_const->total_particle_num)
//            return;
//
//        auto p_i = d_nsParams->particleIndices_cuData[i];
//        auto neib_ind = p_i * d_nsConst->maxNeighborNum;
//
//        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
//            return;
//
//        auto pos_i = d_data->predictPos[p_i];
//        auto dens_i = d_data->density_sph[p_i];
//        Vec3f acc;
//
//        for (unsigned int p_j = d_nsParams->neighbors_cuData[neib_ind], t = 0;
//             p_j != UINT_MAX && t < d_nsConst->maxNeighborNum;
//             ++t, p_j = d_nsParams->neighbors_cuData[neib_ind + t]) {
//
//            if (d_data->mat[p_j] != JL21CT_NONNEWTON)
//                continue;
//
//            auto pos_j = d_data->predictPos[p_j];
//            auto wGrad = cubic_gradient(pos_i - pos_j, d_const->sph_h);
//            auto dens_j = d_data->density_sph[p_j];
//
//            acc += (d_data->stress[p_i] / powf(dens_i, 2) + d_data->stress[p_j] / powf(dens_j, 2)) * wGrad;
//        }
//
//        d_data->acc[p_i] += acc * dens_i;
//    }
//
//    __global__ void
//    computeFluidVisForce(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
//                         NeighborSearchUGParams *d_nsParams) {
//        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
//        if (i >= d_const->total_particle_num)
//            return;
//
//        auto p_i = d_nsParams->particleIndices_cuData[i];
//        auto neib_ind = p_i * d_nsConst->maxNeighborNum;
//
//        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
//            return;
//
//        auto pos_i = d_data->predictPos[p_i];
//        auto vel_i = d_data->vel_mix[p_i];
//        Vec3f acc = {0, 0, 0};
//        float vis = 0.01;
//        float h2_001 = 0.001f * pow(d_const->sph_h, 2);
//
//        for (unsigned int p_j = d_nsParams->neighbors_cuData[neib_ind], t = 0;
//             p_j != UINT_MAX && t < d_nsConst->maxNeighborNum;
//             ++t, p_j = d_nsParams->neighbors_cuData[neib_ind + t]) {
//
//            auto pos_j = d_data->predictPos[p_j];
//            auto vel_j = d_data->vel_mix[p_j];
//
//            auto pi = -vis * min(0.f, dot(vel_i - vel_j, pos_i - pos_j)) /
//                      ((pos_i - pos_j).length() * (pos_i - pos_j).length() + h2_001);
//            acc += -d_data->mass[p_i] * d_data->mass[p_j] * pi * cubic_gradient(pos_i - pos_j, d_const->sph_h);
//        }
//
//        d_data->acc[p_i] += acc;
//    }
//
//    __global__ void
//    computeCohesionForce(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
//                         NeighborSearchUGParams *d_nsParams) {
//        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
//        if (i >= d_const->total_particle_num)
//            return;
//
//        auto p_i = d_nsParams->particleIndices_cuData[i];
//        auto neib_ind = p_i * d_nsConst->maxNeighborNum;
//
//        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
//            return;
//
//        Vec3f acc = {0, 0, 0};
//        auto pos_i = d_data->predictPos[p_i];
//        float gamma = 0.001;
//
//        for (unsigned int p_j = d_nsParams->neighbors_cuData[neib_ind], t = 0;
//             p_j != UINT_MAX && t < d_nsConst->maxNeighborNum;
//             ++t, p_j = d_nsParams->neighbors_cuData[neib_ind + t]) {
//            if (d_data->mat[p_j] != JL21CT_NONNEWTON)
//                continue;
//
//            auto pos_j = d_data->predictPos[p_j];
//            acc += -gamma * d_data->mass[p_i] * d_data->mass[p_j] *
//                   surface_tension_C((pos_i - pos_j).length(), d_const->sph_h) * (pos_i - pos_j) /
//                   (pos_i - pos_j).length();
//        }
//
//        d_data->acc[p_i] += acc;
//    }
//
//    __global__ void
//    computeFluidPressureForce(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
//                              NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams) {
//        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
//        if (i >= d_const->total_particle_num)
//            return;
//
//        auto p_i = d_nsParams->particleIndices_cuData[i];
//        auto neib_ind = p_i * d_nsConst->maxNeighborNum;
//
//        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
//            return;
//
//        auto pos_i = d_data->predictPos[p_i];
//        auto volume_i = d_data->mass[p_i] / d_data->density_sph[p_i];
//        Vec3f acc = {0, 0, 0};
//
//        for (unsigned int p_j = d_nsParams->neighbors_cuData[neib_ind], t = 0;
//             p_j != UINT_MAX && t < d_nsConst->maxNeighborNum;
//             ++t, p_j = d_nsParams->neighbors_cuData[neib_ind + t]) {
//
//            auto pos_j = d_data->predictPos[p_j];
//            auto volume_j = d_data->mass[p_j] / d_data->density_sph[p_j];
//
//            acc -= volume_i * volume_j * (d_data->pressure[p_i] + d_data->pressure[p_j]) *
//                   cubic_gradient(pos_i - pos_j, d_const->sph_h);
//        }
//
//        d_data->acc[p_i] += acc;
//    }
//
//    __global__ void
//    computePhaseVel(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams) {
//        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
//        if (i >= d_const->total_particle_num)
//            return;
//
//        auto p_i = d_nsParams->particleIndices_cuData[i];
//
//        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
//            return;
//
//        // phase 1
//        d_data->vel_k1[p_i] +=
//                (d_data->acc[p_i] / (d_const->rest_density.x * d_const->rest_volume) + d_const->gravity) * d_const->dt;
//
//        // phase 2
//        d_data->vel_k2[p_i] +=
//                (d_data->acc[p_i] / (d_const->rest_density.y * d_const->rest_volume) + d_const->gravity) * d_const->dt;
//    }
//
//    __global__ void
//    computeMixVel(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams) {
//        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
//        if (i >= d_const->total_particle_num)
//            return;
//
//        auto p_i = d_nsParams->particleIndices_cuData[i];
//
//        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
//            return;
//
//        d_data->vel_mix[p_i] = d_data->alpha[p_i].x * d_data->vel_k1[p_i] + d_data->alpha[p_i].y * d_data->vel_k2[p_i];
//        d_data->vel_drift_k1[p_i] = d_data->vel_k1[p_i] - d_data->vel_mix[p_i];
//        d_data->vel_drift_k2[p_i] = d_data->vel_k2[p_i] - d_data->vel_mix[p_i];
//    }
//
//    __global__ void
//    computeDragForce(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams) {
//        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
//        if (i >= d_const->total_particle_num)
//            return;
//
//        auto p_i = d_nsParams->particleIndices_cuData[i];
//
//        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
//            return;
//
//        // phase 1
//        d_data->vel_k1[p_i] -= d_const->kd * d_const->dt * d_data->density_mix[p_i] / d_const->rest_density.x *
//                               (d_data->vel_k1[p_i] - d_data->vel_mix[p_i]);
//
//        // phase 2
//        d_data->vel_k2[p_i] -= d_const->kd * d_const->dt * d_data->density_mix[p_i] / d_const->rest_density.y *
//                               (d_data->vel_k2[p_i] - d_data->vel_mix[p_i]);
//
//    }
//
//    __global__ void
//    advect(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams) {
//        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
//        if (i >= d_const->total_particle_num)
//            return;
//
//        if (d_data->mat[i] != JL21CT_NONNEWTON)
//            return;
//
//        d_data->predictPos[i] += d_data->vel_mix[i] * d_const->dt;
//        d_data->pos[i] = d_data->predictPos[i];
//    }
//
//    __global__ void
//    reset(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams) {
//        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
//        if (i >= d_const->total_particle_num)
//            return;
//
//        if (d_data->mat[i] != JL21CT_NONNEWTON)
//            return;
//
//        d_data->acc[i] *= 0;
//        d_data->lambda[i] = 1;
//        d_data->stress[i] *= 0;
//    }
//
//    __global__ void
//    computePhaseTransport(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
//                          NeighborSearchUGParams *d_nsParams) {
//        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
//        if (i >= d_const->total_particle_num)
//            return;
//
//        auto p_i = d_nsParams->particleIndices_cuData[i];
//        auto neib_ind = p_i * d_nsConst->maxNeighborNum;
//
//        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
//            return;
//
//        auto pos_i = d_data->predictPos[p_i];
//        float h2_001 = 0.001f * pow(d_const->sph_h, 2);
//
//        float d_alpha_k1 = 0;
//        float d_alpha_k2 = 0;
//
//        // phase transport
//        for (unsigned int p_j = d_nsParams->neighbors_cuData[neib_ind], t = 0;
//             p_j != UINT_MAX && t < d_nsConst->maxNeighborNum;
//             ++t, p_j = d_nsParams->neighbors_cuData[neib_ind + t]) {
//
//            if (d_data->mat[p_j] != JL21CT_NONNEWTON)
//                continue;
//
//            auto pos_j = d_data->predictPos[p_j];
//
//            auto lam_ij = min(d_data->lambda[p_i], d_data->lambda[p_j]);
//
//            // compute T_m
//            auto T_m_k1 = -d_const->rest_volume * dot(d_data->alpha[p_i].x * d_data->vel_drift_k1[p_i] +
//                                                      d_data->alpha[p_j].x * d_data->vel_drift_k1[p_j],
//                                                      cubic_gradient(pos_i - pos_j, d_const->sph_h));
//            auto T_d_k1 = -2 * d_const->Cd * d_const->rest_volume * (d_data->alpha[p_i].x - d_data->alpha[p_j].x) *
//                          dot(pos_i - pos_j, cubic_gradient(pos_i - pos_j, d_const->sph_h)) /
//                          ((pos_i - pos_j).length() * (pos_i - pos_j).length() + h2_001);
//            d_alpha_k1 += d_const->dt * lam_ij * (T_m_k1 + T_d_k1);
//
//            auto T_m_k2 = -d_const->rest_volume * dot(d_data->alpha[p_i].y * d_data->vel_drift_k2[p_i] +
//                                                      d_data->alpha[p_j].y * d_data->vel_drift_k2[p_j],
//                                                      cubic_gradient(pos_i - pos_j, d_const->sph_h));
//            auto T_d_k2 = -2 * d_const->Cd * d_const->rest_volume * (d_data->alpha[p_i].y - d_data->alpha[p_j].y) *
//                          dot(pos_i - pos_j, cubic_gradient(pos_i - pos_j, d_const->sph_h)) /
//                          ((pos_i - pos_j).length() * (pos_i - pos_j).length() + h2_001);
//            d_alpha_k2 += d_const->dt * lam_ij * (T_m_k2 + T_d_k2);
//        }
//
//        d_data->d_alpha[p_i] = {d_alpha_k1, d_alpha_k2};
//    }
//
//    __global__ void
//    updateLambda(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams) {
//        uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
//        if (i >= d_const->total_particle_num)
//            return;
//
//        auto p_i = d_nsParams->particleIndices_cuData[i];
//
//        if (d_data->mat[p_i] != JL21CT_NONNEWTON)
//            return;
//
//        auto alpha_tmp = d_data->alpha[p_i] + d_data->d_alpha[p_i];
//        float min_lam = 1.0;
//        if (alpha_tmp.x < 1e-6)
//            min_lam = min(min_lam, -d_data->alpha[p_i].x / d_data->d_alpha[p_i].x);
//        if (alpha_tmp.y < 1e-6)
//            min_lam = min(min_lam, -d_data->alpha[p_i].y / d_data->d_alpha[p_i].y);
//
//        d_data->lambda[p_i] *= fabs(min_lam);
//    }
//
//    __global__ void
//    phaseTransportFinalize(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
//                           NeighborSearchUGParams *d_nsParams) {
//        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
//        if (i >= d_const->total_particle_num)
//            return;
//
//        if (d_data->mat[i] != JL21CT_NONNEWTON)
//            return;
//
//        // update volume fraction
//        d_data->alpha[i] += d_data->d_alpha[i];
//
//        // normalize frac
//        float frac_sum = 1e-6;
//        if (d_data->alpha[i].x < 0)
//            d_data->alpha[i].x = 0;
//        else
//            frac_sum += d_data->alpha[i].x;
//        if (d_data->alpha[i].y < 0)
//            d_data->alpha[i].y = 0;
//        else
//            frac_sum += d_data->alpha[i].y;
//
//        d_data->alpha[i] /= frac_sum;
//
//        if (d_data->alpha[i].x <= 0 && d_data->alpha[i].y <= 0)
//            d_data->alpha[i] = d_data->alpha_last[i];
//
//        // update density_mix and mass
//        d_data->density_mix[i] = dot(d_data->alpha[i], d_const->rest_density);
//        d_data->mass[i] = d_data->density_mix[i] * d_const->rest_volume;
//        d_data->alpha_last[i] = d_data->alpha[i];
//    }
//
//    __global__ void
//    updateColor(JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data, NeighborSearchUGParams *d_nsParams) {
//        unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
//        if (i >= d_const->total_particle_num)
//            return;
//
//        if (d_data->mat[i] != JL21CT_NONNEWTON)
//            return;
//
//        d_data->color[i] = d_data->alpha[i].x * d_const->phase1_color + d_data->alpha[i].y * d_const->phase2_color;
//    }
//
//}
//
///**
// *  host impl
// */
//namespace SoSim {
//
//    __host__ void
//    InitOtherData(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data) {
//        initOtherData<<<h_const.block_num, h_const.thread_num>>>(
//                d_const,
//                d_data);
//    }
//
//    __host__ void
//    EstimateDensity(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
//                    NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams) {
//        estimateDensity<<<h_const.block_num, h_const.thread_num>>>(
//                d_const,
//                d_data,
//                d_nsConst,
//                d_nsParams);
//    }
//
//    __host__ void
//    ComputeFluidVisForce(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
//                         NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams) {
//        computeFluidVisForce<<<h_const.block_num, h_const.thread_num>>>(
//                d_const,
//                d_data,
//                d_nsConst,
//                d_nsParams);
//    }
//
//    __host__ void
//    ComputeFluidPressureForce(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
//                              NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams) {
//        computeFluidPressureForce<<<h_const.block_num, h_const.thread_num>>>(
//                d_const,
//                d_data,
//                d_nsConst,
//                d_nsParams);
//    }
//
//    __host__ void
//    ComputePhaseVel(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
//                    NeighborSearchUGParams *d_nsParams) {
//        computePhaseVel<<<h_const.block_num, h_const.thread_num>>>(
//                d_const,
//                d_data,
//                d_nsParams);
//    }
//
//    __host__ void
//    ComputeMixVel(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
//                  NeighborSearchUGParams *d_nsParams) {
//        computeMixVel<<<h_const.block_num, h_const.thread_num>>>(
//                d_const,
//                d_data,
//                d_nsParams);
//    }
//
//    __host__ void
//    ComputeDragForce(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
//                     NeighborSearchUGParams *d_nsParams) {
//        computeDragForce<<<h_const.block_num, h_const.thread_num>>>(
//                d_const,
//                d_data,
//                d_nsParams);
//    }
//
//    __host__ void
//    Advect(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
//           NeighborSearchUGParams *d_nsParams) {
//        advect<<<h_const.block_num, h_const.thread_num>>>(
//                d_const,
//                d_data,
//                d_nsParams);
//    }
//
//    __host__ void
//    Reset(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
//          NeighborSearchUGParams *d_nsParams) {
//        reset<<<h_const.block_num, h_const.thread_num>>>(
//                d_const,
//                d_data,
//                d_nsParams);
//    }
//
//    __host__ void
//    PhaseTransport(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
//                   NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams) {
//        for (int i = 0; i < 2; ++i) {
//            computePhaseTransport<<<h_const.block_num, h_const.thread_num>>>(
//                    d_const,
//                    d_data,
//                    d_nsConst,
//                    d_nsParams);
//
//            updateLambda<<<h_const.block_num, h_const.thread_num>>>(
//                    d_const,
//                    d_data,
//                    d_nsParams);
//        }
//
//        phaseTransportFinalize<<<h_const.block_num, h_const.thread_num>>>(
//                d_const,
//                d_data,
//                d_nsParams);
//    }
//
//    __host__ void
//    UpdateColor(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
//                NeighborSearchUGParams *d_nsParams) {
//        updateColor<<<h_const.block_num, h_const.thread_num>>>(
//                d_const,
//                d_data,
//                d_nsParams);
//    }
//
//    __host__ void
//    ComputeVelGrad(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
//                   NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams) {
//        computeVelGrad<<<h_const.block_num, h_const.thread_num>>>(
//                d_const,
//                d_data,
//                d_nsConst,
//                d_nsParams);
//    }
//
//    __host__ void
//    ComputeCompliance(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
//                      NeighborSearchUGParams *d_nsParams) {
//        computeCompliance<<<h_const.block_num, h_const.thread_num>>>(
//                d_const,
//                d_data,
//                d_nsParams);
//    }
//
//    __host__ void
//    ComputeConformationTensorForce(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const,
//                                   JL21CTDynamicParams *d_data, NeighborSearchUGConfig *d_nsConst,
//                                   NeighborSearchUGParams *d_nsParams) {
//        updateConformationTensor<<<h_const.block_num, h_const.thread_num>>>(
//                d_const,
//                d_data,
//                d_nsConst,
//                d_nsParams);
//
//        computeElasticStress<<<h_const.block_num, h_const.thread_num>>>(
//                d_const,
//                d_data,
//                d_nsParams);
//
//        computeStressForce<<<h_const.block_num, h_const.thread_num>>>(
//                d_const,
//                d_data,
//                d_nsConst,
//                d_nsParams);
//    }
//
//    __host__ void
//    ComputeCohesionForce(JL21CTConstantParams &h_const, JL21CTConstantParams *d_const, JL21CTDynamicParams *d_data,
//                         NeighborSearchUGConfig *d_nsConst, NeighborSearchUGParams *d_nsParams) {
//        computeCohesionForce<<<h_const.block_num, h_const.thread_num>>>(
//                d_const,
//                d_data,
//                d_nsConst,
//                d_nsParams);
//    }
//
//}