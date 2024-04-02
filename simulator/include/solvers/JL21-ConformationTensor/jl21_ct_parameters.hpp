//
// Created by ADMIN on 2024/3/22.
//

#ifndef SOSIM_JL21_CT_PARAMETERS_HPP
#define SOSIM_JL21_CT_PARAMETERS_HPP

#include <cuda_runtime.h>

#include "core/math/matrix.hpp"
#include "framework/MetaFramework/sim_material.hpp"

namespace SoSim {

    struct JL21CTConstantParams {
        float dt{0};
        float cur_sim_time{0};
        Vec3f gravity{0, -9.8f, 0};
        unsigned block_num;
        unsigned thread_num;

        // JL21
        Vec3f phase1_color;
        Vec3f phase2_color;
        Vec2f rest_density;
        float wc_stiffness;
        float rest_volume;
        float sph_h;
        float kd;
        float Cd;
        int total_particle_num;

        // conformation tensor
        float relation_time;
        float alpha_up_threshold;
        float alpha_low_threshold;
        float compliance_0;
        float compliance_inf;

        // suspension fluid
        float shear_modulus{10000};
        float Y{100};
    };

    struct JL21CTDynamicParams {
    public:
        Vec3f *vel_mix;
        Vec3f *vel_k1;
        Vec3f *vel_k2;
        Vec3f *vel_drift_k1;
        Vec3f *vel_drift_k2;
        Vec3f *acc;
        Vec3f *predictPos;
        Vec3f *pos;
        Vec3f *color;
        Vec2f *alpha;
        Vec2f *alpha_last;
        Vec2f *d_alpha;
        float *lambda;
        float *density_mix;
        float *density_sph;
        float *pressure;
        float *mass;
        float *volume;
        float *compliance;
        float *shear_rate;
        Material *mat;
        bool *active;

        // extra params
        Mat33f *vel_grad;
        Mat33f *elastic_stress;
        Mat33f *stress;
        Mat33f *Q;

    private:
        bool isInit{false};

    public:
        void malloc(int particle_num) {
            freeMemory();

            cudaMalloc((void **) &vel_mix, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &vel_k1, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &vel_k2, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &vel_drift_k1, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &vel_drift_k2, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &acc, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &predictPos, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &pos, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &color, particle_num * sizeof(Vec3f));
            cudaMalloc((void **) &alpha, particle_num * sizeof(Vec2f));
            cudaMalloc((void **) &alpha_last, particle_num * sizeof(Vec2f));
            cudaMalloc((void **) &d_alpha, particle_num * sizeof(Vec2f));
            cudaMalloc((void **) &lambda, particle_num * sizeof(float));
            cudaMalloc((void **) &density_mix, particle_num * sizeof(float));
            cudaMalloc((void **) &density_sph, particle_num * sizeof(float));
            cudaMalloc((void **) &pressure, particle_num * sizeof(float));
            cudaMalloc((void **) &mass, particle_num * sizeof(float));
            cudaMalloc((void **) &volume, particle_num * sizeof(float));
            cudaMalloc((void **) &compliance, particle_num * sizeof(float));
            cudaMalloc((void **) &shear_rate, particle_num * sizeof(float));
            cudaMalloc((void **) &mat, particle_num * sizeof(Material));
            cudaMalloc((void **) &active, particle_num * sizeof(bool));

            cudaMalloc((void **) &vel_grad, particle_num * sizeof(Mat33f));
            cudaMalloc((void **) &elastic_stress, particle_num * sizeof(Mat33f));
            cudaMalloc((void **) &stress, particle_num * sizeof(Mat33f));
            cudaMalloc((void **) &Q, particle_num * sizeof(Mat33f));

            cudaGetLastError();

            isInit = true;
        }

        void freeMemory() {
            if (isInit) {
                cudaFree(vel_mix);
                cudaFree(vel_k1);
                cudaFree(vel_k2);
                cudaFree(vel_drift_k1);
                cudaFree(vel_drift_k2);
                cudaFree(acc);
                cudaFree(predictPos);
                cudaFree(pos);
                cudaFree(color);
                cudaFree(alpha);
                cudaFree(alpha_last);
                cudaFree(d_alpha);
                cudaFree(lambda);
                cudaFree(density_mix);
                cudaFree(density_sph);
                cudaFree(pressure);
                cudaFree(mass);
                cudaFree(volume);
                cudaFree(compliance);
                cudaFree(shear_rate);
                cudaFree(mat);
                cudaFree(active);

                cudaFree(vel_grad);
                cudaFree(elastic_stress);
                cudaFree(stress);
                cudaFree(Q);

                cudaGetLastError();

                isInit = false;
            }
        }
    };

}

#endif //SOSIM_JL21_CT_PARAMETERS_HPP
