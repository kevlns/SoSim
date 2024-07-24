//
// Created by ADMIN on 2024/7/23.
//

#ifndef SOSIM_MATRIX_V2_HPP
#define SOSIM_MATRIX_V2_HPP

#include <cuda_runtime.h>
#include <vector>
#include <iostream>

namespace SoSim {

    /**
     * Matrix_Type: Vec2<T>
     * @tparam T: float, int, double...
     */
    template<typename T>
    struct Vec2 {
        T data[2];

        __host__ __device__
        const size_t size() const {
            return 2;
        }

        __host__ __device__
        Vec2() {
            memset(data, 0, sizeof(data));
        }


        __host__ __device__
        Vec2(const T &n1, const T &n2) {
            data[0] = n1;
            data[1] = n2;
        }

        __host__ __device__
        Vec2(const Vec2<T> &v2) {
            memcpy_s(data, sizeof(data), v2.data, sizeof(data));
        }

        __host__ __device__
        Vec2<T> &operator=(const Vec2<T> &v2) {
            memcpy_s(data, sizeof(data), v2.data, sizeof(data));
            return *this;
        }

        __host__ __device__
        T &operator[](const size_t &index) {
            return data[index];
        }

        __host__ __device__
        const T &operator[](const size_t &index) const {
            return data[index];
        }

        __host__ __device__
        inline Vec2<T> operator+(const Vec2<T> &v) {
            return Vec2<T>(data[0] + v[0], data[1] + v[1]);
        }

        // -
        __host__ __device__
        inline Vec2<T> operator-(const Vec2<T> &v) {
            return Vec2<T>(data[0] - v[0], data[1] - v[1]);
        }

        __host__ __device__
        inline Vec2<T> operator*(const T &value) {
            return Vec2<T>(data[0] * value, data[1] * value);
        }

        __host__ __device__
        inline Vec2<T> operator/(const T &value) {
            return Vec2<T>(data[0] / value, data[1] / value);
        }

        // +=
        __host__ __device__
        inline Vec2<T> &operator+=(const Vec2<T> &v) {
            for (int i = 0; i < 2; i++)
                data[i] += v[i];
            return *this;
        }

        // -=
        __host__ __device__
        inline Vec2<T> &operator-=(const Vec2<T> &v) {
            for (int i = 0; i < 2; i++)
                data[i] -= v[i];
            return *this;
        }

        // *=
        __host__ __device__
        inline Vec2<T> &operator*=(const T &value) {
            for (int i = 0; i < 2; i++)
                data[i] *= value;
            return *this;
        }

        // /=
        __host__ __device__
        inline Vec2<T> &operator/=(const T &value) {
            for (int i = 0; i < 2; i++)
                data[i] /= value;
            return *this;
        }

        // normalized
        __host__ __device__
        Vec2<T> normalized() const {
            float norm = sqrt(data[0] * data[0] + data[1] * data[1]);
            if (norm == 0)
                return Vec2<T>(0, 0);
            return Vec2<T>(data[0] / norm, data[1] / norm);
        }

        // norm
        __host__ __device__
        T norm() const {
            return sqrt(data[0] * data[0] + data[1] * data[1]);
        }

        // cout
        __host__ __device__
        friend std::ostream &operator<<(std::ostream &os, const Vec2<T> &v2) {
            os << "[" << v2.data[0] << ", " << v2.data[1] << "]";
            return os;
        }
    };

    /**
     * Matrix_Type: Vec3<T>
     * @tparam T: float, int, double...
     */
    template<typename T>
    struct Vec3 {
        T data[3];

        __host__ __device__
        const size_t size() const {
            return 3;
        }

        __host__ __device__
        Vec3() {
            memset(data, 0, sizeof(data));
        }

        __host__ __device__
        Vec3(const T &n1, const T &n2, const T &n3) {
            data[0] = n1;
            data[1] = n2;
            data[2] = n3;
        }

        __host__ __device__
        Vec3(const Vec3<T> &v3) {
            memcpy_s(data, sizeof(data), v3.data, sizeof(data));
        }

        __host__ __device__
        Vec3<T> &operator=(const Vec3<T> &v3) {
            memcpy_s(data, sizeof(data), v3.data, sizeof(data));
            return *this;
        }

        __host__ __device__
        T &operator[](const size_t &index) {
            return data[index];
        }

        __host__ __device__
        const T &operator[](const size_t &index) const {
            return data[index];
        }

        __host__ __device__
        inline Vec3<T> operator+(const Vec3<T> &v) {
            return Vec3<T>(data[0] + v[0], data[1] + v[1], data[2] + v[2]);
        }

        // -
        __host__ __device__
        inline Vec3<T> operator-(const Vec3<T> &v) {
            return Vec3<T>(data[0] - v[0], data[1] - v[1], data[2] - v[2]);
        }

        __host__ __device__
        inline Vec3<T> operator*(const T &value) {
            return Vec3<T>(data[0] * value, data[1] * value, data[2] * value);
        }

        __host__ __device__
        inline Vec3<T> operator/(const T &value) {
            return Vec3<T>(data[0] / value, data[1] / value, data[2] / value);
        }

        // +=
        __host__ __device__
        inline Vec3<T> &operator+=(const Vec3<T> &v) {
            for (int i = 0; i < 3; i++)
                data[i] += v[i];
            return *this;
        }

        // -=
        __host__ __device__
        inline Vec3<T> &operator-=(const Vec3<T> &v) {
            for (int i = 0; i < 3; i++)
                data[i] -= v[i];
            return *this;
        }

        // *=
        __host__ __device__
        inline Vec3<T> &operator*=(const T &value) {
            for (int i = 0; i < 3; i++)
                data[i] *= value;
            return *this;
        }

        // /=
        __host__ __device__
        inline Vec3<T> &operator/=(const T &value) {
            for (int i = 0; i < 3; i++)
                data[i] /= value;
            return *this;
        }

        // normalized
        __host__ __device__
        Vec3<T> normalized() const {
            float norm = sqrt(data[0] * data[0] + data[1] * data[1] + data[2] * data[2]);
            if (norm == 0)
                return Vec3<T>(0, 0, 0);
            return Vec3<T>(data[0] / norm, data[1] / norm, data[2] / norm);
        }

        // norm
        __host__ __device__
        T norm() const {
            return sqrt(data[0] * data[0] + data[1] * data[1] + data[2] * data[2]);
        }

        // cout
        __host__ __device__
        friend std::ostream &operator<<(std::ostream &os, const Vec3<T> &v3) {
            os << "[" << v3.data[0] << ", " << v3.data[1] << ", " << v3.data[2] << "]";
            return os;
        }
    };

    // 4维向量
    template<typename T>
    struct Vec4 {
        T data[4];

        __host__ __device__
        const size_t size() const {
            return 4;
        }

        __host__ __device__
        Vec4() {
            memset(data, 0, sizeof(data));
        }

        __host__ __device__
        Vec4(const T &n1, const T &n2, const T &n3, const T &n4) {
            data[0] = n1;
            data[1] = n2;
            data[2] = n3;
            data[3] = n4;
        }

        __host__ __device__
        Vec4(const Vec4<T> &v4) {
            memcpy_s(data, sizeof(data), v4.data, sizeof(data));
        }

        __host__ __device__
        Vec4<T> &operator=(const Vec4<T> &v4) {
            memcpy_s(data, sizeof(data), v4.data, sizeof(data));
            return *this;
        }

        __host__ __device__
        T &operator[](const size_t &index) {
            return data[index];
        }

        __host__ __device__
        const T &operator[](const size_t &index) const {
            return data[index];
        }

        // +
        __host__ __device__
        inline Vec4<T> operator+(const Vec4<T> &v) {
            return Vec4<T>(data[0] + v[0], data[1] + v[1], data[2] + v[2], data[3] + v[3]);
        }

        // -
        __host__ __device__
        inline Vec4<T> operator-(const Vec4<T> &v) {
            return Vec4<T>(data[0] - v[0], data[1] - v[1], data[2] - v[2], data[3] - v[3]);
        }

        __host__ __device__
        inline Vec4<T> operator*(const T &value) {
            return Vec4<T>(data[0] * value, data[1] * value, data[2] * value, data[3] * value);
        }

        __host__ __device__
        inline Vec4<T> operator/(const T &value) {
            return Vec4<T>(data[0] / value, data[1] / value, data[2] / value, data[3] / value);
        }

        // +=
        __host__ __device__
        inline Vec4<T> &operator+=(const Vec4<T> &v) {
            for (int i = 0; i < 4; i++)
                data[i] += v[i];
            return *this;
        }

        // -=
        __host__ __device__
        inline Vec4<T> &operator-=(const Vec4<T> &v) {
            for (int i = 0; i < 4; i++)
                data[i] -= v[i];
            return *this;
        }

        // *=
        __host__ __device__
        inline Vec4<T> &operator*=(const T &value) {
            for (int i = 0; i < 4; i++)
                data[i] *= value;
            return *this;
        }

        // /=
        __host__ __device__
        inline Vec4<T> &operator/=(const T &value) {
            for (int i = 0; i < 4; i++)
                data[i] /= value;
            return *this;
        }

        // normalized
        __host__ __device__
        Vec4<T> normalized() const {
            float norm = sqrt(data[0] * data[0] + data[1] * data[1] + data[2] * data[2] + data[3] * data[3]);
            if (norm == 0)
                return Vec4<T>(0, 0, 0, 0);
            return Vec4<T>(data[0] / norm, data[1] / norm, data[2] / norm, data[3] / norm);
        }

        // norm
        __host__ __device__
        T norm() const {
            return sqrt(data[0] * data[0] + data[1] * data[1] + data[2] * data[2] + data[3] * data[3]);
        }

        // cout
        __host__ __device__
        friend std::ostream &operator<<(std::ostream &os, const Vec4<T> &v4) {
            os << "[" << v4.data[0] << ", " << v4.data[1] << ", " << v4.data[2] << ", " << v4.data[3] << "]";
            return os;
        }
    };

    // 2x2矩阵
    template<typename T>
    struct Mat22 {

        __host__ __device__
        const size_t size() const {
            return 4;
        }

        __host__ __device__
        Mat22() {
            memset(data, 0, sizeof(data));
        }

        __host__ __device__
        Mat22(const T &n1, const T &n2, const T &n3, const T &n4) {
            data[0] = n1;
            data[1] = n2;
            data[2] = n3;
            data[3] = n4;
        }

        __host__ __device__
        Mat22(const Mat22<T> &m22) {
            memcpy_s(data, sizeof(data), m22.data, sizeof(data));
        }

        __host__ __device__
        Mat22<T> &operator=(const Mat22<T> &m22) {
            memcpy_s(data, sizeof(data), m22.data, sizeof(data));
            return *this;
        }

        // at(row, col)
        __host__ __device__
        T &at(const size_t &row, const size_t &col) {
            return data[row * 2 + col];
        }

        // eye()
        __host__ __device__
        static Mat22<T> eye() {
            Mat22<T> m22;
            m22.data[0] = 1;
            m22.data[3] = 1;
            return m22;
        }

        // transpose
        __host__ __device__
        Mat22<T> transpose() const {
            Mat22<T> m22;
            m22.data[0] = data[0];
            m22.data[1] = data[2];
            m22.data[2] = data[1];
            m22.data[3] = data[3];
            return m22;
        }

        // trace
        __host__ __device__
        T trace() const {
            return data[0] + data[3];
        }

        // zero
        __host__ __device__
        static Mat22<T> zero() {
            Mat22<T> m22;
            memset(m22.data, 0, sizeof(m22.data));
            return m22;
        }

        // +
        __host__ __device__
        inline Mat22<T> operator+(const Mat22<T> &m) {
            return Mat22<T>(data[0] + m.data[0], data[1] + m.data[1], data[2] + m.data[2], data[3] + m.data[3]);
        }

        // -
        __host__ __device__
        inline Mat22<T> operator-(const Mat22<T> &m) {
            return Mat22<T>(data[0] - m.data[0], data[1] - m.data[1], data[2] - m.data[2], data[3] - m.data[3]);
        }

        // *
        __host__ __device__
        inline Mat22<T> operator*(const Mat22<T> &m) {
            return Mat22<T>(data[0] * m.data[0] + data[1] * m.data[2], data[0] * m.data[1] + data[1] * m.data[3],
                            data[2] * m.data[0] + data[3] * m.data[2], data[2] * m.data[1] + data[3] * m.data[3]);
        }

        // *
        __host__ __device__
        inline Mat22<T> operator*(const T &value) {
            return Mat22<T>(data[0] * value, data[1] * value, data[2] * value, data[3] * value);
        }

        // *
        __host__ __device__
        inline Mat22<T> operator*(const Vec2<T> &v) {
            return Mat22<T>(data[0] * v[0] + data[1] * v[1], data[2] * v[0] + data[3] * v[1]);
        }

        // /
        __host__ __device__
        inline Mat22<T> operator/(const T &value) {
            return Mat22<T>(data[0] / value, data[1] / value, data[2] / value, data[3] / value);
        }

        // *=
        __host__ __device__
        inline Mat22<T> &operator*=(const T &value) {
            for (int i = 0; i < 4; i++)
                data[i] *= value;
            return *this;
        }

        // *=
        __host__ __device__
        inline Mat22<T> &operator*=(const Mat22<T> &m) {
            T temp[4];
            temp[0] = data[0] * m.data[0] + data[1] * m.data[2];
            temp[1] = data[0] * m.data[1] + data[1] * m.data[3];
            temp[2] = data[2] * m.data[0] + data[3] * m.data[2];
            temp[3] = data[2] * m.data[1] + data[3] * m.data[3];
            memcpy_s(data, sizeof(data), temp, sizeof(temp));
            return *this;
        }

        // *=
        __host__ __device__
        inline Mat22<T> &operator*=(const Vec2<T> &v) {
            T temp[2];
            temp[0] = data[0] * v[0] + data[1] * v[1];
            temp[1] = data[2] * v[0] + data[3] * v[1];
            memcpy_s(data, sizeof(data), temp, sizeof(temp));
            return *this;
        }

        // +=
        __host__ __device__
        inline Mat22<T> &operator+=(const Mat22<T> &m) {
            for (int i = 0; i < 4; i++)
                data[i] += m.data[i];
            return *this;
        }

        // -=
        __host__ __device__
        inline Mat22<T> &operator-=(const Mat22<T> &m) {
            for (int i = 0; i < 4; i++)
                data[i] -= m.data[i];
            return *this;
        }

        // /=
        __host__ __device__
        inline Mat22<T> &operator/=(const T &value) {
            for (int i = 0; i < 4; i++)
                data[i] /= value;
            return *this;
        }

        // cout
        __host__ __device__
        friend std::ostream &operator<<(std::ostream &os, const Mat22<T> &m22) {
            os << "[" << m22.data[0] << ", " << m22.data[1] << "\n";
            os << " " << m22.data[2] << ", " << m22.data[3] << "]";
            return os;
        }

    private:
        T data[4];
    };

    // 3x3矩阵
    template<typename T>
    struct Mat33 {

        __host__ __device__
        const size_t size() const {
            return 9;
        }

        __host__ __device__
        Mat33() {
            memset(data, 0, sizeof(data));
        }

        __host__ __device__
        Mat33(const T &n1, const T &n2, const T &n3, const T &n4, const T &n5, const T &n6, const T &n7, const T &n8,
              const T &n9) {
            data[0] = n1;
            data[1] = n2;
            data[2] = n3;
            data[3] = n4;
            data[4] = n5;
            data[5] = n6;
            data[6] = n7;
            data[7] = n8;
            data[8] = n9;
        }

        __host__ __device__
        Mat33(const Mat33<T> &m33) {
            memcpy_s(data, sizeof(data), m33.data, sizeof(data));
        }

        __host__ __device__
        Mat33<T> &operator=(const Mat33<T> &m33) {
            memcpy_s(data, sizeof(data), m33.data, sizeof(data));
            return *this;
        }

        __host__ __device__
        T &at(const size_t &row, const size_t &col) {
            return data[row * 3 + col];
        }

        // eye()
        __host__ __device__
        static Mat33<T> eye() {
            Mat33<T> m33;
            m33.data[0] = 1;
            m33.data[4] = 1;
            m33.data[8] = 1;
            return m33;
        }

        // transpose
        __host__ __device__
        Mat33<T> transpose() {
            Mat33<T> m33;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    m33.data[i * 3 + j] = data[j * 3 + i];
                }
            }
            return m33;
        }

        // trace
        __host__ __device__
        T trace() {
            return data[0] + data[4] + data[8];
        }

        // zero
        __host__ __device__
        static Mat33<T> zero() {
            Mat33<T> m33;
            memset(m33.data, 0, sizeof(m33.data));
            return m33;
        }

        // +
        __host__ __device__
        inline Mat33<T> operator+(const Mat33<T> &m) {
            Mat33<T> m33;
            for (int i = 0; i < 9; i++) {
                m33.data[i] = data[i] + m.data[i];
            }
            return m33;
        }

        // -
        __host__ __device__
        inline Mat33<T> operator-(const Mat33<T> &m33) {
            Mat33<T> m;
            for (int i = 0; i < 9; i++) {
                m.data[i] = data[i] - m33.data[i];
            }
            return m;
        }

        // *
        __host__ __device__
        Mat33<T> operator*(const Mat33<T> &m33) {
            Mat33<T> m;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    m.data[i * 3 + j] = 0;
                    for (int k = 0; k < 3; k++) {
                        m.data[i * 3 + j] += data[i * 3 + k] * m33.data[k * 3 + j];
                    }
                }
            }
            return m;
        }

        // *
        __host__ __device__
        inline Mat33<T> operator*(const T &value) {
            Mat33<T> m;
            for (int i = 0; i < 9; i++) {
                m.data[i] = data[i] * value;
            }
            return m;
        }

        // *
        __host__ __device__
        inline Mat33<T> operator*(const Vec3<T> &v) {
            Mat33<T> m;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    m.data[i * 3 + j] = 0;
                    for (int k = 0; k < 3; k++) {
                        m.data[i * 3 + j] += data[i * 3 + k] * v.data[k];
                    }
                }
            }
            return m;
        }

        // *=
        __host__ __device__
        inline Mat33<T> &operator*=(const T &value) {
            for (int i = 0; i < 9; i++) {
                data[i] *= value;
            }
            return *this;
        }

        // *=
        __host__ __device__
        inline Mat33<T> &operator*=(const Mat33<T> &m33) {
            Mat33<T> m;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    m.data[i * 3 + j] = 0;
                    for (int k = 0; k < 3; k++) {
                        m.data[i * 3 + j] += data[i * 3 + k] * m33.data[k * 3 + j];
                    }
                }
            }
            *this = m;
        }

        // *=
        __host__ __device__
        inline Mat33<T> &operator*=(const Vec3<T> &v) {
            Mat33<T> m;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    m.data[i * 3 + j] = 0;
                    for (int k = 0; k < 3; k++) {
                        m.data[i * 3 + j] += data[i * 3 + k] * v.data[k];
                    }
                }
            }
            *this = m;
        }

        // mat / value
        __host__ __device__
        inline Mat33<T> operator/(const T &value) {
            Mat33<T> m;
            for (int i = 0; i < 9; i++) {
                m.data[i] = data[i] / value;
            }
            return m;
        }

        // +=
        __host__ __device__
        inline Mat33<T> &operator+=(const Mat33<T> &m) {
            for (int i = 0; i < 9; i++) {
                data[i] += m.data[i];
            }
            return *this;
        }

        // -=
        __host__ __device__
        inline Mat33<T> &operator-=(const Mat33<T> &m) {
            for (int i = 0; i < 9; i++) {
                data[i] -= m.data[i];
            }
            return *this;
        }

        // /=
        __host__ __device__
        inline Mat33<T> &operator/=(const T &value) {
            for (int i = 0; i < 9; i++) {
                data[i] /= value;
            }
            return *this;
        }

        // cout
        __host__ __device__
        friend std::ostream &operator<<(std::ostream &os, const Mat33<T> &m33) {
            os << "[" << m33.data[0] << ", " << m33.data[1] << ", " << m33.data[2] << "\n";
            os << " " << m33.data[3] << ", " << m33.data[4] << ", " << m33.data[5] << "\n";
            os << " " << m33.data[6] << ", " << m33.data[7] << ", " << m33.data[8] << "]";
            return os;
        }

    private:
        T data[9];
    };

    // 4x4矩阵
    template<typename T>
    struct Mat44 {

        __host__ __device__
        const size_t size() const {
            return 16;
        }

        __host__ __device__
        Mat44() {
            memset(data, 0, sizeof(data));
        }

        __host__ __device__
        Mat44(const T& n1, const T& n2, const T& n3, const T& n4,
              const T& n5, const T& n6, const T& n7, const T& n8,
              const T& n9, const T& n10, const T& n11, const T& n12,
              const T& n13, const T& n14, const T& n15, const T& n16){
            data[0] = n1;
            data[1] = n2;
            data[2] = n3;
            data[3] = n4;
            data[4] = n5;
            data[5] = n6;
            data[6] = n7;
            data[7] = n8;
            data[8] = n9;
            data[9] = n10;
            data[10] = n11;
            data[11] = n12;
            data[12] = n13;
            data[13] = n14;
            data[14] = n15;
            data[15] = n16;
        }

        __host__ __device__
        Mat44(const Mat44<T> &m44) {
            memcpy_s(data, sizeof(data), m44.data, sizeof(data));
        }

        __host__ __device__
        Mat44<T> &operator=(const Mat44<T> &m44) {
            memcpy_s(data, sizeof(data), m44.data, sizeof(data));
            return *this;
        }

        // at(row, col)
        __host__ __device__
        T &at(const size_t &row, const size_t &col) {
            return data[row * 4 + col];
        }

        // eye()
        __host__ __device__
        static Mat44<T> eye() {
            Mat44<T> m44;
            m44.data[0] = 1;
            m44.data[5] = 1;
            m44.data[10] = 1;
            m44.data[15] = 1;
            return m44;
        }

        __host__ __device__
        static Mat44<T> zero() {
            Mat44<T> m44;
            memset(m44.data, 0, sizeof(m44.data));
            return m44;
        }

        // transpose
        __host__ __device__
        Mat44<T> transpose() {
            Mat44<T> m;
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    m.data[i * 4 + j] = data[j * 4 + i];
                }
            }
            return m;
        }

        // trace
        __host__ __device__
        T trace() {
            return data[0] + data[5] + data[10] + data[15];
        }

        // +
        __host__ __device__
        Mat44<T> operator+(const Mat44<T> &m44) {
            Mat44<T> m;
            for (int i = 0; i < 16; i++) {
                m.data[i] = data[i] + m44.data[i];
            }
            return m;
        }

        // -
        __host__ __device__
        Mat44<T> operator-(const Mat44<T> &m44) {
            Mat44<T> m;
            for (int i = 0; i < 16; i++) {
                m.data[i] = data[i] - m44.data[i];
            }
            return m;
        }

        // *
        __host__ __device__
        Mat44<T> operator*(const Mat44<T> &m44) {
            Mat44<T> m;
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    m.data[i * 4 + j] = 0;
                    for (int k = 0; k < 4; k++) {
                        m.data[i * 4 + j] += data[i * 4 + k] * m44.data[k * 4 + j];
                    }
                }
            }
        }

        __host__ __device__
        Mat44<T> operator*(const T &value) {
            Mat44<T> m;
            for (int i = 0; i < 16; i++) {
                m.data[i] = data[i] * value;
            }
            return m;
        }

        // *
        __host__ __device__
        Mat44<T> operator*(const Vec4<T> &v) {
            Mat44<T> m;
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    m.data[i * 4 + j] = 0;
                    for (int k = 0; k < 4; k++) {
                        m.data[i * 4 + j] += data[i * 4 + k] * v.data[k];
                    }
                }
            }
            return m;
        }

        // *=
        __host__ __device__
        Mat44<T> &operator*=(const T &value) {
            for (int i = 0; i < 16; i++) {
                data[i] *= value;
            }
            return *this;
        }

        // *=
        __host__ __device__
        Mat44<T> &operator*=(const Mat44<T> &m44) {
            Mat44<T> m;
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    m.data[i * 4 + j] = 0;
                    for (int k = 0; k < 4; k++) {
                        m.data[i * 4 + j] += data[i * 4 + k] * m44.data[k * 4 + j];
                    }
                }
            }
            *this = m;
        }

        // *=
        __host__ __device__
        Mat44<T> &operator*=(const Vec4<T> &v) {
            Mat44<T> m;
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    m.data[i * 4 + j] = 0;
                    for (int k = 0; k < 4; k++) {
                        m.data[i * 4 + j] += data[i * 4 + k] * v.data[k];
                    }
                }
            }
            *this = m;
        }

        // /
        __host__ __device__
        Mat44<T> operator/(const T &value) {
            Mat44<T> m;
            for (int i = 0; i < 16; i++) {
                m.data[i] = data[i] / value;
            }
            return m;
        }

        // /=
        __host__ __device__
        Mat44<T> &operator/=(const T &value) {
            for (int i = 0; i < 16; i++) {
                data[i] /= value;
            }
            return *this;
        }

        // +=
        __host__ __device__
        Mat44<T> &operator+=(const Mat44<T> &m44) {
            for (int i = 0; i < 16; i++) {
                data[i] += m44.data[i];
            }
            return *this;
        }

        // -=
        __host__ __device__
        Mat44<T> &operator-=(const Mat44<T> &m44) {
            for (int i = 0; i < 16; i++) {
                data[i] -= m44.data[i];
            }
            return *this;
        }

        // cout
        __host__ __device__
        friend std::ostream &operator<<(std::ostream &os, const Mat44<T> &m44) {
            os << "[" << m44.data[0] << ", " << m44.data[1] << ", " << m44.data[2] << ", " << m44.data[3] << "\n";
            os << " " << m44.data[4] << ", " << m44.data[5] << ", " << m44.data[6] << "," << m44.data[7] << "\n";
            os << " " << m44.data[8] << ", " << m44.data[9] << ", " << m44.data[10] << "," << m44.data[11] << "\n";
            os << " " << m44.data[12] << ", " << m44.data[13] << ", " << m44.data[14] << "," << m44.data[15] << "]";
            return os;
        }

    private:
        T data[16];
    };


    /**
     * Operators
     */

    // Vec2
    template<typename T>
    __host__ __device__
    inline Vec2<T> operator*(T &&value, Vec2<T> &v) {
        return v * value;
    }

    // dot
    template<typename T>
    __host__ __device__
    inline T dot(const Vec2<T> &v1, const Vec2<T> &v2) {
        return v1[0] * v2[0] + v1[1] * v2[1];
    }

    // outer product
    template<typename T>
    __host__ __device__
    inline Mat22<T> outer(const Vec2<T> &v1, const Vec2<T> &v2) {
        Mat22<T> m22;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                m22.at(i, j) = v1[i] * v2[j];
            }
        }
        return m22;
    }

    // Vec3
    template<typename T>
    __host__ __device__
    inline T dot(const Vec3<T> &v1, const Vec3<T> &v2) {
        return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    }

    template<typename T>
    __host__ __device__
    inline Mat33<T> outer(const Vec3<T> &v1, const Vec3<T> &v2) {
        Mat33<T> m33;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                m33.at(i, j) = v1[i] * v2[j];
            }
        }
        return m33;
    }

    template<typename T>
    __host__ __device__
    inline Vec3<T> operator*(T &&value, Vec3<T> &v) {
        return v * value;
    }

    template<typename T>
    __host__ __device__
    inline Vec3<T> cross(const Vec3<T> &v1, const Vec3<T> &v2) {
        return Vec3<T>(v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0]);
    }

    // Vec4
    template<typename T>
    __host__ __device__
    inline T dot(const Vec4<T> &v1, const Vec4<T> &v2) {
        return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3];
    }

    template<typename T>
    __host__ __device__
    inline Vec4<T> operator*(T &&value, Vec4<T> &v) {
        return v * value;
    }

    template<typename T>
    __host__ __device__
    inline Mat44<T> outer(const Vec4<T> &v1, const Vec4<T> &v2) {
        Mat44<T> m44;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                m44.at(i, j) = v1[i] * v2[j];
            }
        }
        return m44;
    }

    template<typename T>
    __host__ __device__
    inline Mat22<T> operator*(T &&value, Mat22<T> &m) {
        return m * value;
    }

    template<typename T>
    __host__ __device__
    inline Mat33<T> operator*(T &&value, Mat33<T> &m) {
        return m * value;
    }

    template<typename T>
    __host__ __device__
    inline Mat44<T> operator*(T &&value, Mat44<T> &m) {
        return m * value;
    }
}

#endif //SOSIM_MATRIX_V2_HPP
