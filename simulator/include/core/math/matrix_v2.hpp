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
    private:
        T _data_[2];

    public:
        __host__ __device__
        const size_t size() const {
            return 2;
        }

        __host__ __device__
        Vec2() {
            memset(_data_, 0, sizeof(_data_));
        }


        __host__ __device__
        Vec2(const T &n1, const T &n2) {
            _data_[0] = n1;
            _data_[1] = n2;
        }

        __host__ __device__
        Vec2(const Vec2<T> &v2) {
            memcpy_s(_data_, sizeof(_data_), v2._data_, sizeof(_data_));
        }

        __host__ __device__
        Vec2<T> &operator=(const Vec2<T> &v2) {
            memcpy_s(_data_, sizeof(_data_), v2._data_, sizeof(_data_));
            return *this;
        }

        __host__ __device__
        T &operator[](const size_t &index) {
            return _data_[index];
        }

        __host__ __device__
        const T &operator[](const size_t &index) const {
            return _data_[index];
        }

        // data()
        __host__ __device__
        T *data() {
            return _data_;
        }

        __host__ __device__
        inline Vec2<T> operator+(const Vec2<T> &v) {
            return Vec2<T>(_data_[0] + v[0], _data_[1] + v[1]);
        }

        // -
        __host__ __device__
        inline Vec2<T> operator-(const Vec2<T> &v) {
            return Vec2<T>(_data_[0] - v[0], _data_[1] - v[1]);
        }

        __host__ __device__
        inline Vec2<T> operator*(const T &value) {
            return Vec2<T>(_data_[0] * value, _data_[1] * value);
        }

        __host__ __device__
        inline Vec2<T> operator/(const T &value) {
            return Vec2<T>(_data_[0] / value, _data_[1] / value);
        }

        // +=
        __host__ __device__
        inline Vec2<T> &operator+=(const Vec2<T> &v) {
            for (int i = 0; i < 2; i++)
                _data_[i] += v[i];
            return *this;
        }

        // -=
        __host__ __device__
        inline Vec2<T> &operator-=(const Vec2<T> &v) {
            for (int i = 0; i < 2; i++)
                _data_[i] -= v[i];
            return *this;
        }

        // *=
        __host__ __device__
        inline Vec2<T> &operator*=(const T &value) {
            for (int i = 0; i < 2; i++)
                _data_[i] *= value;
            return *this;
        }

        // /=
        __host__ __device__
        inline Vec2<T> &operator/=(const T &value) {
            for (int i = 0; i < 2; i++)
                _data_[i] /= value;
            return *this;
        }

        // normalized
        __host__ __device__
        Vec2<T> normalized() const {
            float norm = sqrt(_data_[0] * _data_[0] + _data_[1] * _data_[1]);
            if (norm == 0)
                return Vec2<T>(0, 0);
            return Vec2<T>(_data_[0] / norm, _data_[1] / norm);
        }

        // norm
        __host__ __device__
        T norm() const {
            return sqrt(_data_[0] * _data_[0] + _data_[1] * _data_[1]);
        }

        // cout
        __host__ __device__
        friend std::ostream &operator<<(std::ostream &os, const Vec2<T> &v2) {
            os << "[" << v2._data_[0] << ", " << v2._data_[1] << "]";
            return os;
        }
    };

    /**
     * Matrix_Type: Vec3<T>
     * @tparam T: float, int, double...
     */
    template<typename T>
    struct Vec3 {
    private:
        T _data_[3];

    public:
        __host__ __device__
        const size_t size() const {
            return 3;
        }

        __host__ __device__
        Vec3() {
            memset(_data_, 0, sizeof(_data_));
        }

        __host__ __device__
        Vec3(const T &n1, const T &n2, const T &n3) {
            _data_[0] = n1;
            _data_[1] = n2;
            _data_[2] = n3;
        }

        __host__ __device__
        Vec3(const Vec3<T> &v3) {
            memcpy_s(_data_, sizeof(_data_), v3._data_, sizeof(_data_));
        }

        __host__ __device__
        Vec3<T> &operator=(const Vec3<T> &v3) {
            memcpy_s(_data_, sizeof(_data_), v3._data_, sizeof(_data_));
            return *this;
        }

        __host__ __device__
        T &operator[](const size_t &index) {
            return _data_[index];
        }

        __host__ __device__
        const T &operator[](const size_t &index) const {
            return _data_[index];
        }

        // data()
        __host__ __device__
        T *data() {
            return _data_;
        }

        __host__ __device__
        inline Vec3<T> operator+(const Vec3<T> &v) {
            return Vec3<T>(_data_[0] + v[0], _data_[1] + v[1], _data_[2] + v[2]);
        }

        // -
        __host__ __device__
        inline Vec3<T> operator-(const Vec3<T> &v) {
            return Vec3<T>(_data_[0] - v[0], _data_[1] - v[1], _data_[2] - v[2]);
        }

        __host__ __device__
        inline Vec3<T> operator*(const T &value) {
            return Vec3<T>(_data_[0] * value, _data_[1] * value, _data_[2] * value);
        }

        __host__ __device__
        inline Vec3<T> operator/(const T &value) {
            return Vec3<T>(_data_[0] / value, _data_[1] / value, _data_[2] / value);
        }

        // +=
        __host__ __device__
        inline Vec3<T> &operator+=(const Vec3<T> &v) {
            for (int i = 0; i < 3; i++)
                _data_[i] += v[i];
            return *this;
        }

        // -=
        __host__ __device__
        inline Vec3<T> &operator-=(const Vec3<T> &v) {
            for (int i = 0; i < 3; i++)
                _data_[i] -= v[i];
            return *this;
        }

        // *=
        __host__ __device__
        inline Vec3<T> &operator*=(const T &value) {
            for (int i = 0; i < 3; i++)
                _data_[i] *= value;
            return *this;
        }

        // /=
        __host__ __device__
        inline Vec3<T> &operator/=(const T &value) {
            for (int i = 0; i < 3; i++)
                _data_[i] /= value;
            return *this;
        }

        // normalized
        __host__ __device__
        Vec3<T> normalized() const {
            float norm = sqrt(_data_[0] * _data_[0] + _data_[1] * _data_[1] + _data_[2] * _data_[2]);
            if (norm == 0)
                return Vec3<T>(0, 0, 0);
            return Vec3<T>(_data_[0] / norm, _data_[1] / norm, _data_[2] / norm);
        }

        // norm
        __host__ __device__
        T norm() const {
            return sqrt(_data_[0] * _data_[0] + _data_[1] * _data_[1] + _data_[2] * _data_[2]);
        }

        // cout
        __host__ __device__
        friend std::ostream &operator<<(std::ostream &os, const Vec3<T> &v3) {
            os << "[" << v3._data_[0] << ", " << v3._data_[1] << ", " << v3._data_[2] << "]";
            return os;
        }
    };

    // 4维向量
    template<typename T>
    struct Vec4 {
    private:
        T _data_[4];

    public:
        __host__ __device__
        const size_t size() const {
            return 4;
        }

        __host__ __device__
        Vec4() {
            memset(_data_, 0, sizeof(_data_));
        }

        __host__ __device__
        Vec4(const T &n1, const T &n2, const T &n3, const T &n4) {
            _data_[0] = n1;
            _data_[1] = n2;
            _data_[2] = n3;
            _data_[3] = n4;
        }

        __host__ __device__
        Vec4(const Vec4<T> &v4) {
            memcpy_s(_data_, sizeof(_data_), v4._data_, sizeof(_data_));
        }

        __host__ __device__
        Vec4<T> &operator=(const Vec4<T> &v4) {
            memcpy_s(_data_, sizeof(_data_), v4._data_, sizeof(_data_));
            return *this;
        }

        __host__ __device__
        T &operator[](const size_t &index) {
            return _data_[index];
        }

        __host__ __device__
        const T &operator[](const size_t &index) const {
            return _data_[index];
        }

        // data()
        __host__ __device__
        T *data() {
            return _data_;
        }

        // +
        __host__ __device__
        inline Vec4<T> operator+(const Vec4<T> &v) {
            return Vec4<T>(_data_[0] + v[0], _data_[1] + v[1], _data_[2] + v[2], _data_[3] + v[3]);
        }

        // -
        __host__ __device__
        inline Vec4<T> operator-(const Vec4<T> &v) {
            return Vec4<T>(_data_[0] - v[0], _data_[1] - v[1], _data_[2] - v[2], _data_[3] - v[3]);
        }

        __host__ __device__
        inline Vec4<T> operator*(const T &value) {
            return Vec4<T>(_data_[0] * value, _data_[1] * value, _data_[2] * value, _data_[3] * value);
        }

        __host__ __device__
        inline Vec4<T> operator/(const T &value) {
            return Vec4<T>(_data_[0] / value, _data_[1] / value, _data_[2] / value, _data_[3] / value);
        }

        // +=
        __host__ __device__
        inline Vec4<T> &operator+=(const Vec4<T> &v) {
            for (int i = 0; i < 4; i++)
                _data_[i] += v[i];
            return *this;
        }

        // -=
        __host__ __device__
        inline Vec4<T> &operator-=(const Vec4<T> &v) {
            for (int i = 0; i < 4; i++)
                _data_[i] -= v[i];
            return *this;
        }

        // *=
        __host__ __device__
        inline Vec4<T> &operator*=(const T &value) {
            for (int i = 0; i < 4; i++)
                _data_[i] *= value;
            return *this;
        }

        // /=
        __host__ __device__
        inline Vec4<T> &operator/=(const T &value) {
            for (int i = 0; i < 4; i++)
                _data_[i] /= value;
            return *this;
        }

        // normalized
        __host__ __device__
        Vec4<T> normalized() const {
            float norm = sqrt(_data_[0] * _data_[0] + _data_[1] * _data_[1] + _data_[2] * _data_[2] + _data_[3] * _data_[3]);
            if (norm == 0)
                return Vec4<T>(0, 0, 0, 0);
            return Vec4<T>(_data_[0] / norm, _data_[1] / norm, _data_[2] / norm, _data_[3] / norm);
        }

        // norm
        __host__ __device__
        T norm() const {
            return sqrt(_data_[0] * _data_[0] + _data_[1] * _data_[1] + _data_[2] * _data_[2] + _data_[3] * _data_[3]);
        }

        // cout
        __host__ __device__
        friend std::ostream &operator<<(std::ostream &os, const Vec4<T> &v4) {
            os << "[" << v4._data_[0] << ", " << v4._data_[1] << ", " << v4._data_[2] << ", " << v4._data_[3] << "]";
            return os;
        }
    };

    // 2x2矩阵
    template<typename T>
    struct Mat22 {
    private:
        T _data_[4];
        
    public:

        __host__ __device__
        const size_t size() const {
            return 4;
        }

        __host__ __device__
        Mat22() {
            memset(_data_, 0, sizeof(_data_));
        }

        __host__ __device__
        Mat22(const T &n1, const T &n2, const T &n3, const T &n4) {
            _data_[0] = n1;
            _data_[1] = n2;
            _data_[2] = n3;
            _data_[3] = n4;
        }

        __host__ __device__
        Mat22(const Mat22<T> &m22) {
            memcpy_s(_data_, sizeof(_data_), m22._data_, sizeof(_data_));
        }

        __host__ __device__
        Mat22<T> &operator=(const Mat22<T> &m22) {
            memcpy_s(_data_, sizeof(_data_), m22._data_, sizeof(_data_));
            return *this;
        }

        // at(row, col)
        __host__ __device__
        T &at(const size_t &row, const size_t &col) {
            return _data_[row * 2 + col];
        }

        // data()
        __host__ __device__
        T *data() {
            return _data_;
        }

        // eye()
        __host__ __device__
        static Mat22<T> eye() {
            Mat22<T> m22;
            m22._data_[0] = 1;
            m22._data_[3] = 1;
            return m22;
        }

        // transpose
        __host__ __device__
        Mat22<T> transpose() const {
            Mat22<T> m22;
            m22._data_[0] = _data_[0];
            m22._data_[1] = _data_[2];
            m22._data_[2] = _data_[1];
            m22._data_[3] = _data_[3];
            return m22;
        }

        // trace
        __host__ __device__
        T trace() const {
            return _data_[0] + _data_[3];
        }

        // zero
        __host__ __device__
        static Mat22<T> zero() {
            Mat22<T> m22;
            memset(m22._data_, 0, sizeof(m22._data_));
            return m22;
        }

        // +
        __host__ __device__
        inline Mat22<T> operator+(const Mat22<T> &m) {
            return Mat22<T>(_data_[0] + m._data_[0], _data_[1] + m._data_[1], _data_[2] + m._data_[2], _data_[3] + m._data_[3]);
        }

        // -
        __host__ __device__
        inline Mat22<T> operator-(const Mat22<T> &m) {
            return Mat22<T>(_data_[0] - m._data_[0], _data_[1] - m._data_[1], _data_[2] - m._data_[2], _data_[3] - m._data_[3]);
        }

        // *
        __host__ __device__
        inline Mat22<T> operator*(const Mat22<T> &m) {
            return Mat22<T>(_data_[0] * m._data_[0] + _data_[1] * m._data_[2], _data_[0] * m._data_[1] + _data_[1] * m._data_[3],
                            _data_[2] * m._data_[0] + _data_[3] * m._data_[2], _data_[2] * m._data_[1] + _data_[3] * m._data_[3]);
        }

        // *
        __host__ __device__
        inline Mat22<T> operator*(const T &value) {
            return Mat22<T>(_data_[0] * value, _data_[1] * value, _data_[2] * value, _data_[3] * value);
        }

        // *
        __host__ __device__
        inline Mat22<T> operator*(const Vec2<T> &v) {
            return Mat22<T>(_data_[0] * v[0] + _data_[1] * v[1], _data_[2] * v[0] + _data_[3] * v[1]);
        }

        // /
        __host__ __device__
        inline Mat22<T> operator/(const T &value) {
            return Mat22<T>(_data_[0] / value, _data_[1] / value, _data_[2] / value, _data_[3] / value);
        }

        // *=
        __host__ __device__
        inline Mat22<T> &operator*=(const T &value) {
            for (int i = 0; i < 4; i++)
                _data_[i] *= value;
            return *this;
        }

        // *=
        __host__ __device__
        inline Mat22<T> &operator*=(const Mat22<T> &m) {
            T temp[4];
            temp[0] = _data_[0] * m._data_[0] + _data_[1] * m._data_[2];
            temp[1] = _data_[0] * m._data_[1] + _data_[1] * m._data_[3];
            temp[2] = _data_[2] * m._data_[0] + _data_[3] * m._data_[2];
            temp[3] = _data_[2] * m._data_[1] + _data_[3] * m._data_[3];
            memcpy_s(_data_, sizeof(_data_), temp, sizeof(temp));
            return *this;
        }

        // *=
        __host__ __device__
        inline Mat22<T> &operator*=(const Vec2<T> &v) {
            T temp[2];
            temp[0] = _data_[0] * v[0] + _data_[1] * v[1];
            temp[1] = _data_[2] * v[0] + _data_[3] * v[1];
            memcpy_s(_data_, sizeof(_data_), temp, sizeof(temp));
            return *this;
        }

        // +=
        __host__ __device__
        inline Mat22<T> &operator+=(const Mat22<T> &m) {
            for (int i = 0; i < 4; i++)
                _data_[i] += m._data_[i];
            return *this;
        }

        // -=
        __host__ __device__
        inline Mat22<T> &operator-=(const Mat22<T> &m) {
            for (int i = 0; i < 4; i++)
                _data_[i] -= m._data_[i];
            return *this;
        }

        // /=
        __host__ __device__
        inline Mat22<T> &operator/=(const T &value) {
            for (int i = 0; i < 4; i++)
                _data_[i] /= value;
            return *this;
        }

        // cout
        __host__ __device__
        friend std::ostream &operator<<(std::ostream &os, const Mat22<T> &m22) {
            os << "[" << m22._data_[0] << ", " << m22._data_[1] << "\n";
            os << " " << m22._data_[2] << ", " << m22._data_[3] << "]";
            return os;
        }
    };

    // 3x3矩阵
    template<typename T>
    struct Mat33 {
    private:
        T _data_[9];

    public:
        __host__ __device__
        const size_t size() const {
            return 9;
        }

        __host__ __device__
        Mat33() {
            memset(_data_, 0, sizeof(_data_));
        }

        __host__ __device__
        Mat33(const T &n1, const T &n2, const T &n3, const T &n4, const T &n5, const T &n6, const T &n7, const T &n8,
              const T &n9) {
            _data_[0] = n1;
            _data_[1] = n2;
            _data_[2] = n3;
            _data_[3] = n4;
            _data_[4] = n5;
            _data_[5] = n6;
            _data_[6] = n7;
            _data_[7] = n8;
            _data_[8] = n9;
        }

        __host__ __device__
        Mat33(const Mat33<T> &m33) {
            memcpy_s(_data_, sizeof(_data_), m33._data_, sizeof(_data_));
        }

        __host__ __device__
        Mat33<T> &operator=(const Mat33<T> &m33) {
            memcpy_s(_data_, sizeof(_data_), m33._data_, sizeof(_data_));
            return *this;
        }

        __host__ __device__
        T &at(const size_t &row, const size_t &col) {
            return _data_[row * 3 + col];
        }

        // data()
        __host__ __device__
        T *data() {
            return _data_;
        }

        // eye()
        __host__ __device__
        static Mat33<T> eye() {
            Mat33<T> m33;
            m33._data_[0] = 1;
            m33._data_[4] = 1;
            m33._data_[8] = 1;
            return m33;
        }

        // transpose
        __host__ __device__
        Mat33<T> transpose() {
            Mat33<T> m33;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    m33._data_[i * 3 + j] = _data_[j * 3 + i];
                }
            }
            return m33;
        }

        // trace
        __host__ __device__
        T trace() {
            return _data_[0] + _data_[4] + _data_[8];
        }

        // zero
        __host__ __device__
        static Mat33<T> zero() {
            Mat33<T> m33;
            memset(m33._data_, 0, sizeof(m33._data_));
            return m33;
        }

        // +
        __host__ __device__
        inline Mat33<T> operator+(const Mat33<T> &m) {
            Mat33<T> m33;
            for (int i = 0; i < 9; i++) {
                m33._data_[i] = _data_[i] + m._data_[i];
            }
            return m33;
        }

        // -
        __host__ __device__
        inline Mat33<T> operator-(const Mat33<T> &m33) {
            Mat33<T> m;
            for (int i = 0; i < 9; i++) {
                m._data_[i] = _data_[i] - m33._data_[i];
            }
            return m;
        }

        // *
        __host__ __device__
        Mat33<T> operator*(const Mat33<T> &m33) {
            Mat33<T> m;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    m._data_[i * 3 + j] = 0;
                    for (int k = 0; k < 3; k++) {
                        m._data_[i * 3 + j] += _data_[i * 3 + k] * m33._data_[k * 3 + j];
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
                m._data_[i] = _data_[i] * value;
            }
            return m;
        }

        // *
        __host__ __device__
        inline Mat33<T> operator*(const Vec3<T> &v) {
            Mat33<T> m;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    m._data_[i * 3 + j] = 0;
                    for (int k = 0; k < 3; k++) {
                        m._data_[i * 3 + j] += _data_[i * 3 + k] * v._data_[k];
                    }
                }
            }
            return m;
        }

        // *=
        __host__ __device__
        inline Mat33<T> &operator*=(const T &value) {
            for (int i = 0; i < 9; i++) {
                _data_[i] *= value;
            }
            return *this;
        }

        // *=
        __host__ __device__
        inline Mat33<T> &operator*=(const Mat33<T> &m33) {
            Mat33<T> m;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    m._data_[i * 3 + j] = 0;
                    for (int k = 0; k < 3; k++) {
                        m._data_[i * 3 + j] += _data_[i * 3 + k] * m33._data_[k * 3 + j];
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
                    m._data_[i * 3 + j] = 0;
                    for (int k = 0; k < 3; k++) {
                        m._data_[i * 3 + j] += _data_[i * 3 + k] * v._data_[k];
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
                m._data_[i] = _data_[i] / value;
            }
            return m;
        }

        // +=
        __host__ __device__
        inline Mat33<T> &operator+=(const Mat33<T> &m) {
            for (int i = 0; i < 9; i++) {
                _data_[i] += m._data_[i];
            }
            return *this;
        }

        // -=
        __host__ __device__
        inline Mat33<T> &operator-=(const Mat33<T> &m) {
            for (int i = 0; i < 9; i++) {
                _data_[i] -= m._data_[i];
            }
            return *this;
        }

        // /=
        __host__ __device__
        inline Mat33<T> &operator/=(const T &value) {
            for (int i = 0; i < 9; i++) {
                _data_[i] /= value;
            }
            return *this;
        }

        // cout
        __host__ __device__
        friend std::ostream &operator<<(std::ostream &os, const Mat33<T> &m33) {
            os << "[" << m33._data_[0] << ", " << m33._data_[1] << ", " << m33._data_[2] << "\n";
            os << " " << m33._data_[3] << ", " << m33._data_[4] << ", " << m33._data_[5] << "\n";
            os << " " << m33._data_[6] << ", " << m33._data_[7] << ", " << m33._data_[8] << "]";
            return os;
        }
    };

    // 4x4矩阵
    template<typename T>
    struct Mat44 {

    private:
        T _data_[16];

    public:

        __host__ __device__
        const size_t size() const {
            return 16;
        }

        __host__ __device__
        Mat44() {
            memset(_data_, 0, sizeof(_data_));
        }

        __host__ __device__
        Mat44(const T& n1, const T& n2, const T& n3, const T& n4,
              const T& n5, const T& n6, const T& n7, const T& n8,
              const T& n9, const T& n10, const T& n11, const T& n12,
              const T& n13, const T& n14, const T& n15, const T& n16){
            _data_[0] = n1;
            _data_[1] = n2;
            _data_[2] = n3;
            _data_[3] = n4;
            _data_[4] = n5;
            _data_[5] = n6;
            _data_[6] = n7;
            _data_[7] = n8;
            _data_[8] = n9;
            _data_[9] = n10;
            _data_[10] = n11;
            _data_[11] = n12;
            _data_[12] = n13;
            _data_[13] = n14;
            _data_[14] = n15;
            _data_[15] = n16;
        }

        __host__ __device__
        Mat44(const Mat44<T> &m44) {
            memcpy_s(_data_, sizeof(_data_), m44._data_, sizeof(_data_));
        }

        __host__ __device__
        Mat44<T> &operator=(const Mat44<T> &m44) {
            memcpy_s(_data_, sizeof(_data_), m44._data_, sizeof(_data_));
            return *this;
        }

        // at(row, col)
        __host__ __device__
        T &at(const size_t &row, const size_t &col) {
            return _data_[row * 4 + col];
        }

        // data()
        __host__ __device__
        T *data() {
            return _data_;
        }

        // eye()
        __host__ __device__
        static Mat44<T> eye() {
            Mat44<T> m44;
            m44._data_[0] = 1;
            m44._data_[5] = 1;
            m44._data_[10] = 1;
            m44._data_[15] = 1;
            return m44;
        }

        __host__ __device__
        static Mat44<T> zero() {
            Mat44<T> m44;
            memset(m44._data_, 0, sizeof(m44._data_));
            return m44;
        }

        // transpose
        __host__ __device__
        Mat44<T> transpose() {
            Mat44<T> m;
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    m._data_[i * 4 + j] = _data_[j * 4 + i];
                }
            }
            return m;
        }

        // trace
        __host__ __device__
        T trace() {
            return _data_[0] + _data_[5] + _data_[10] + _data_[15];
        }

        // +
        __host__ __device__
        Mat44<T> operator+(const Mat44<T> &m44) {
            Mat44<T> m;
            for (int i = 0; i < 16; i++) {
                m._data_[i] = _data_[i] + m44._data_[i];
            }
            return m;
        }

        // -
        __host__ __device__
        Mat44<T> operator-(const Mat44<T> &m44) {
            Mat44<T> m;
            for (int i = 0; i < 16; i++) {
                m._data_[i] = _data_[i] - m44._data_[i];
            }
            return m;
        }

        // *
        __host__ __device__
        Mat44<T> operator*(const Mat44<T> &m44) {
            Mat44<T> m;
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    m._data_[i * 4 + j] = 0;
                    for (int k = 0; k < 4; k++) {
                        m._data_[i * 4 + j] += _data_[i * 4 + k] * m44._data_[k * 4 + j];
                    }
                }
            }
        }

        __host__ __device__
        Mat44<T> operator*(const T &value) {
            Mat44<T> m;
            for (int i = 0; i < 16; i++) {
                m._data_[i] = _data_[i] * value;
            }
            return m;
        }

        // *
        __host__ __device__
        Mat44<T> operator*(const Vec4<T> &v) {
            Mat44<T> m;
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    m._data_[i * 4 + j] = 0;
                    for (int k = 0; k < 4; k++) {
                        m._data_[i * 4 + j] += _data_[i * 4 + k] * v._data_[k];
                    }
                }
            }
            return m;
        }

        // *=
        __host__ __device__
        Mat44<T> &operator*=(const T &value) {
            for (int i = 0; i < 16; i++) {
                _data_[i] *= value;
            }
            return *this;
        }

        // *=
        __host__ __device__
        Mat44<T> &operator*=(const Mat44<T> &m44) {
            Mat44<T> m;
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    m._data_[i * 4 + j] = 0;
                    for (int k = 0; k < 4; k++) {
                        m._data_[i * 4 + j] += _data_[i * 4 + k] * m44._data_[k * 4 + j];
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
                    m._data_[i * 4 + j] = 0;
                    for (int k = 0; k < 4; k++) {
                        m._data_[i * 4 + j] += _data_[i * 4 + k] * v._data_[k];
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
                m._data_[i] = _data_[i] / value;
            }
            return m;
        }

        // /=
        __host__ __device__
        Mat44<T> &operator/=(const T &value) {
            for (int i = 0; i < 16; i++) {
                _data_[i] /= value;
            }
            return *this;
        }

        // +=
        __host__ __device__
        Mat44<T> &operator+=(const Mat44<T> &m44) {
            for (int i = 0; i < 16; i++) {
                _data_[i] += m44._data_[i];
            }
            return *this;
        }

        // -=
        __host__ __device__
        Mat44<T> &operator-=(const Mat44<T> &m44) {
            for (int i = 0; i < 16; i++) {
                _data_[i] -= m44._data_[i];
            }
            return *this;
        }

        // cout
        __host__ __device__
        friend std::ostream &operator<<(std::ostream &os, const Mat44<T> &m44) {
            os << "[" << m44._data_[0] << ", " << m44._data_[1] << ", " << m44._data_[2] << ", " << m44._data_[3] << "\n";
            os << " " << m44._data_[4] << ", " << m44._data_[5] << ", " << m44._data_[6] << "," << m44._data_[7] << "\n";
            os << " " << m44._data_[8] << ", " << m44._data_[9] << ", " << m44._data_[10] << "," << m44._data_[11] << "\n";
            os << " " << m44._data_[12] << ", " << m44._data_[13] << ", " << m44._data_[14] << "," << m44._data_[15] << "]";
            return os;
        }
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
