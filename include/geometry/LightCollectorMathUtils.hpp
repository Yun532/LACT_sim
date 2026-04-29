#ifndef _MATH_UTILS_H_
#define _MATH_UTILS_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <iostream>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace MathUtils {
    constexpr double kcutoff = 1e-5;

    template <typename... T, typename = std::enable_if_t<
                                 ((std::is_convertible_v<T, double>)&&...)>>
    inline constexpr double magtitude(T... t) {
        return sqrt(((t * t) + ...));
    }

    inline auto quadratic_equation_root(double a, double b, double c) {
        if (fabs(a) < 1e-5) {
            return std::make_tuple(1.0, -c / b, -c / b);
        }
        double delta = b * b - 4 * a * c;
        double min_root = 0.0;
        double max_root = 0.0;

        if (delta >= 0) {
            min_root = (-b - sqrt(delta)) / (2 * a);
            max_root = (-b + sqrt(delta)) / (2 * a);
        }
        return std::make_tuple(delta,
                               min_root <= max_root ? min_root : max_root,
                               max_root >= min_root ? max_root : min_root);
    }
    template <size_t n, size_t m>
    struct DMatrix;

    template <size_t n>
    DMatrix<n, n> identity_matrix();

    template <size_t n, size_t m>
    struct DMatrix {

        DMatrix() { data_.resize(n_ * m_, 0); }
        // DMatrix(const std::vector<double> &vec) { data_ = vec; }
        DMatrix(std::vector<double> &&vec) : n_(n), m_(m) {
            data_ = std::move(vec);
        }
        DMatrix(std::initializer_list<double> &&values) : data_(values) {}

        DMatrix(const DMatrix &other) {
            n_ = other.m_;
            m_ = other.m_;
            data_ = other.data_;
        }
        DMatrix &operator=(const DMatrix &other) {
            n_ = other.m_;
            m_ = other.m_;
            data_ = other.data_;
            return *this;
        }

        DMatrix(DMatrix &&other) noexcept {
            n_ = other.n_;
            m_ = other.m_;
            data_ = std::move(other.data_);
        }

        DMatrix &operator=(DMatrix &&other) {
            n_ = other.n_;
            m_ = other.m_;
            data_ = std::move(other.data_);
            return *this;
        }

        double operator[](size_t index) const {
            assert(index < data_.size());
            return data_[index];
        }

        double &operator()(size_t row, size_t column) {
            assert(row < n_ && column < m_);
            return data_[row * m_ + column];
        }
        double operator()(size_t row, size_t column) const {
            assert(row < n_ && column < m_);
            return data_[row * m_ + column];
        }

        DMatrix<n, m> &operator+(const DMatrix<n, m> &other) {
            for (size_t i = 0; i < data_.size(); i++) {
                data_[i] += other.data_[i];
            }
            return *this;
        }

        DMatrix<n, m> &operator-(const DMatrix<n, m> &other) {
            for (size_t i = 0; i < data_.size(); i++) {
                data_[i] -= other.data_[i];
            }
            return *this;
        }

        template <size_t _m>
        auto product(const DMatrix<m, _m> &matrix) const  {
            DMatrix<n, _m> result;
            for (size_t i = 0; i < n_; i++) {
                for (size_t j = 0; j < _m; j++) {
                    for (size_t k = 0; k < m_; k++) {
                        result(i, j) += (*this)(i, k) * matrix(k, j);
                    }
                }
            }
            return result;
        }

        DMatrix<n, m> &set_data_row(size_t row, const DMatrix<1, m> &value) {
            for (size_t i = 0; i < m_; i++) {
                data_[row * m_ + i] = value(1, i);
            }
            return *this;
        }

        DMatrix<n, m> &set_data_column(size_t column,
                                       const DMatrix<n, 1> &value) {
            for (size_t i = 0; i < n_; i++) {
                data_[i * m_ + column] = value(i, 1);
            }
            return *this;
        }

        DMatrix<n, m> &swap_row(size_t r1, size_t r2) {
            assert(r1 < n_ && r2 < n_);
            for (size_t i = 0; i < m_; i++) {
                std::swap(data_[r1 * m_ + i], data_[r2 * m_ + i]);
            }
            return *this;
        }

        DMatrix<1, m> get_data_row(size_t row) const {
            std::vector<double> result(m_, 0);
            if (row < n_) {
                for (size_t i = 0; i < m_; i++) {
                    result.at(i) = this->data_[row * m_ + i];
                }
            }
            return DMatrix<1, m>(std::move(result));
        }

        DMatrix<n, 1> get_data_column(size_t column) const {
            std::vector<double> result(n_, 0);
            if (column < m_) {
                for (size_t i = 0; i < n_; i++) {
                    result.at(i) = data_[i * m_ + column];
                }
            }
            return DMatrix<n, 1>(std::move(result));
        }

        DMatrix<m, n> transpose() const {
            DMatrix<m, n> result;
            for (size_t i = 0; i < n_; i++) {
                for (size_t j = 0; j < m_; j++) {
                    result(j, i) = data_[i * m_ + j];
                }
            }
            return result;
        }

        // 先将最大的元素挑出来作为每行的主元，同时调整置换矩阵，然后作 LU
        // 分解
        auto plu_decompose() const -> std::tuple<DMatrix<n, n>, DMatrix<n, n>,
                                                 DMatrix<n, n>, bool, size_t> {
            static_assert((n == m), "must be square matrix");
            DMatrix up_triangle_matrix(*this);
            auto down_triangle_matrix = identity_matrix<n>();
            auto permutation_matrix = identity_matrix<n>();
            bool is_singular = false;
            size_t number = 0;
            for (size_t k = 0; k < n_; k++) {
                auto max_pivot = std::fabs(up_triangle_matrix(k, k));
                auto max_row = k;
                for (size_t i = k + 1; i < n_; i++) {
                    auto element = std::fabs(up_triangle_matrix(i, k));
                    if (element > max_pivot) {
                        max_pivot = element;
                        max_row = i;
                    }
                }
                if (max_row != k) {
                    up_triangle_matrix.swap_row(k, max_row);
                    permutation_matrix.swap_row(k, max_row);
                    number++;
                    for (size_t j = 0; j < k; j++) {
                        std::swap(down_triangle_matrix(k, j),
                                  down_triangle_matrix(max_row, j));
                    }
                }
                for (size_t i = k + 1; i < n_; i++) {
                    // 如果从主元开始往下的所有元素都是
                    // 0，那么这个矩阵是奇异的（前面已经将最大的元素挑出来作为主元）
                    if (std::fabs(up_triangle_matrix(k, k)) < kcutoff) {
                        up_triangle_matrix(k, k) = 0.0;
                        k = n_;
                        is_singular = true;
                        break;
                    }
                    down_triangle_matrix(i, k) =
                        up_triangle_matrix(i, k) / up_triangle_matrix(k, k);
                    for (size_t j = k; j < n_; j++) {
                        up_triangle_matrix(i, j) -= down_triangle_matrix(i, k) *
                                                    up_triangle_matrix(k, j);
                    }
                }
            }
            // 最后一行如果全是零，也说明是奇异矩阵
            auto index = n_ - 1;
            size_t equal_zero = 0;
            for (size_t i = 0; i < n_; i++) {
                if (std::fabs(up_triangle_matrix(index, i)) < kcutoff) {
                    equal_zero++;
                }
            }
            if (equal_zero == n_) {
                is_singular = true;
            }
            return std::make_tuple(permutation_matrix, up_triangle_matrix,
                                   down_triangle_matrix, is_singular, number);
        }

        // 行列式即 LU 分解得到上三角矩阵，然后对角线元素相乘
        inline double determinant() {
            if (!(det_ < 0.0)) {
                return det_;
            }
            static_assert((n == m), "must be square matrix");
            double result = 0.0;
            auto [permutation_matrix, up_triangle_matrix, down_triangle_matrix,
                  is_singular, number] = this->plu_decompose();
            if (!is_singular) {
                result = (number % 2 == 1) ? -1 : 1;
                for (size_t i = 0; i < n_; i++) {
                    result *= up_triangle_matrix(i, i);
                }
            }
            det_ = result;
            return result;
        }

        // 前向替换法，解方程组（适用于下三角矩阵）
        template <size_t _m>
        DMatrix<n, _m> forward_substitute(const DMatrix<n, _m> &b) {
            static_assert((n == m), "triangle matrix must be square");
            DMatrix<n, _m> result(b);
            for (size_t i = 0; i < n_; i++) {
                // 遍历 b 的所有列向量，backwards substitute 相同
                for (size_t k = 0; k < _m; k++) {
                    for (size_t j = 0; j < i; j++) {
                        result(i, k) -= result(j, k) * (*this)(i, j);
                    }
                    result(i, k) /= (*this)(i, i);
                }
            }
            return result;
        }

        // 后向替换法，解方程组（适用于上三角矩阵）
        template <size_t _m>
        DMatrix<n, _m> backward_substitute(const DMatrix<n, _m> &b) {
            static_assert((n == m), "triangle matrix must be square");
            DMatrix<n, _m> result(b);
            for (int i = n_ - 1; i >= 0; i--) {
                for (size_t k = 0; k < _m; k++) {
                    for (size_t j = i + 1; j < n_; j++) {
                        result(i, k) -= result(j, k) * (*this)(i, j);
                    }
                    result(i, k) /= (*this)(i, i);
                }
            }
            return result;
        }

        DMatrix<n, m> inverse() {
            auto [permutation_matrix, up_triangle_matrix, down_triangle_matrix,
                  is_singular, number] = this->plu_decompose();
            DMatrix<n, m> result;
            if (!is_singular) {
                auto y =
                    down_triangle_matrix.forward_substitute(permutation_matrix);
                result = std::move(up_triangle_matrix.backward_substitute(y));
            }
            return result;
        }

        void print() {
            std::cout << "\n";
            for (size_t i = 0; i < n_; i++) {
                std::cout << "[ ";
                for (size_t j = 0; j < m_; j++) {
                    std::cout << data_[i * m_ + j] << " ";
                }
                std::cout << "]"
                          << "\n";
            }
            std::cout << std::endl;
        }

        inline size_t get_size() { return data_.size(); }

        inline size_t get_row_size() { return n_; }

        inline size_t get_column_size() { return m_; }

    private:
        // 行和列的 index 都是从 0 开始
        size_t n_ = n;
        size_t m_ = m;
        double det_ = -1.0;
        std::vector<double> data_;
    };

    template <size_t n>
    DMatrix<n, n> identity_matrix() {
        DMatrix<n, n> result;
        for (size_t i = 0; i < n; i++) {
            result(i, i) = 1;
        }
        return result;
    }

    using D3Vecter = DMatrix<3, 1>;
    template <>
    class DMatrix<3, 1> {
    private:
        std::vector<double> data_ = {0.0, 0.0, 0.0};
        size_t n_ = 3;
        size_t m_ = 1;

    public:
        std::reference_wrapper<double> x_ = data_[0];
        std::reference_wrapper<double> y_ = data_[1];
        std::reference_wrapper<double> z_ = data_[2];

        DMatrix<3, 1>() = default;
        DMatrix<3, 1>(double x, double y, double z) : data_{x, y, z} {}
        DMatrix<3, 1>(const std::vector<double> &data) {
            data_.resize(3, 0);
            for (size_t i = 0; i < 3; i++) {
                data_[i] = data[i];
            }
        }
        DMatrix<3, 1>(const DMatrix<3, 1> &v)
            : data_(v.data_), x_(data_[0]), y_(data_[1]), z_(data_[2]) {}
        DMatrix<3, 1>(DMatrix<3, 1> &&v) noexcept
            : data_(std::move(v.data_)), x_(data_[0]), y_(data_[1]),
              z_(data_[2]) {}

        DMatrix<3, 1> &operator=(const DMatrix<3, 1> &v) {
            data_ = v.data_;
            x_ = data_[0];
            y_ = data_[1];
            z_ = data_[2];
            return *this;
        }
        DMatrix<3, 1> &operator=(DMatrix<3, 1> &&v) {
            data_ = std::move(v.data_);
            x_ = data_[0];
            y_ = data_[1];
            z_ = data_[2];
            return *this;
        }

        double &operator[](size_t n) {
            assert(n < 3);
            return data_[n];
        }
        double operator[](size_t n) const {
            assert(n < 3);
            return data_[n];
        }
        double &operator()(size_t n, size_t _) { return (*this)[n]; }
        double operator()(size_t n, size_t _) const { return (*this)[n]; }

        D3Vecter operator+(const D3Vecter &v) const {
            return D3Vecter(x_ + v.x_, y_ + v.y_, z_ + v.z_);
        }

        D3Vecter operator-(const D3Vecter &v) const {
            return D3Vecter(x_ - v.x_, y_ - v.y_, z_ - v.z_);
        }

        D3Vecter operator*(double scalar) const {
            return D3Vecter(x_ * scalar, y_ * scalar, z_ * scalar);
        }
        D3Vecter operator/(double scalar) const {
            return D3Vecter(x_ / scalar, y_ / scalar, z_ / scalar);
        }

        friend D3Vecter operator*(double scalar, const D3Vecter &v) {
            return D3Vecter(v.x_ * scalar, v.y_ * scalar, v.z_ * scalar);
        }

        D3Vecter &normalize_vector() {
            double norm = magtitude(x_, y_, z_);
            if (!((std::fabs(norm - 1.0)) < kcutoff)) {
                for (size_t i = 0; i < 3; i++) {
                    data_[i] /= norm;
                }
            }
            return *this;
        };

        D3Vecter &change_direction() {
            for (size_t i = 0; i < 3; i++) {
                data_[i] = -data_[i];
            }
            return *this;
        }

        inline double get_magtitude() const {
            return magtitude(x_, y_, z_);
        }

        inline double dot_product(const D3Vecter &v) const {
            return x_ * v.x_ + y_ * v.y_ + z_ * v.z_;
        }

        inline D3Vecter cross_product(const D3Vecter &v) const {
            return D3Vecter(y_ * v.z_ - z_ * v.y_, z_ * v.x_ - x_ * v.z_,
                            x_ * v.y_ - y_ * v.x_);
        }

        // 在向量的前面乘一个矩阵
        template <size_t _n>
        inline D3Vecter pre_multiply(const DMatrix<_n, 3> &m) const {
            D3Vecter result;
            for (size_t i = 0; i < _n; i++) {
                for (size_t j = 0; j < 3; j++) {
                    result[i] += m(i, j) * data_[j];
                }
            }
            return result;
        }

        // 在向量的后面乘一个矩阵
        template <size_t _m>
        inline DMatrix<3, _m> post_multiply(const DMatrix<1, _m> &m) const{
            DMatrix<3, _m> result;
            for (size_t i = 0; i < 3; i++) {
                for (size_t j = 0; j < _m; j++) {
                    result(i, j) = data_[i] * m[j];
                }
            }
            return result;
        }

        template <size_t _m>
        inline DMatrix<3, _m> post_multiply(const DMatrix<_m, 1> &m) const {
            DMatrix<3, _m> result;
            for (size_t i = 0; i < 3; i++) {
                for (size_t j = 0; j < _m; j++) {
                    result(i, j) = data_[i] * m[j];
                }
            }
            return result;
        }

        inline void print() const {
            std::cout << "\n[ " << x_ << " " << y_ << " " << z_ << " ]"
                      << std::endl;
        }

        inline bool is_normal_value(){
            for(size_t i = 0;i<3;i++){
                if (std::isnan(data_[i])){
                    return false;
                }
            }
            return true;
        }
    };

} // namespace MathUtils

#endif