#ifndef _SQUARE_CONE_H_
#define _SQUARE_CONE_H_

// Square light-collector ray tracing model.
// 原始光收集器算法作者：刘伟。
// Integrated here as an optional camera collector model for LACT_sim.
#include "LightCollectorMathUtils.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <functional>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

namespace Cone {
    constexpr double kcutoff = 1e-5;

    using DirectionVecter = MathUtils::D3Vecter;
    using Position = MathUtils::D3Vecter;
    using MathUtils::magtitude;
    using MathUtils::quadratic_equation_root;

    inline double get_distance(const Position &p1, const Position &p2) {
        return (p1 - p2).get_magtitude();
    }

    class Surface {
    public:
        Surface() = default;
        virtual ~Surface() = default;
        virtual DirectionVecter get_normal_vector(const Position &p) const = 0;
        virtual std::pair<Position, bool>
        get_intersect_position(const Position &p,
                               const DirectionVecter &d) const = 0;
        // 根据反射角等于入射角，可以通过向量的计算得到反射光的方向：
        // R=I-2*(I \dot N)N
        auto get_reflect_direction(const Position &p, DirectionVecter &v) {
            auto ray_direction = v;
            DirectionVecter result;
            ray_direction.normalize_vector();
            auto normal_vector = get_normal_vector(p);
            auto product = ray_direction.dot_product(normal_vector);
            auto reflect = (ray_direction - (2 * product * normal_vector));
            return reflect.normalize_vector();
        }

    protected:
        Surface(const Surface &s) = default;
        Surface(Surface &&s) = default;
        Surface &operator=(const Surface &s) = default;
        Surface &operator=(Surface &&s) = default;
    };

    class Material {
    public:
        Material() = default;
        virtual ~Material() = default;
        virtual bool is_reflective(double theta) = 0;
        virtual DirectionVecter
        get_reflect_direction(const Position &pos,
                              const DirectionVecter &v) = 0;
        virtual double get_reflect_intensity(double theta) = 0;
    };

    class ParabolicCylindricalSurface : public Surface {
    public:
        ParabolicCylindricalSurface() = delete;
        ParabolicCylindricalSurface(const MathUtils::DMatrix<3, 3> &transform)
            : transform_(std::make_unique<MathUtils::DMatrix<3, 3>>(
                  std::move(transform))) {
            inverse_transform_ = std::make_unique<MathUtils::DMatrix<3, 3>>(
                std::move(transform_->inverse()));
        }

        std::pair<Position, bool>
        get_intersect_position(const Position &pos,
                               const DirectionVecter &vec) const override {
            auto p = inverse_transform_->product(pos);
            auto v = inverse_transform_->product(vec);
            v.normalize_vector();
            // 平行于柱面方向的光线，不会与柱面相交
            if ((fabs(v.x_) < 1e-6) && (fabs(v.z_) < 1e-6)) {
                return std::make_pair(pos, false);
            }
            bool is_intersect = false;
            Position result;

            auto [delta, roots] = solve_equation(p.x_, v.x_, p.z_, v.z_);
            std::vector<std::pair<double, Position>> positions;
            if (delta >= 0) {
                for (auto root : roots) {
                    auto theta = acos(root);
                    // 限制 theta 的范围
                    if (theta <= theta_min_ || theta >= theta_max_) {
                        continue;
                    }
                    auto x = para_a_ * sin(theta - para_b_) / (1 - cos(theta)) -
                             para_c_;
                    auto z = para_a_ * cos(theta - para_b_) / (1 - cos(theta));
                    auto t = (fabs(v.z_) > 1e-6) ? (z - p.z_) / v.z_
                                                 : (x - p.x_) / v.x_;
                    auto y = p.y_ + t * v.y_;
                    auto distance = get_distance(p, {x, y, z});
                    // 因为是变换到 x 正半轴，所以 x 必须大于
                    // 0，由于相交和向量的方向 t 必须大于
                    // 0，最后距离判断是保证解出来的是对的
                    if ((x < 0) || (t <= 1e-6) ||
                        (fabs(distance - fabs(t)) > 1e-4)) {
                        continue;
                    }
                    positions.push_back(
                        std::make_pair(distance, Position(x, y, z)));
                }
            }
            if (!positions.empty()) {
                auto result_ptr = std::min_element(
                    positions.begin(), positions.end(),
                    [](const std::pair<double, Position> &p1,
                       const std::pair<double, Position> &p2) {
                        return std::fabs(p1.first) < std::fabs(p2.first);
                    });
                result = result_ptr->second;
                is_intersect = (result.x_ > 0 && result.z_ > 0);
            }
            return std::make_pair(result.pre_multiply(*transform_),
                                  is_intersect);
        }

        // 求解抛物柱面的法向量。
        // x=a*sin(theta-b)/(1-cos(theta))+c;z=a*cos(theta-b)/(1-cos(theta))。
        // 对 t 和 自由变量 y 求导，得到两个切向量，再求叉积即可得到法向量。
        DirectionVecter get_normal_vector(const Position &pos) const override {
            auto p = inverse_transform_->product(pos);
            auto theta = (fabs(p.z_) > 1e-6)
                             ? atan(((p.x_ + para_c_) / p.z_)) + para_b_
                             : M_PI / 2 + para_b_;
            auto sin_theta = sin(theta);
            auto cos_theta = cos(theta);
            auto x = para_a_ * cos(para_b_ - theta) / (1 - cos_theta) +
                     para_a_ * sin(para_b_ - theta) * sin_theta /
                         ((1 - cos_theta) * (1 - cos_theta));
            auto z = para_a_ * sin(para_b_ - theta) / (1 - cos_theta) -
                     para_a_ * cos(para_b_ - theta) * sin_theta /
                         ((1 - cos_theta) * (1 - cos_theta));
            DirectionVecter vector_1{0.0, 1.0, 0.0};
            DirectionVecter vector_2{x, 0, z};
            auto n = vector_1.cross_product(vector_2)
                         .normalize_vector()
                         .pre_multiply(*transform_);
            return n;
        }

    private:
        // 求解抛物柱面和光线的交点。抛物柱面的方程为一般写为 z=ax^2（其中 z
        // 和 x 都是坐标值）而这里用的是极坐标下的参数方程的一个变形：
        // x=a*sin(t-b)/(1-cos(t))+c;z=a*cos(t-b)/(1-cos(t))。
        // 光线的方程是直线的参数方程即 x=p+t*v（其中 x 和 p
        // 是坐标向量，v 是光线的单位方向向量，这里的 t
        // 是另一个参数）。将后者代入前者，解方程即得参量 t
        // 的值，最后代入光线方程即可得到交点的坐标。
        // 这个函数即是解联立的方程组，
        // 返回值是抛物线的参数方程中的参数 t（角度）。
        std::tuple<double, std::vector<double>>
        solve_equation(double pos_x, double direct_x, double pos_y,
                       double direct_y) const {
            double result = 0.0;
            auto a = pos_x * direct_y - pos_y * direct_x + para_c_ * direct_y;
            auto b =
                para_a_ * (sin(para_b_) * direct_y + cos(para_b_) * direct_x);
            auto c =
                para_a_ * (cos(para_b_) * direct_y - sin(para_b_) * direct_x);
            auto [delta, min_root, max_root] =
                quadratic_equation_root(((b - a) * (b - a) + c * c),
                                        (2 * a * (b - a)), ((a * a) - (c * c)));
            std::vector<double> roots{min_root, max_root};
            return std::make_tuple(delta, roots);
        }

    private:
        std::unique_ptr<MathUtils::DMatrix<3, 3>> transform_;
        std::unique_ptr<MathUtils::DMatrix<3, 3>> inverse_transform_;
        // 抛物线的参数方程中的参数范围和其他参数，参数方程：
        // x=a*sin(t-b)/(1-cos(t))+c;z=a*cos(t-b)/(1-cos(t))
        double theta_min_ = 1.3242;
        double theta_max_ = 2.2328;
        double para_a_ = 24.2216;
        double para_b_ = 0.6621;
        double para_c_ = 7.5;
    };

    class Bezier2CylindricalSurface : public Surface {
    public:
        MathUtils::DMatrix<3, 3> transform_;
        MathUtils::DMatrix<3, 3> inverse_transform_;
        // 贝塞尔曲线的矩阵表示，X=T^{tran}*B*P，其中 T^{tran} = [1,t,t^2]，B
        // 是二项式系数矩阵，P 是控制点控制向量。
        MathUtils::D3Vecter control_x_{0, 0, 0};
        MathUtils::D3Vecter control_z_{0, 0, 0};
        MathUtils::DMatrix<3, 3> bezier_matrix_{1, 0, 0, -2, 2, 0, 1, -2, 1};

        Bezier2CylindricalSurface() = delete;
        Bezier2CylindricalSurface(const MathUtils::DMatrix<3, 3> &transform,
                                  const MathUtils::D3Vecter &control_x,
                                  const MathUtils::D3Vecter &control_z)
            : transform_(transform), control_x_(control_x),
              control_z_(control_z) {
            inverse_transform_ = (std::move(transform_.inverse()));
        }

        std::pair<Position, bool>
        get_intersect_position(const Position &pos,
                               const DirectionVecter &vec) const override {
            auto p = pos.pre_multiply(inverse_transform_);
            auto v = vec.pre_multiply(inverse_transform_).normalize_vector();
            // 平行于柱面方向的光线，不会与柱面相交
            if ((fabs(v.x_) < 1e-6) && (fabs(v.z_) < 1e-6)) {
                return std::make_pair(pos, false);
            }
            Position result;
            bool is_intersect = false;

            auto x = control_x_.pre_multiply(bezier_matrix_);
            auto z = control_z_.pre_multiply(bezier_matrix_);
            auto a = x[2] * v.z_ - z[2] * v.x_;
            auto b = x[1] * v.z_ - z[1] * v.x_;
            auto c = x[0] * v.z_ - z[0] * v.x_ - p.x_ * v.z_ + p.z_ * v.x_;
            // a == 0 的时候返回的大小根都是一样的，delta == 1.0
            auto [delta, min_root, max_root] = quadratic_equation_root(a, b, c);
            std::vector<std::pair<double, Position>> positions;
            if (delta > 0) {
                for (auto root : {min_root, max_root}) {
                    // t 的范围一定是 0~1
                    if (root < 0 || root > 1) {
                        continue;
                    }
                    MathUtils::D3Vecter t_vector(1, root, root * root);
                    auto pos_x = x.dot_product(t_vector);
                    auto pos_z = z.dot_product(t_vector);
                    auto s = (fabs(v.z_) > 1e-6) ? (pos_z - p.z_) / v.z_
                                                 : (pos_x - p.x_) / v.x_;
                    auto pos_y = p.y_ + s * v.y_;
                    auto distance = get_distance(p, {pos_x, pos_y, pos_z});

                    if ((pos_x < 0) || (s <= 1e-6) ||
                        (fabs(distance - fabs(s)) > 1e-4)) {
                        continue;
                    }
                    positions.emplace_back(std::make_pair(
                        distance, Position(pos_x, pos_y, pos_z)));
                }
            }
            if (!positions.empty()) {
                auto result_ptr = std::min_element(
                    positions.begin(), positions.end(),
                    [](const std::pair<double, Position> &p1,
                       const std::pair<double, Position> &p2) {
                        return std::fabs(p1.first) < std::fabs(p2.first);
                    });
                result = result_ptr->second;
                is_intersect = (result.x_ > 0 && result.z_ > 0 &&
                                (fabs(result.y_) < result.x_));
            }
            return std::make_pair(result.pre_multiply(transform_),
                                  is_intersect);
        }

        DirectionVecter get_normal_vector(const Position &pos) const override {
            auto p = pos.pre_multiply(inverse_transform_);
            DirectionVecter result;
            auto x = bezier_matrix_.product(control_x_);
            auto z = bezier_matrix_.product(control_z_);

            auto [deltaz, min_rootz, max_rootz] =
                quadratic_equation_root(z[2], z[1], z[0] - p.z_);
            auto [deltax, min_rootx, max_rootx] =
                quadratic_equation_root(x[2], x[1], x[0] - p.x_);
            std::vector<double> roots{min_rootz, max_rootz, min_rootx,
                                      max_rootx};
            // 有两个根相等，说明是切线对应的那个点
            auto it =
                std::find_if(roots.begin(), roots.end(), [&roots](double num) {
                    return std::count_if(
                               roots.begin(), roots.end(), [num](double other) {
                                   return std::fabs(num - other) < 1e-5;
                               }) > 1;
                });

            auto t = *it;
            MathUtils::D3Vecter t_vector(0, 1, t * 2);
            auto dx = x.dot_product(t_vector);
            auto dz = z.dot_product(t_vector);
            DirectionVecter vector_1{0.0, 1.0, 0.0};
            DirectionVecter vector_2{dx, 0, dz};
            result = vector_1.cross_product(vector_2)
                         .normalize_vector()
                         .pre_multiply(transform_);
            return result;
        }
    };

    class Plane : public Surface {
    public:
        Plane() = default;
        Plane(const DirectionVecter &normal_vector)
            : normal_vector_(normal_vector) {
            center_ = Position{0, 0, 0};
        }
        Plane(const DirectionVecter &x, const DirectionVecter &y,
              const Position &p)
            : center_(p), x_axis_(x), y_axis_(y),
              normal_vector_(x.cross_product(y)) {}
        Plane(std::initializer_list<double> &&lst) : normal_vector_(lst) {}
        Plane(const DirectionVecter &normal_vector, const Position &center)
            : normal_vector_(normal_vector), center_(center) {}
        Plane(double width, double length) : width_(width), length_(length) {}
        Plane(double width, double length, const Position &center,
              const DirectionVecter &x, const DirectionVecter &y)
            : width_(width), length_(length), center_(center), x_axis_(x),
              y_axis_(y),
              normal_vector_(x.cross_product(y).normalize_vector()) {}

        ~Plane() = default;

        Plane &set_center(const Position &p) {
            center_ = p;
            return *this;
        }

        DirectionVecter get_normal_vector(const Position &p) const override {
            return normal_vector_;
        }
        std::pair<Position, bool>
        get_intersect_position(const Position &p,
                               const DirectionVecter &v) const override {
            // 分数 0/0 会导致无穷大，这里处理一下
            auto up = normal_vector_.dot_product(p - center_);
            auto down = normal_vector_.dot_product(v);
            double t = 0;
            bool is_intersect = false;
            Position result;
            if (!(std::fabs(down) < 1e-6)) {
                t = -(up / down);
                result = p + (v * t);
                auto dv = (result - center_);
                auto x = dv.dot_product(x_axis_) / x_axis_.get_magtitude();
                auto y = dv.dot_product(y_axis_) / y_axis_.get_magtitude();
                // auto distance = get_distance(result, center_);
                is_intersect =
                    ((fabs(x) < length_ / 2) && (fabs(y) < width_ / 2));
            } else {
                auto x = (p - center_).dot_product(x_axis_) /
                         x_axis_.get_magtitude();
                auto y = (p - center_).dot_product(y_axis_) /
                         y_axis_.get_magtitude();
                is_intersect =
                    ((fabs(x) < length_ / 2) && (fabs(y) < width_ / 2));
            }
            return std::make_pair(result, is_intersect);
        }

        std::pair<bool, Position>
        is_in_plane(const Position &p, const DirectionVecter &v,
                    std::function<bool(double, double)> fn = nullptr) const {
            auto [pos, is_intersect] = get_intersect_position(p, v);
            auto dv = (pos - center_);
            auto x = dv.dot_product(x_axis_) / x_axis_.get_magtitude();
            auto y = dv.dot_product(y_axis_) / y_axis_.get_magtitude();
            auto is_out = ((is_intersect) && (fabs(pos.x_) < length_ / 2) &&
                           (fabs(pos.y_) < width_ / 2));

            if (fn == nullptr) {
                is_out = is_intersect && (fabs(pos.x_) < length_ / 2) &&
                         (fabs(pos.y_) < width_ / 2);
            } else {
                is_out = is_intersect && fn(x, y);
            }
            return std::make_pair(is_out, pos);
        }
        DirectionVecter normal_vector_{0, 0, 1};
        DirectionVecter x_axis_{1, 0, 0};
        DirectionVecter y_axis_{0, 1, 0};
        Position center_{0, 0, 0};
        // x 轴定为长度方向，y 轴定为宽度方向
        double length_{15.0};
        double width_{15.0};
    };

    class MirrorReflectMaterial : public Material {
    public:
        MirrorReflectMaterial() = default;
        ~MirrorReflectMaterial() = default;
        inline bool is_reflective(double theta) override { return true; }

        inline double get_reflect_intensity(double theta) override {
            return 1.0;
        }

        inline DirectionVecter
        get_reflect_direction(const Position &pos,
                              const DirectionVecter &vec) override {
            return vec;
        }
    };

    class LambertianReflectMaterial : public Material {
    public:
        LambertianReflectMaterial() = default;
        ~LambertianReflectMaterial() = default;
        inline bool is_reflective(double theta) override { return true; }

        inline double get_reflect_intensity(double theta) override {
            return fabs(cos(theta));
        }

        inline DirectionVecter
        get_reflect_direction(const Position &pos,
                              const DirectionVecter &vec) override {
            return vec;
        }
    };

    class TrueReflectMaterial : public Material {
    public:
        TrueReflectMaterial() = default;
        ~TrueReflectMaterial() = default;
        inline bool is_reflective(double theta) override { return true; }

        inline double get_reflect_intensity(double theta) override {
            double a1 = 2.81219e-5;
            double f0 = 0.90565;
            double t1 = -11.0723;
            auto reflectivity = a1 * exp(-(((theta)) / t1)) + f0;
            return reflectivity;
        }

        inline DirectionVecter
        get_reflect_direction(const Position &pos,
                              const DirectionVecter &vec) override {
            return vec;
        }
    };

    class TableReflectMaterial : public Material {
    public:
        explicit TableReflectMaterial(std::vector<std::pair<double, double>> points)
            : points_(std::move(points)) {
            std::sort(points_.begin(), points_.end(),
                      [](const auto &a, const auto &b) { return a.first < b.first; });
        }

        inline bool is_reflective(double theta) override {
            return get_reflect_intensity(theta) > 0.0;
        }

        inline double get_reflect_intensity(double theta) override {
            if (points_.empty()) return 1.0;
            if (theta <= points_.front().first) return points_.front().second;
            if (theta >= points_.back().first) return points_.back().second;
            auto it = std::lower_bound(
                points_.begin(), points_.end(), theta,
                [](const auto &point, double value) { return point.first < value; });
            const auto &hi = *it;
            const auto &lo = *(it - 1);
            const double t = (theta - lo.first) / (hi.first - lo.first);
            return lo.second + t * (hi.second - lo.second);
        }

        inline DirectionVecter
        get_reflect_direction(const Position &pos,
                              const DirectionVecter &vec) override {
            return vec;
        }

    private:
        std::vector<std::pair<double, double>> points_;
    };

    class SquareCone {
    public:
        SquareCone() {
            // 两个平面都是默认法向量为 z 轴正方向
            // 抛物面是算出来这个高度的具体位置，所以不加 start
            out_plane_ = std::make_unique<Plane>(cone_out_, cone_out_);
            in_plane_ = std::make_unique<Plane>(cone_in_, cone_in_);
            in_plane_->set_center({0, 0, height_});
            out_plane_->set_center({0, 0, height_start_});
            material_ = std::make_unique<MirrorReflectMaterial>();
        };
        SquareCone(double cone_in, double cone_out, double height,
                   double height_start)
            : cone_in_(cone_in), cone_out_(cone_out), height_(height),
              height_start_(height_start) {
            // 这个地方的相对高度是知道的，所以加上一开始的。
            out_plane_ = std::make_unique<Plane>(cone_out_, cone_out_);
            in_plane_ = std::make_unique<Plane>(cone_in_, cone_in_);
            in_plane_->set_center({0, 0, height_ + height_start_});
            out_plane_->set_center({0, 0, height_start_});
            material_ = std::make_unique<MirrorReflectMaterial>();
        }

        ~SquareCone() = default;

        SquareCone(const SquareCone &s) = delete;
        SquareCone(SquareCone &&s) = delete;
        SquareCone &operator=(const SquareCone &s) = delete;
        SquareCone &operator=(SquareCone &&s) = delete;

        SquareCone &set_material(std::unique_ptr<Material> &&material) {
            material_ = std::move(material);
            return *this;
        }

        SquareCone &add_surface(std::unique_ptr<Surface> &&surface) {
            surfaces_.emplace_back(std::move(surface));
            return *this;
        }

        std::pair<std::tuple<size_t, double, Position>, bool>
        get_intersect_position(const Position &pos,
                               DirectionVecter &vec) const {
            auto v = vec;
            v.normalize_vector();
            std::vector<std::tuple<size_t, double, Position>> positions;
            bool is_intersect = false;
            std::tuple<size_t, double, Position> result;

            for (size_t index = 0; index < surfaces_.size(); index++) {
                auto [res, _is_intersect] =
                    surfaces_[index]->get_intersect_position(pos, v);
                auto distance = get_distance(res, pos);
                if (!_is_intersect || distance < 1e-6) {
                    continue;
                }
                positions.push_back(std::make_tuple(index, distance, res));
            }

            if (!positions.empty()) {
                // 如果在对角线上需要返回固定的顺序 0->1->2->3->0 返回第二个面
                // 这里特殊处理第三个，因为第三个后面是 0
                if (positions.size() == 2 &&
                    fabs(fabs(pos.x_) - fabs(pos.y_)) < 1e-6 &&
                    std::get<0>(positions[0]) == 0 &&
                    std::get<0>(positions[1]) == 3) {
                    std::swap(positions[0], positions[1]);
                }
                auto result_ptr = std::min_element(
                    positions.begin(), positions.end(),
                    [](const std::tuple<size_t, double, Position> &p1,
                       const std::tuple<size_t, double, Position> &p2) {
                        auto d1 = std::get<1>(p1);
                        auto d2 = std::get<1>(p2);
                        return d1 < d2 || fabs(d1 - d2) < 1e-6;
                    });
                result = *result_ptr;
                is_intersect = true;
            }
            return std::make_pair(result, is_intersect);
        }

        DirectionVecter get_reflect_direction(const Position &p,
                                              DirectionVecter &v,
                                              size_t index) const {
            DirectionVecter result;
            if (fabs(fabs(p.x_) - fabs(p.y_)) < 1e-6) {
                auto i1 = index;
                auto i2 = (index == 0) ? 3 : index - 1;
                auto norm1 = surfaces_[i1]->get_normal_vector(p);
                auto norm2 = surfaces_[i2]->get_normal_vector(p);
                auto norm = (norm1 + norm2).normalize_vector();
                auto product = v.dot_product(norm);
                result = (v - 2 * product * norm).normalize_vector();
            } else {
                result = surfaces_[index]->get_reflect_direction(p, v);
            }
            return result;
        }

        std::tuple<double, bool, Position, DirectionVecter, int>
        ray_trace_impl(const Position &pos, const DirectionVecter &vec,
                       const double in_intensity) const {
            Position p = pos;
            DirectionVecter v = vec;
            v.normalize_vector();
            double intensity = in_intensity;
            size_t number = 0;
            // 检查是否直接出去
            auto [is_on, postion] = out_plane_->is_in_plane(p, v);
            bool is_out =
                (is_on && (out_plane_->normal_vector_.dot_product(v) < 0));
            p = is_out ? postion : p;

            // 检查能否在入口平面上
            if (!is_out) {
                auto [is_on_in_plane, _] = in_plane_->is_in_plane(p, v);
                bool is_out_in_plane =
                    is_on_in_plane &&
                    (in_plane_->normal_vector_.dot_product(v) > 0);

                if (!is_out_in_plane && is_on_in_plane) {
                    bool is_outside = false;
                    do {
                        auto [pos_tuple, is_intersect] =
                            get_intersect_position(p, v);
                        auto [index, distance, position] = pos_tuple;
                        if (!is_intersect) {
                            // 意外情况默认从上面离开，一般是浮点数误差造成的。
                            is_out = false;
                            break;
                        }
                        // 检查反射率
                        auto reflect_direcion =
                            get_reflect_direction(position, v, index);
                        auto theta =
                            90 -
                            57.3 * acos(reflect_direcion.dot_product(v)) / 2.0;
                        auto is_reflective = material_->is_reflective(theta);
                        if (!is_reflective) {
                            is_out = false;
                            break;
                        }
                        intensity *= material_->get_reflect_intensity(theta);
                        //  反射光线方向检查
                        auto reflect = material_->get_reflect_direction(
                            position, reflect_direcion);
                        // 检查是否在离开
                        auto [_is_on_out_plane, _p] =
                            out_plane_->is_in_plane(position, reflect);
                        auto [_is_on_in_plane, _] =
                            in_plane_->is_in_plane(position, reflect);
                        is_out_in_plane =
                            _is_on_in_plane &&
                            (in_plane_->normal_vector_.dot_product(reflect) >
                             0);
                        is_out = _is_on_out_plane &&
                                 (out_plane_->normal_vector_.dot_product(
                                      reflect) < 0);
                        is_outside = is_out || is_out_in_plane;
                        p = (is_out) ? _p : position;
                        v = reflect;
                        number++;
                    } while (!is_outside && (number <= 500));
                }
                intensity = (number == 0) ? 0.0 : intensity;
            }
            return std::make_tuple(intensity, is_out, p, v, number);
        }

        inline bool ray_trace(Position &pos, DirectionVecter &vec,
                              double &intensity) const {
            auto [i, is_out, p, v, n] = ray_trace_impl(pos, vec, intensity);
            intensity = i;
            pos = p;
            vec = v;
            return is_out;
        }

        std::vector<std::unique_ptr<Surface>> surfaces_;
        std::unique_ptr<Plane> out_plane_;
        std::unique_ptr<Plane> in_plane_;
        std::unique_ptr<Material> material_;

        double cone_in_ = 24.3992;
        double cone_out_ = 15.0014;
        double height_ = 25.2718;
        double height_start_ = 0.00144497;
    };

    inline SquareCone *create_cone_parabolic_surface() {
        auto cone = new SquareCone;
        cone->add_surface(std::make_unique<ParabolicCylindricalSurface>(
            MathUtils::DMatrix<3, 3>{1, 0, 0, 0, 1, 0, 0, 0, 1}));
        cone->add_surface(std::make_unique<ParabolicCylindricalSurface>(
            MathUtils::DMatrix<3, 3>{0, 1, 0, 1, 0, 0, 0, 0, 1}));
        cone->add_surface(std::make_unique<ParabolicCylindricalSurface>(
            MathUtils::DMatrix<3, 3>{-1, 0, 0, 0, 1, 0, 0, 0, 1}));
        cone->add_surface(std::make_unique<ParabolicCylindricalSurface>(
            MathUtils::DMatrix<3, 3>{0, 1, 0, -1, 0, 0, 0, 0, 1}));
        return cone;
    }

    inline SquareCone *create_cone_bezier_surface(double cone_in,
                                                  double cone_out,
                                                  double height,
                                                  double height_start = 0.0) {
        auto cone = new SquareCone(cone_in, cone_out, height, height_start);
        const double half_in = 0.5 * cone_in;
        const double half_out = 0.5 * cone_out;
        MathUtils::D3Vecter control_x{
            half_in,
            half_in * (11.630 / 12.2),
            half_out
        };
        MathUtils::D3Vecter control_z{
            height,
            height * (11.2841 / 29.7),
            0
        };
        cone->add_surface(std::make_unique<Bezier2CylindricalSurface>(
            MathUtils::DMatrix<3, 3>{1, 0, 0, 0, 1, 0, 0, 0, 1}, control_x,
            control_z));

        cone->add_surface(std::make_unique<Bezier2CylindricalSurface>(
            MathUtils::DMatrix<3, 3>{0, 1, 0, 1, 0, 0, 0, 0, 1}, control_x,
            control_z));

        cone->add_surface(std::make_unique<Bezier2CylindricalSurface>(
            MathUtils::DMatrix<3, 3>{-1, 0, 0, 0, 1, 0, 0, 0, 1}, control_x,
            control_z));

        cone->add_surface(std::make_unique<Bezier2CylindricalSurface>(
            MathUtils::DMatrix<3, 3>{0, 1, 0, -1, 0, 0, 0, 0, 1}, control_x,
            control_z));

        return cone;
    }

    inline SquareCone *create_cone_bezier_surface() {
        return create_cone_bezier_surface(12.2 * 2, 6.5 * 2, 29.7, 0.0);
    }

    inline SquareCone *create_cone_plane() {
        auto cone = new SquareCone(12.2 * 2, 12.2 * 2, 29.7, 0);
        cone->add_surface(std::make_unique<Plane>(
            12.2, 29.7, Position{12.2, 0, 29.7 / 2}, DirectionVecter{0, 0, 1},
            DirectionVecter{0, 1, 0}));
        cone->add_surface(std::make_unique<Plane>(
            29.7, 12.2 * 2, Position{-12.2, 0, 29.7 / 2},
            DirectionVecter{0, 1, 0}, DirectionVecter{0, 0, 1}));
        cone->add_surface(std::make_unique<Plane>(
            29.7, 12.2 * 2, Position{0, 12.2, 29.7 / 2},
            DirectionVecter{1, 0, 0}, DirectionVecter{0, 0, 1}));
        cone->add_surface(std::make_unique<Plane>(
            12.2 * 2, 29.7, Position{0, -12.2, 29.7 / 2},
            DirectionVecter{0, 0, 1}, DirectionVecter{1, 0, 0}));

        return cone;
    }

} // namespace Cone
#endif
