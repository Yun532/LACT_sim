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

    struct SurfaceIntersection {
        bool hit = false;
        double ray_distance = 0.0;
        double surface_parameter = 0.0;
        Position position{0.0, 0.0, 0.0};
        DirectionVecter normal{0.0, 0.0, 0.0};
    };

    struct ConeIntersection {
        bool hit = false;
        size_t surface_index = 0;
        double ray_distance = 0.0;
        Position position{0.0, 0.0, 0.0};
        DirectionVecter normal{0.0, 0.0, 0.0};
    };

    struct RayTraceResult {
        double intensity = 0.0;
        bool exits_collector = false;
        Position position{0.0, 0.0, 0.0};
        DirectionVecter direction{0.0, 0.0, 0.0};
        int reflections = 0;
        double path_length = 0.0;
        bool reflection_limit_reached = false;
    };

    class Surface {
    public:
        Surface() = default;
        virtual ~Surface() = default;
        virtual SurfaceIntersection
        intersect(const Position &p, const DirectionVecter &d) const = 0;

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

        SurfaceIntersection
        intersect(const Position &pos,
                  const DirectionVecter &vec) const override {
            auto p = inverse_transform_->product(pos);
            auto v = inverse_transform_->product(vec);
            v.normalize_vector();
            // 平行于柱面方向的光线，不会与柱面相交
            if ((fabs(v.x_) < 1e-6) && (fabs(v.z_) < 1e-6)) {
                return {};
            }

            auto [delta, roots] = solve_equation(p.x_, v.x_, p.z_, v.z_);
            struct Candidate {
                double ray_distance = 0.0;
                double theta = 0.0;
                Position position;
            };
            std::vector<Candidate> positions;
            if (delta >= 0) {
                for (auto root : roots) {
                    if (!std::isfinite(root) || root < -1.0 || root > 1.0) {
                        continue;
                    }
                    root = std::clamp(root, -1.0, 1.0);
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
                    positions.push_back({distance, theta, Position(x, y, z)});
                }
            }
            if (positions.empty()) {
                return {};
            }
            const auto result_ptr = std::min_element(
                positions.begin(), positions.end(),
                [](const Candidate &p1, const Candidate &p2) {
                    return p1.ray_distance < p2.ray_distance;
                });
            if (!(result_ptr->position.x_ > 0 &&
                  result_ptr->position.z_ > 0)) {
                return {};
            }
            SurfaceIntersection result;
            result.hit = true;
            result.ray_distance = result_ptr->ray_distance;
            result.surface_parameter = result_ptr->theta;
            result.position = result_ptr->position.pre_multiply(*transform_);
            result.normal = normalAtTheta(result_ptr->theta);
            return result;
        }

        // 求解抛物柱面的法向量。
        // x=a*sin(theta-b)/(1-cos(theta))+c;z=a*cos(theta-b)/(1-cos(theta))。
        // 对 t 和 自由变量 y 求导，得到两个切向量，再求叉积即可得到法向量。
        DirectionVecter normalAtTheta(double theta) const {
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
            return n.normalize_vector();
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

        SurfaceIntersection
        intersect(const Position &pos,
                  const DirectionVecter &vec) const override {
            auto p = pos.pre_multiply(inverse_transform_);
            auto v = vec.pre_multiply(inverse_transform_).normalize_vector();
            // 平行于柱面方向的光线，不会与柱面相交
            if ((fabs(v.x_) < 1e-6) && (fabs(v.z_) < 1e-6)) {
                return {};
            }

            auto x = control_x_.pre_multiply(bezier_matrix_);
            auto z = control_z_.pre_multiply(bezier_matrix_);
            auto a = x[2] * v.z_ - z[2] * v.x_;
            auto b = x[1] * v.z_ - z[1] * v.x_;
            auto c = x[0] * v.z_ - z[0] * v.x_ - p.x_ * v.z_ + p.z_ * v.x_;
            // a == 0 的时候返回的大小根都是一样的，delta == 1.0
            auto [delta, min_root, max_root] = quadratic_equation_root(a, b, c);
            struct Candidate {
                double ray_distance = 0.0;
                double parameter = 0.0;
                Position position;
            };
            std::vector<Candidate> positions;
            if (delta >= 0) {
                for (auto root : {min_root, max_root}) {
                    // t 的范围一定是 0~1
                    if (!std::isfinite(root) || root < 0 || root > 1) {
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
                    positions.push_back(
                        {distance, root, Position(pos_x, pos_y, pos_z)});
                }
            }
            if (positions.empty()) {
                return {};
            }
            const auto result_ptr = std::min_element(
                positions.begin(), positions.end(),
                [](const Candidate &p1, const Candidate &p2) {
                    return p1.ray_distance < p2.ray_distance;
                });
            if (!(result_ptr->position.x_ > 0 &&
                  result_ptr->position.z_ > 0 &&
                  fabs(result_ptr->position.y_) < result_ptr->position.x_)) {
                return {};
            }

            auto x_coeff = control_x_.pre_multiply(bezier_matrix_);
            auto z_coeff = control_z_.pre_multiply(bezier_matrix_);
            MathUtils::D3Vecter derivative_vector(
                0, 1, result_ptr->parameter * 2);
            auto dx = x_coeff.dot_product(derivative_vector);
            auto dz = z_coeff.dot_product(derivative_vector);
            DirectionVecter vector_1{0.0, 1.0, 0.0};
            DirectionVecter vector_2{dx, 0, dz};
            auto normal = vector_1.cross_product(vector_2);
            if (!std::isfinite(normal.get_magtitude()) ||
                normal.get_magtitude() <= 1e-12) {
                return {};
            }

            SurfaceIntersection result;
            result.hit = true;
            result.ray_distance = result_ptr->ray_distance;
            result.surface_parameter = result_ptr->parameter;
            result.position = result_ptr->position.pre_multiply(transform_);
            result.normal =
                normal.normalize_vector().pre_multiply(transform_).normalize_vector();
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

        SurfaceIntersection
        intersect(const Position &p,
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
                result = p;
                auto x = (p - center_).dot_product(x_axis_) /
                         x_axis_.get_magtitude();
                auto y = (p - center_).dot_product(y_axis_) /
                         y_axis_.get_magtitude();
                is_intersect = std::fabs(up) < 1e-9 &&
                    ((fabs(x) < length_ / 2) && (fabs(y) < width_ / 2));
            }
            SurfaceIntersection intersection;
            intersection.hit = is_intersect && t >= -1e-9;
            intersection.ray_distance = intersection.hit
                ? get_distance(p, result)
                : 0.0;
            intersection.position = result;
            intersection.normal = normal_vector_;
            return intersection;
        }

        std::pair<bool, Position>
        is_in_plane(const Position &p, const DirectionVecter &v,
                    std::function<bool(double, double)> fn = nullptr) const {
            auto intersection = intersect(p, v);
            const auto &pos = intersection.position;
            const bool is_intersect = intersection.hit;
            auto dv = (pos - center_);
            auto x = dv.dot_product(x_axis_) / x_axis_.get_magtitude();
            auto y = dv.dot_product(y_axis_) / y_axis_.get_magtitude();
            if (fn == nullptr) {
                return std::make_pair(
                    is_intersect && (fabs(x) < length_ / 2) &&
                        (fabs(y) < width_ / 2),
                    pos);
            } else {
                return std::make_pair(is_intersect && fn(x, y), pos);
            }
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
            return std::clamp(reflectivity, 0.0, 1.0);
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

        ConeIntersection get_intersect_position(const Position &pos,
                                                DirectionVecter &vec) const {
            auto v = vec;
            v.normalize_vector();
            std::vector<ConeIntersection> positions;

            for (size_t index = 0; index < surfaces_.size(); index++) {
                const auto intersection = surfaces_[index]->intersect(pos, v);
                if (!intersection.hit ||
                    !std::isfinite(intersection.ray_distance) ||
                    intersection.ray_distance < 1e-6) {
                    continue;
                }
                ConeIntersection candidate;
                candidate.hit = true;
                candidate.surface_index = index;
                candidate.ray_distance = intersection.ray_distance;
                candidate.position = intersection.position;
                candidate.normal = intersection.normal;
                positions.push_back(candidate);
            }

            if (positions.empty()) {
                return {};
            }

            const auto result_ptr = std::min_element(
                positions.begin(), positions.end(),
                [](const ConeIntersection &p1, const ConeIntersection &p2) {
                    if (p1.ray_distance != p2.ray_distance) {
                        return p1.ray_distance < p2.ray_distance;
                    }
                    return p1.surface_index < p2.surface_index;
                });
            ConeIntersection result = *result_ptr;

            // At a collector seam/corner the mathematical surface normal is
            // not unique. Average every wall normal at the same first
            // intersection, using a deterministic tolerance and index tie
            // break above.
            DirectionVecter seam_normal{0.0, 0.0, 0.0};
            for (const auto &candidate : positions) {
                if (std::fabs(candidate.ray_distance -
                              result.ray_distance) < 1e-6) {
                    seam_normal = seam_normal + candidate.normal;
                }
            }
            const double seam_norm = seam_normal.get_magtitude();
            if (std::isfinite(seam_norm) && seam_norm > 1e-12) {
                result.normal = seam_normal.normalize_vector();
            }
            return result;
        }

        DirectionVecter get_reflect_direction(const DirectionVecter &v,
                                              DirectionVecter normal) const {
            auto ray = v;
            ray.normalize_vector();
            normal.normalize_vector();
            const auto product = ray.dot_product(normal);
            auto result = ray - 2 * product * normal;
            result.normalize_vector();
            return result;
        }

        RayTraceResult
        ray_trace_impl(const Position &pos, const DirectionVecter &vec,
                       const double in_intensity) const {
            Position p = pos;
            DirectionVecter v = vec;
            v.normalize_vector();
            double intensity = in_intensity;
            size_t number = 0;
            double path_length = 0.0;
            bool reflection_limit_reached = false;
            constexpr size_t max_reflections = 500;
            // 检查是否直接出去
            auto [is_on, position] = out_plane_->is_in_plane(p, v);
            bool is_out =
                (is_on && (out_plane_->normal_vector_.dot_product(v) < 0));
            if (is_out) {
                path_length += get_distance(p, position);
                p = position;
            }

            // 检查能否在入口平面上
            if (!is_out) {
                auto [is_on_in_plane, input_position] =
                    in_plane_->is_in_plane(p, v);
                bool is_out_in_plane =
                    is_on_in_plane &&
                    (in_plane_->normal_vector_.dot_product(v) > 0);

                if (!is_out_in_plane && is_on_in_plane) {
                    bool is_outside = false;
                    while (!is_outside && number < max_reflections) {
                        const auto intersection =
                            get_intersect_position(p, v);
                        if (!intersection.hit) {
                            is_out = false;
                            break;
                        }
                        path_length += intersection.ray_distance;
                        // 检查反射率
                        auto reflect_direction =
                            get_reflect_direction(v, intersection.normal);
                        const double reflected_dot = std::clamp(
                            reflect_direction.dot_product(v), -1.0, 1.0);
                        auto theta = 90 -
                            (180.0 / M_PI) * acos(reflected_dot) / 2.0;
                        auto is_reflective = material_->is_reflective(theta);
                        if (!is_reflective) {
                            is_out = false;
                            break;
                        }
                        intensity *= material_->get_reflect_intensity(theta);
                        if (!std::isfinite(intensity) || intensity <= 0.0) {
                            intensity = 0.0;
                            is_out = false;
                            break;
                        }
                        //  反射光线方向检查
                        auto reflect = material_->get_reflect_direction(
                            intersection.position, reflect_direction);
                        reflect.normalize_vector();
                        // 检查是否在离开
                        auto [_is_on_out_plane, out_position] =
                            out_plane_->is_in_plane(
                                intersection.position, reflect);
                        auto [_is_on_in_plane, in_position] =
                            in_plane_->is_in_plane(
                                intersection.position, reflect);
                        is_out_in_plane =
                            _is_on_in_plane &&
                            (in_plane_->normal_vector_.dot_product(reflect) >
                             0);
                        is_out = _is_on_out_plane &&
                                  (out_plane_->normal_vector_.dot_product(
                                       reflect) < 0);
                        is_outside = is_out || is_out_in_plane;
                        if (is_out) {
                            path_length +=
                                get_distance(intersection.position, out_position);
                            p = out_position;
                        } else if (is_out_in_plane) {
                            path_length +=
                                get_distance(intersection.position, in_position);
                            p = in_position;
                        } else {
                            p = intersection.position;
                        }
                        v = reflect;
                        number++;
                    }
                    reflection_limit_reached =
                        !is_outside && number >= max_reflections;
                }
                intensity = (number == 0) ? 0.0 : intensity;
            }
            RayTraceResult result;
            result.intensity = intensity;
            result.exits_collector = is_out;
            result.position = p;
            result.direction = v;
            result.reflections = static_cast<int>(number);
            result.path_length = path_length;
            result.reflection_limit_reached = reflection_limit_reached;
            return result;
        }

        inline bool ray_trace(Position &pos, DirectionVecter &vec,
                              double &intensity) const {
            const auto result = ray_trace_impl(pos, vec, intensity);
            intensity = result.intensity;
            pos = result.position;
            vec = result.direction;
            return result.exits_collector;
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
