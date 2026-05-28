#include "optics/OpticalTracer.hpp"
#include <cmath>
#include <limits>
#include <algorithm>
#include <functional>
#include <random>
#include "optics/Reflection.hpp"

namespace {
constexpr double EPS = 1e-14;

bool insideAperture(const Vec3& rel, const Vec3& normal, const MirrorTile& tile)
{
    Vec3 n = normal.normalized();
    Vec3 u = tile.aperture_u_axis.norm2() > 0.0
        ? tile.aperture_u_axis.normalized()
        : ((std::abs(n.z) < 0.9) ? Vec3{0.0, 0.0, 1.0}
                                 : Vec3{0.0, 1.0, 0.0}).cross(n).normalized();
    Vec3 v = tile.aperture_v_axis.norm2() > 0.0
        ? tile.aperture_v_axis.normalized()
        : n.cross(u).normalized();

    double x = rel.dot(u);
    double y = rel.dot(v);
    double c = std::cos(tile.aperture_rotation_rad);
    double s = std::sin(tile.aperture_rotation_rad);
    double xr = c * x + s * y;
    double yr = -s * x + c * y;

    if (tile.aperture_shape == ApertureShape::Circular) {
        return xr * xr + yr * yr <= tile.aperture_size1 * tile.aperture_size1 + 1e-12;
    }
    if (tile.aperture_shape == ApertureShape::Hexagon) {
        double apothem = 0.5 * tile.aperture_size1;
        double q1 = std::abs(xr);
        double q2 = std::abs(0.5 * xr + 0.8660254037844386 * yr);
        double q3 = std::abs(-0.5 * xr + 0.8660254037844386 * yr);
        return std::max(q1, std::max(q2, q3)) <= apothem + 1e-12;
    }
    if (tile.aperture_shape == ApertureShape::Square) {
        double half = 0.5 * tile.aperture_size1;
        return std::abs(xr) <= half + 1e-12 && std::abs(yr) <= half + 1e-12;
    }

    return false;
}

std::uint64_t mixSeed(std::uint64_t seed, std::uint64_t value)
{
    seed ^= value + 0x9E3779B97F4A7C15ULL + (seed << 6) + (seed >> 2);
    return seed;
}

std::uint64_t hashDouble(double value)
{
    return static_cast<std::uint64_t>(std::hash<double>{}(value));
}

Vec3 perturbDirection(const Vec3& direction, double sigma_rad, std::uint64_t seed)
{
    if (sigma_rad <= 0.0) {
        return direction.normalized();
    }

    Vec3 w = direction.normalized();
    Vec3 ref = (std::abs(w.z) < 0.9) ? Vec3{0.0, 0.0, 1.0}
                                     : Vec3{0.0, 1.0, 0.0};
    Vec3 u = ref.cross(w).normalized();
    Vec3 v = w.cross(u).normalized();

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double theta = sigma_rad * std::sqrt(-2.0 * std::log(std::max(1e-16, 1.0 - uniform(rng))));
    theta = std::min(theta, 3.14159265358979323846);
    double phi = 2.0 * 3.14159265358979323846 * uniform(rng);

    return (w * std::cos(theta) +
            u * (std::sin(theta) * std::cos(phi)) +
            v * (std::sin(theta) * std::sin(phi))).normalized();
}
}

OpticalTracer::OpticalTracer(double speed_of_light_m_per_ns,
                             double reflect_direction_sigma_rad,
                             std::uint64_t random_seed)
    : speed_of_light_m_per_ns_(speed_of_light_m_per_ns),
      reflect_direction_sigma_rad_(reflect_direction_sigma_rad),
      random_seed_(random_seed)
{
}

OpticalSurfaceHit OpticalTracer::traceToPlane(const Photon& photon,
                                              const MirrorLayout& mirrors,
                                              const OutputPlane& plane_in,
                                              const OpticalEfficiency& eff) const
{
    OpticalSurfaceHit hit;

    if (mirrors.empty()) return hit;

    OutputPlane plane = plane_in;
    if (plane.u_axis.norm2() <= 0.0 || plane.v_axis.norm2() <= 0.0) {
        plane.buildLocalFrame();
    } else {
        plane.normal = plane.normal.normalized();
        plane.u_axis = plane.u_axis.normalized();
        plane.v_axis = plane.v_axis.normalized();
    }

    double best_t = std::numeric_limits<double>::max();
    const MirrorTile* best_tile = nullptr;
    MirrorIntersection best_sol;

    for (const auto& tile : mirrors.tiles()) {
        auto sol = intersectMirror(photon.pos, photon.dir, tile);
        if (!sol.has_value()) continue;

        if (sol->t < best_t) {
            best_t = sol->t;
            best_tile = &tile;
            best_sol = *sol;
        }
    }

    if (!best_tile) {
        return hit;
    }

    hit.hit_mirror = true;
    hit.mirror_id = best_tile->id;
    hit.mirror_point = best_sol.point;

    Vec3 n = best_sol.normal.normalized();
    Vec3 out_dir = reflectDirection(photon.dir, n);
    std::uint64_t scatter_seed = random_seed_;
    scatter_seed = mixSeed(scatter_seed, static_cast<std::uint64_t>(best_tile->id + 1));
    scatter_seed = mixSeed(scatter_seed, hashDouble(photon.pos.x));
    scatter_seed = mixSeed(scatter_seed, hashDouble(photon.pos.y));
    scatter_seed = mixSeed(scatter_seed, hashDouble(photon.pos.z));
    scatter_seed = mixSeed(scatter_seed, hashDouble(photon.dir.x));
    scatter_seed = mixSeed(scatter_seed, hashDouble(photon.dir.y));
    scatter_seed = mixSeed(scatter_seed, hashDouble(photon.dir.z));
    out_dir = perturbDirection(out_dir, reflect_direction_sigma_rad_, scatter_seed);

    auto surf_sol = intersectOutputPlane(best_sol.point, out_dir, plane);
    if (!surf_sol.has_value()) {
        return hit;
    }

    const auto& [ts, surf_p] = *surf_sol;
    hit.hit_surface = true;
    hit.surface_point = surf_p;
    hit.out_dir = out_dir;

    Vec3 rel = surf_p - plane.point;
    hit.u_m = rel.dot(plane.u_axis);
    hit.v_m = rel.dot(plane.v_axis);

    double total_path_m = best_sol.t + ts;
    hit.time_ns = photon.time_ns + total_path_m / speed_of_light_m_per_ns_;
    hit.wavelength_nm = photon.wavelength_nm;
    hit.weight = photon.weight;

    double cosang = std::clamp(std::abs(out_dir.dot(plane.normal.normalized())), 0.0, 1.0);
    double incidence_angle = std::acos(cosang);
    hit.relative_efficiency = best_tile->reflectivity_scale *
                              eff.total(photon.wavelength_nm, incidence_angle);

    return hit;
}

OpticalSurfaceHit OpticalTracer::traceBackprojectedToPlane(
    const Photon& photon,
    const MirrorLayout& mirrors,
    const OutputPlane& plane_in,
    const OpticalEfficiency& eff) const
{
    OpticalSurfaceHit hit;

    if (mirrors.empty()) return hit;

    OutputPlane plane = plane_in;
    if (plane.u_axis.norm2() <= 0.0 || plane.v_axis.norm2() <= 0.0) {
        plane.buildLocalFrame();
    } else {
        plane.normal = plane.normal.normalized();
        plane.u_axis = plane.u_axis.normalized();
        plane.v_axis = plane.v_axis.normalized();
    }

    const Vec3 reverse_dir = photon.dir * -1.0;
    double best_t = std::numeric_limits<double>::max();
    const MirrorTile* best_tile = nullptr;
    MirrorIntersection best_sol;

    for (const auto& tile : mirrors.tiles()) {
        auto sol = intersectMirror(photon.pos, reverse_dir, tile);
        if (!sol.has_value()) continue;

        if (sol->t < best_t) {
            best_t = sol->t;
            best_tile = &tile;
            best_sol = *sol;
        }
    }

    if (!best_tile) {
        return hit;
    }

    hit.hit_mirror = true;
    hit.mirror_id = best_tile->id;
    hit.mirror_point = best_sol.point;

    Vec3 n = best_sol.normal.normalized();
    Vec3 out_dir = reflectDirection(photon.dir, n);
    std::uint64_t scatter_seed = random_seed_;
    scatter_seed = mixSeed(scatter_seed, static_cast<std::uint64_t>(best_tile->id + 1));
    scatter_seed = mixSeed(scatter_seed, hashDouble(photon.pos.x));
    scatter_seed = mixSeed(scatter_seed, hashDouble(photon.pos.y));
    scatter_seed = mixSeed(scatter_seed, hashDouble(photon.pos.z));
    scatter_seed = mixSeed(scatter_seed, hashDouble(photon.dir.x));
    scatter_seed = mixSeed(scatter_seed, hashDouble(photon.dir.y));
    scatter_seed = mixSeed(scatter_seed, hashDouble(photon.dir.z));
    out_dir = perturbDirection(out_dir, reflect_direction_sigma_rad_, scatter_seed);

    auto surf_sol = intersectOutputPlane(best_sol.point, out_dir, plane);
    if (!surf_sol.has_value()) {
        return hit;
    }

    const auto& [ts, surf_p] = *surf_sol;
    hit.hit_surface = true;
    hit.surface_point = surf_p;
    hit.out_dir = out_dir;

    Vec3 rel = surf_p - plane.point;
    hit.u_m = rel.dot(plane.u_axis);
    hit.v_m = rel.dot(plane.v_axis);

    hit.time_ns = photon.time_ns + (ts - best_sol.t) / speed_of_light_m_per_ns_;
    hit.wavelength_nm = photon.wavelength_nm;
    hit.weight = photon.weight;

    double cosang = std::clamp(std::abs(out_dir.dot(plane.normal.normalized())), 0.0, 1.0);
    double incidence_angle = std::acos(cosang);
    hit.relative_efficiency = best_tile->reflectivity_scale *
                              eff.total(photon.wavelength_nm, incidence_angle);

    return hit;
}

std::optional<OpticalTracer::MirrorIntersection>
OpticalTracer::intersectMirror(const Vec3& p0, const Vec3& d, const MirrorTile& tile)
{
    switch (tile.type) {
        case SurfaceType::Planar:
            return intersectPlaneDisk(p0, d, tile);
        case SurfaceType::Spherical:
            return intersectSphericalFacet(p0, d, tile);
        case SurfaceType::Parabolic:
            return intersectParaboloid(p0, d, tile);
        default:
            return std::nullopt;
    }
}

std::optional<OpticalTracer::MirrorIntersection>
OpticalTracer::intersectPlaneDisk(const Vec3& p0, const Vec3& d, const MirrorTile& tile)
{
    Vec3 n = tile.normal.normalized();
    double denom = d.dot(n);

    if (std::abs(denom) < EPS) {
        return std::nullopt;
    }

    double t = (tile.center - p0).dot(n) / denom;
    if (t <= 0.0) {
        return std::nullopt;
    }

    Vec3 p = p0 + d * t;
    Vec3 dp = p - tile.center;

    if (dp.norm() > tile.aperture_radius || !insideAperture(dp, n, tile)) {
        return std::nullopt;
    }

    MirrorIntersection out;
    out.t = t;
    out.point = p;
    out.normal = n;
    return out;
}

std::optional<OpticalTracer::MirrorIntersection>
OpticalTracer::intersectSphericalFacet(const Vec3& p0, const Vec3& d, const MirrorTile& tile)
{
    // 这里直接使用镜片曲率半径 R，而不是先从 focal length 换算。
    double R = tile.radius_of_curvature;
    if (!std::isfinite(R) || R <= 0.0) {
        return std::nullopt;
    }

    Vec3 n0 = tile.normal.normalized();

    // 球心位置：沿镜片中心法向反方向退 R
    Vec3 sphere_center = tile.center + n0 * R;

    // 解射线与球面的交点
    Vec3 oc = p0 - sphere_center;

    double A = d.dot(d);
    double B = 2.0 * oc.dot(d);
    double C = oc.dot(oc) - R * R;

    double disc = B * B - 4.0 * A * C;
    if (disc < 0.0) {
        return std::nullopt;
    }

    double sqrt_disc = std::sqrt(std::max(0.0, disc));
    double t1 = (-B - sqrt_disc) / (2.0 * A);
    double t2 = (-B + sqrt_disc) / (2.0 * A);

    auto try_root = [&](double t) -> std::optional<MirrorIntersection> {
        if (t <= 0.0) return std::nullopt;

        Vec3 p = p0 + d * t;
        Vec3 rel = p - tile.center;

        // 在镜片中心切平面上做横向裁切
        Vec3 tangential = rel - n0 * rel.dot(n0);
        double rho = tangential.norm();
        if (rho > tile.aperture_radius || !insideAperture(rel, n0, tile)) {
            return std::nullopt;
        }

        // 只保留镜片中心附近的小球冠，排除同一完整球面的远端交点。
        double axial = rel.dot(n0);
        double sag_limit = R - std::sqrt(std::max(0.0, R * R - tile.aperture_radius * tile.aperture_radius));
        if (axial < -1e-9 || axial > sag_limit + 1e-9) {
            return std::nullopt;
        }

        MirrorIntersection out;
        out.t = t;
        out.point = p;
        out.normal = (sphere_center - p).normalized();
        return out;
    };

    auto sol1 = try_root(t1);
    auto sol2 = try_root(t2);

    if (!sol1 && !sol2) return std::nullopt;
    if (sol1 && !sol2) return sol1;
    if (!sol1 && sol2) return sol2;

    double a1 = std::abs((sol1->point - tile.center).dot(n0));
    double a2 = std::abs((sol2->point - tile.center).dot(n0));

    return (a1 < a2) ? sol1 : sol2;
}

std::optional<OpticalTracer::MirrorIntersection>
OpticalTracer::intersectParaboloid(const Vec3& p0, const Vec3& d, const MirrorTile& tile)
{
    // 连续理想抛物面测试
    // 这里 tile.radius_of_curvature 不适用，仍然需要一个抛物面焦距参数。
    // 为了兼容当前结构，先用 focal_length() 辅助函数。
    double f = tile.focal_length();
    if (!std::isfinite(f) || f <= 0.0) {
        return std::nullopt;
    }

    Vec3 axis = tile.normal.normalized();
    Vec3 u = tile.aperture_u_axis.norm2() > 0.0
        ? tile.aperture_u_axis.normalized()
        : ((std::abs(axis.z) < 0.9) ? Vec3{0.0, 0.0, 1.0}
                                    : Vec3{0.0, 1.0, 0.0}).cross(axis).normalized();
    Vec3 v = tile.aperture_v_axis.norm2() > 0.0
        ? tile.aperture_v_axis.normalized()
        : axis.cross(u).normalized();

    Vec3 q0 = p0 - tile.center;

    double x0 = q0.dot(u);
    double y0 = q0.dot(v);
    double z0 = q0.dot(axis);

    double dx = d.dot(u);
    double dy = d.dot(v);
    double dz = d.dot(axis);

    double A = dx * dx + dy * dy;
    double B = 2.0 * (x0 * dx + y0 * dy - 2.0 * f * dz);
    double C = x0 * x0 + y0 * y0 - 4.0 * f * z0;

    double t = -1.0;

    if (std::abs(A) < EPS) {
        if (std::abs(B) < EPS) {
            return std::nullopt;
        }
        t = -C / B;
        if (t <= 0.0) {
            return std::nullopt;
        }
    } else {
        double disc = B * B - 4.0 * A * C;
        if (disc < 0.0) {
            return std::nullopt;
        }

        double sqrt_disc = std::sqrt(std::max(0.0, disc));
        double t1 = (-B - sqrt_disc) / (2.0 * A);
        double t2 = (-B + sqrt_disc) / (2.0 * A);

        bool ok1 = (t1 > 0.0);
        bool ok2 = (t2 > 0.0);

        if (!ok1 && !ok2) {
            return std::nullopt;
        } else if (ok1 && ok2) {
            t = std::min(t1, t2);
        } else {
            t = ok1 ? t1 : t2;
        }
    }

    Vec3 p = p0 + d * t;
    Vec3 q = p - tile.center;
    double x = q.dot(u);
    double y = q.dot(v);
    double z = q.dot(axis);

    double r2 = x * x + y * y;
    if (std::sqrt(r2) > tile.aperture_radius || !insideAperture(q, tile.normal, tile)) {
        return std::nullopt;
    }

    if (z < -1e-9) {
        return std::nullopt;
    }

    Vec3 n_global = (u * (2.0 * x) + v * (2.0 * y) - axis * (4.0 * f)).normalized();

    MirrorIntersection out;
    out.t = t;
    out.point = p;
    out.normal = n_global;
    return out;
}

std::optional<std::pair<double, Vec3>>
OpticalTracer::intersectOutputPlane(const Vec3& p0, const Vec3& d, const OutputPlane& plane)
{
    Vec3 n = plane.normal.normalized();
    double denom = d.dot(n);

    if (std::abs(denom) < EPS) {
        return std::nullopt;
    }

    double t = (plane.point - p0).dot(n) / denom;
    if (t <= 0.0) {
        return std::nullopt;
    }

    Vec3 p = p0 + d * t;
    return std::make_optional(std::make_pair(t, p));
}
