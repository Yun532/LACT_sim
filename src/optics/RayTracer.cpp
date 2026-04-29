#include "optics/RayTracer.hpp"
#include <cmath>
#include <limits>
#include <algorithm>
#include "optics/Reflection.hpp"

namespace {
constexpr double C_M_PER_NS = 0.299792458;
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
}

OpticalHit RayTracer::trace(const Photon& photon,
                            const MirrorLayout& mirrors,
                            const CameraGeometry& camera,
                            const OpticalEfficiency& eff) const
{
    OpticalHit hit;

    if (mirrors.empty()) return hit;

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

    auto cam_sol = intersectPlaneZ(best_sol.point, out_dir, camera.focal_surface_z);
    if (!cam_sol.has_value()) {
        return hit;
    }

    const auto& [tc, cam_p] = *cam_sol;
    hit.hit_camera = true;
    hit.camera_point = cam_p;

    hit.xcam = cam_p.x;
    hit.ycam = cam_p.y;
    hit.sxcam = out_dir.x;
    hit.sycam = out_dir.y;

    double total_path_m = best_sol.t + tc;
    hit.time_to_camera_ns = photon.time_ns + total_path_m / C_M_PER_NS;

    hit.pixel_id = camera.findNearestPixel(cam_p.x, cam_p.y);

    double cosang = std::clamp(std::abs(out_dir.z), 0.0, 1.0);
    double incidence_angle = std::acos(cosang);
    hit.relative_efficiency = eff.total(photon.wavelength_nm, incidence_angle);

    hit.accepted = (hit.pixel_id >= 0 && hit.relative_efficiency > 0.0);
    return hit;
}

std::optional<RayTracer::MirrorIntersection>
RayTracer::intersectMirror(const Vec3& p0, const Vec3& d, const MirrorTile& tile)
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

std::optional<RayTracer::MirrorIntersection>
RayTracer::intersectPlaneDisk(const Vec3& p0, const Vec3& d, const MirrorTile& tile)
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

    // 平面 facet：先直接用中心到命中点距离做口径裁切
    if (dp.norm() > tile.aperture_radius || !insideAperture(dp, n, tile)) {
        return std::nullopt;
    }

    MirrorIntersection out;
    out.t = t;
    out.point = p;
    out.normal = n;
    return out;
}

std::optional<RayTracer::MirrorIntersection>
RayTracer::intersectSphericalFacet(const Vec3& p0, const Vec3& d, const MirrorTile& tile)
{
    // --------------------------------------------------
    // 中文说明：
    // 把每块 facet 看成一个“球面镜片的小球冠”。
    //
    // 已知：
    //   1) 镜片中心点 tile.center
    //   2) 镜片中心法向 tile.normal
    //   3) 镜片焦距 tile.focal_length
    //
    // 球面镜近轴关系：R = 2 f
    // 这里用 tile.normal 来定义球心方向：
    //   sphere_center = tile.center - n * R
    //
    // 这样在镜片中心附近，球面法向与 tile.normal 一致。
    // --------------------------------------------------


    Vec3 n0 = tile.normal.normalized();
    double R = tile.radius_of_curvature;
    
    Vec3 sphere_center = tile.center + n0 * R;

    // 解射线与球面的交点：
    // | p0 + t d - c |^2 = R^2
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

    // --------------------------------------------------
    // 对完整球体会有两个交点。
    // 我们只接受“落在镜片小球冠区域内”的那个点。
    // 这里用相对于镜片中心切平面的横向距离做裁切。
    // 同时，如果两个根都合法，就选“更靠近镜片中心顶点”的那个。
    // --------------------------------------------------
    auto try_root = [&](double t) -> std::optional<MirrorIntersection> {
        if (t <= 0.0) return std::nullopt;

        Vec3 p = p0 + d * t;
        Vec3 rel = p - tile.center;

        // 在镜片中心切平面上做横向投影
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

        // 球面的几何法向
        out.normal = (sphere_center - p).normalized();
        return out;
    };

    auto sol1 = try_root(t1);
    auto sol2 = try_root(t2);

    if (!sol1 && !sol2) {
        return std::nullopt;
    }
    if (sol1 && !sol2) {
        return sol1;
    }
    if (!sol1 && sol2) {
        return sol2;
    }

    // 两个都合法时，选更靠近镜片顶点（tile.center）的那一个
    double a1 = std::abs((sol1->point - tile.center).dot(n0));
    double a2 = std::abs((sol2->point - tile.center).dot(n0));

    return (a1 < a2) ? sol1 : sol2;
}

std::optional<RayTracer::MirrorIntersection>
RayTracer::intersectParaboloid(const Vec3& p0, const Vec3& d, const MirrorTile& tile)
{
    // 连续理想抛物面测试：
    // 顶点在 tile.center，轴沿 +z
    double f = tile.radius_of_curvature;
    if (f <= 0.0) {
        return std::nullopt;
    }

    Vec3 q0 = p0 - tile.center;

    double x0 = q0.x;
    double y0 = q0.y;
    double z0 = q0.z;

    double dx = d.x;
    double dy = d.y;
    double dz = d.z;

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

    double r2 = q.x * q.x + q.y * q.y;
    if (std::sqrt(r2) > tile.aperture_radius || !insideAperture(q, tile.normal, tile)) {
        return std::nullopt;
    }

    // F = x^2 + y^2 - 4 f z = 0
    Vec3 n_local{2.0 * q.x, 2.0 * q.y, -4.0 * f};

    MirrorIntersection out;
    out.t = t;
    out.point = p;
    out.normal = n_local.normalized();
    return out;
}

std::optional<std::pair<double, Vec3>>
RayTracer::intersectPlaneZ(const Vec3& p0, const Vec3& d, double zplane)
{
    if (std::abs(d.z) < EPS) {
        return std::nullopt;
    }

    double t = (zplane - p0.z) / d.z;
    if (t <= 0.0) {
        return std::nullopt;
    }

    Vec3 p = p0 + d * t;
    return std::make_optional(std::make_pair(t, p));
}
