#include <cmath>
#include <iostream>

#include "core/Photon.hpp"
#include "geometry/MirrorLayout.hpp"
#include "optics/OpticalEfficiency.hpp"
#include "optics/OpticalTracer.hpp"
#include "optics/OutputPlane.hpp"

namespace {

bool nearlyEqual(double a, double b, double tol = 1e-9) {
    return std::abs(a - b) <= tol;
}

bool check(bool condition, const char* label) {
    if (!condition) {
        std::cerr << "FAILED: " << label << "\n";
        return false;
    }
    return true;
}

Photon makePhoton(const Vec3& pos, const Vec3& dir) {
    Photon p;
    p.pos = pos;
    p.dir = dir;
    p.normalizeDirection();
    return p;
}

OutputPlane makeZPlane(double z) {
    OutputPlane plane;
    plane.point = {0.0, 0.0, z};
    plane.normal = {0.0, 0.0, 1.0};
    plane.buildLocalFrame();
    return plane;
}

} // namespace

int main() {
    OpticalTracer tracer;
    OpticalEfficiency eff;
    OutputPlane focal_plane = makeZPlane(5.0);

    bool ok = true;

    {
        MirrorTile mirror;
        mirror.id = 1;
        mirror.center = {0.0, 0.0, 0.0};
        mirror.normal = {0.0, 0.0, 1.0};
        mirror.aperture_radius = 1.0;
        mirror.type = SurfaceType::Planar;

        MirrorLayout layout;
        layout.addTile(mirror);

        auto hit = tracer.traceToPlane(makePhoton({0.2, -0.1, 10.0}, {0.0, 0.0, -1.0}),
                                       layout, focal_plane, eff);

        ok &= check(hit.hit_mirror, "planar mirror is hit");
        ok &= check(hit.hit_surface, "planar reflected ray reaches output plane");
        ok &= check(hit.mirror_id == 1, "planar mirror id is preserved");
        ok &= check(nearlyEqual(hit.mirror_point.x, 0.2), "planar mirror x intersection");
        ok &= check(nearlyEqual(hit.mirror_point.y, -0.1), "planar mirror y intersection");
        ok &= check(nearlyEqual(hit.mirror_point.z, 0.0), "planar mirror z intersection");
        ok &= check(nearlyEqual(hit.surface_point.x, 0.2), "planar output x");
        ok &= check(nearlyEqual(hit.surface_point.y, -0.1), "planar output y");
        ok &= check(nearlyEqual(hit.surface_point.z, 5.0), "planar output z");
    }

    {
        MirrorTile mirror;
        mirror.id = 2;
        mirror.center = {0.0, 0.0, 0.0};
        mirror.normal = {0.0, 0.0, 1.0};
        mirror.aperture_radius = 0.5;
        mirror.radius_of_curvature = 10.0;
        mirror.type = SurfaceType::Spherical;

        MirrorLayout layout;
        layout.addTile(mirror);

        auto hit = tracer.traceToPlane(makePhoton({0.0, 0.0, 30.0}, {0.0, 0.0, -1.0}),
                                       layout, focal_plane, eff);

        ok &= check(hit.hit_mirror, "spherical center ray hits mirror");
        ok &= check(hit.hit_surface, "spherical center ray reaches output plane");
        ok &= check(hit.mirror_id == 2, "spherical mirror id is preserved");
        ok &= check(nearlyEqual(hit.mirror_point.z, 0.0), "spherical center ray uses near cap");
        ok &= check(nearlyEqual(hit.surface_point.x, 0.0), "spherical center output x");
        ok &= check(nearlyEqual(hit.surface_point.y, 0.0), "spherical center output y");
        ok &= check(nearlyEqual(hit.surface_point.z, 5.0), "spherical center output z");
    }

    {
        MirrorTile near_mirror;
        near_mirror.id = 3;
        near_mirror.center = {-1.1, -1.65, 0.410065};
        near_mirror.normal = {0.112327, 0.168491, 0.979282};
        near_mirror.aperture_radius = 0.22;
        near_mirror.radius_of_curvature = 10.0;
        near_mirror.type = SurfaceType::Spherical;

        MirrorTile far_cap_trap;
        far_cap_trap.id = 4;
        far_cap_trap.center = {1.1, 1.65, 0.410065};
        far_cap_trap.normal = {-0.112327, -0.168491, 0.979282};
        far_cap_trap.aperture_radius = 0.22;
        far_cap_trap.radius_of_curvature = 10.0;
        far_cap_trap.type = SurfaceType::Spherical;

        MirrorLayout layout;
        layout.addTile(near_mirror);
        layout.addTile(far_cap_trap);

        auto hit = tracer.traceToPlane(makePhoton({-1.1, -1.65, 30.0}, {0.0, 0.0, -1.0}),
                                       layout, focal_plane, eff);

        ok &= check(hit.hit_mirror, "spherical far-cap regression hits a mirror");
        ok &= check(hit.hit_surface, "spherical far-cap regression reaches output plane");
        ok &= check(hit.mirror_id == 3, "spherical far-cap regression rejects opposite far cap");
        ok &= check(hit.mirror_point.z < 1.0, "spherical far-cap regression uses physical cap");
    }

    if (ok) {
        std::cout << "Raytrace primitive checks passed\n";
        return 0;
    }

    return 1;
}
