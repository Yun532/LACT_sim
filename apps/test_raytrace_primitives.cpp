#include <cmath>
#include <iostream>

#include "app/OpticalSimCommon.hpp"
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

    {
        MirrorTile mirror;
        mirror.id = 5;
        mirror.center = {0.0, 0.0, 0.0};
        mirror.normal = {0.0, 0.0, 1.0};
        mirror.aperture_radius = 1.0;
        mirror.type = SurfaceType::Planar;

        MirrorLayout layout;
        layout.addTile(mirror);

        auto backside = tracer.traceToPlane(
            makePhoton({0.0, 0.0, -1.0}, {0.0, 0.0, 1.0}),
            layout, focal_plane, eff);
        ok &= check(!backside.hit_mirror,
                    "mirror backside incidence must not reflect");

        auto signed_front = tracer.traceBackprojectedToPlane(
            makePhoton({0.0, 0.0, -1.0}, {0.0, 0.0, -1.0}),
            layout, focal_plane, eff);
        ok &= check(signed_front.hit_mirror && signed_front.hit_surface,
                    "signed EventIO line may reach a front-facing mirror at negative t");

        auto signed_backside = tracer.traceBackprojectedToPlane(
            makePhoton({0.0, 0.0, 1.0}, {0.0, 0.0, 1.0}),
            layout, focal_plane, eff);
        ok &= check(!signed_backside.hit_mirror,
                    "signed-line tracing must not enable backside reflection");
    }

    {
        MirrorTile upstream;
        upstream.id = 7;
        upstream.center = {0.0, 0.0, 1.0};
        upstream.normal = {0.0, 0.0, 1.0};
        upstream.aperture_radius = 1.0;
        upstream.type = SurfaceType::Planar;

        MirrorTile downstream = upstream;
        downstream.id = 8;
        downstream.center = {0.0, 0.0, 0.0};

        MirrorLayout layout;
        layout.addTile(upstream);
        layout.addTile(downstream);

        auto anchor_between = tracer.traceBackprojectedToPlane(
            makePhoton({0.0, 0.0, 0.2}, {0.0, 0.0, -1.0}),
            layout, focal_plane, eff);
        auto anchor_above = tracer.traceBackprojectedToPlane(
            makePhoton({0.0, 0.0, 2.0}, {0.0, 0.0, -1.0}),
            layout, focal_plane, eff);
        ok &= check(anchor_between.hit_mirror && anchor_between.mirror_id == 7,
                    "signed line must select the first mirror, not nearest |t|");
        ok &= check(anchor_above.hit_mirror && anchor_above.mirror_id == 7,
                    "signed-line mirror selection must be anchor invariant");
        ok &= check(nearlyEqual(anchor_between.mirror_point.z,
                                anchor_above.mirror_point.z),
                    "anchor translation must preserve the mirror intersection");
    }

    {
        MirrorTile mirror;
        mirror.id = 6;
        mirror.center = {0.0, 0.0, 0.0};
        mirror.normal = {std::sqrt(0.75), 0.0, 0.5};
        mirror.aperture_radius = 1.0;
        mirror.type = SurfaceType::Planar;

        MirrorLayout layout;
        layout.addTile(mirror);
        auto reflected_away = tracer.traceToPlane(
            makePhoton({0.0, 0.0, 1.0}, {0.0, 0.0, -1.0}),
            layout, focal_plane, eff);
        ok &= check(reflected_away.hit_mirror,
                    "front-side ray should hit the tilted mirror");
        ok &= check(!reflected_away.hit_surface,
                    "reflection directed away from camera must not be back-intersected");
    }

    {
        OpticalTracer scattered_tracer(0.299792458, 0.01, 987654321ULL);
        MirrorTile mirror;
        mirror.id = 9;
        mirror.center = {0.0, 0.0, 0.0};
        mirror.normal = {0.0, 0.0, 1.0};
        mirror.aperture_radius = 2.0;
        mirror.type = SurfaceType::Planar;
        MirrorLayout base_layout;
        base_layout.addTile(mirror);
        auto base_plane = makeZPlane(5.0);
        auto base_photon = makePhoton({0.2, -0.1, 10.0},
                                      {0.0, 0.0, -1.0});
        base_photon.random_stream_id = 11223344ULL;
        const auto base_hit = scattered_tracer.traceToPlane(
            base_photon, base_layout, base_plane, eff);

        auto independent_photon = base_photon;
        independent_photon.random_stream_id += 1ULL;
        const auto independent_hit = scattered_tracer.traceToPlane(
            independent_photon, base_layout, base_plane, eff);
        const auto replay_hit = scattered_tracer.traceToPlane(
            base_photon, base_layout, base_plane, eff);

        ok &= check(independent_hit.hit_surface && replay_hit.hit_surface,
                    "roughness independence test rays should reach the output plane");
        ok &= check(!nearlyEqual(base_hit.out_dir.x, independent_hit.out_dir.x, 0.0) ||
                        !nearlyEqual(base_hit.out_dir.y, independent_hit.out_dir.y, 0.0) ||
                        !nearlyEqual(base_hit.out_dir.z, independent_hit.out_dir.z, 0.0),
                    "different photon stream ids must draw different roughness directions");
        ok &= check(nearlyEqual(base_hit.out_dir.x, replay_hit.out_dir.x, 0.0) &&
                        nearlyEqual(base_hit.out_dir.y, replay_hit.out_dir.y, 0.0) &&
                        nearlyEqual(base_hit.out_dir.z, replay_hit.out_dir.z, 0.0),
                    "the same photon stream id must reproduce its roughness direction");

        const Vec3 translation{13.0, -21.0, 7.0};
        mirror.center += translation;
        MirrorLayout translated_layout;
        translated_layout.addTile(mirror);
        auto translated_plane = base_plane;
        translated_plane.point += translation;
        auto translated_photon = base_photon;
        translated_photon.pos += translation;
        const auto translated_hit = scattered_tracer.traceToPlane(
            translated_photon, translated_layout, translated_plane, eff);

        ok &= check(base_hit.hit_surface && translated_hit.hit_surface,
                    "translated roughness test rays should reach the output plane");
        ok &= check(nearlyEqual(base_hit.out_dir.x, translated_hit.out_dir.x, 0.0) &&
                        nearlyEqual(base_hit.out_dir.y, translated_hit.out_dir.y, 0.0) &&
                        nearlyEqual(base_hit.out_dir.z, translated_hit.out_dir.z, 0.0),
                    "roughness random draw must be coordinate-translation invariant");
        ok &= check(nearlyEqual(base_hit.u_m, translated_hit.u_m) &&
                        nearlyEqual(base_hit.v_m, translated_hit.v_m) &&
                        nearlyEqual(base_hit.time_ns, translated_hit.time_ns),
                    "translated optical system should preserve local hit and travel time");
    }

    {
        lact::ObstructionMask obstruction;
        obstruction.enabled = true;
        obstruction.mode = "primitives";
        obstruction.check_incoming = true;
        obstruction.check_reflected = true;

        lact::ObstructionMask::Primitive camera_body;
        camera_body.type = "cylinder";
        camera_body.role = "camera_body";
        camera_body.p0 = {0.0, 0.0, -11.5};
        camera_body.p1 = {0.0, 0.0, -10.5};
        camera_body.radius_m = 0.1;
        obstruction.primitives.push_back(camera_body);

        ok &= check(lact::incomingRayBlockedByObstruction(
                        {0.0, 0.0, -10.0}, {0.0, 0.0, 1.0}, obstruction),
                    "incoming upstream ray must find an obstruction behind its anchor");
        ok &= check(!lact::incomingRayBlockedByObstruction(
                        {0.0, 0.0, -10.0}, {0.0, 0.0, -1.0}, obstruction),
                    "incoming upstream ray must preserve photon direction");
        ok &= check(lact::incomingSegmentBlockedByObstruction(
                        {0.0, 0.0, -12.0}, {0.0, 0.0, -10.0}, obstruction),
                    "incoming leg must not be inferred from the sign of local dz");
        ok &= check(!lact::segmentBlockedByObstruction(
                        {0.0, 0.0, -10.0}, {0.0, 0.0, -12.0}, obstruction),
                    "camera body must not self-block the reflected leg");

        obstruction.primitives.front().role = "support_strut";
        ok &= check(lact::segmentBlockedByObstruction(
                        {0.0, 0.0, -10.0}, {0.0, 0.0, -12.0}, obstruction),
                    "support strut must be checked on the reflected leg");
    }

    if (ok) {
        std::cout << "Raytrace primitive checks passed\n";
        return 0;
    }

    return 1;
}
