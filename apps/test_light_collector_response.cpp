#include "app/OpticalSimCommon.hpp"

#include <cmath>
#include <iostream>

using namespace lact;

namespace {

bool check(bool cond, const std::string& msg)
{
    if (!cond) {
        std::cerr << msg << "\n";
        return false;
    }
    return cond;
}

CameraGeometry onePixelCamera()
{
    CameraGeometry camera;
    CameraPixel p;
    p.id = 7;
    p.center = {0.0, 0.0, 0.0};
    p.size = 0.0244;
    p.shape = PixelShape::Square;
    camera.addPixel(p);
    return camera;
}

OutputPlane outputPlane()
{
    OutputPlane plane;
    plane.point = {0.0, 0.0, -8.0};
    plane.normal = {0.0, 0.0, -1.0};
    plane.u_axis = {1.0, 0.0, 0.0};
    plane.v_axis = {0.0, 1.0, 0.0};
    return plane;
}

OpticalSurfaceHit makeHit(double u, double v, const Vec3& out_dir)
{
    OpticalSurfaceHit hit;
    hit.hit_surface = true;
    hit.u_m = u;
    hit.v_m = v;
    hit.out_dir = out_dir.normalized();
    hit.weight = 1.0;
    hit.relative_efficiency = 1.0;
    hit.wavelength_nm = 400.0;
    hit.time_ns = 12.0;
    return hit;
}

} // namespace

int main()
{
    const auto camera = onePixelCamera();
    const auto plane = outputPlane();
    SipmConfig sipm;
    sipm.size_m = 0.0130;

    ElectronicsResponse electronics;

    bool ok = true;
    auto direct = makeHit(0.0, 0.0, {0.0, 0.0, 1.0});
    applyCameraResponse(camera, nullptr, plane, sipm, electronics, direct);
    ok &= check(direct.hit_camera, "direct center hit should enter the pixel");
    ok &= check(direct.accepted, "direct center hit should be accepted");
    ok &= check(direct.pixel_id == 7, "direct center hit should map to pixel id 7");
    ok &= check(std::abs(direct.relative_efficiency - 1.0) < 1e-12,
                "direct center hit should keep ideal electronics weight");

    auto outside = makeHit(0.02, 0.0, {0.0, 0.0, 1.0});
    applyCameraResponse(camera, nullptr, plane, sipm, electronics, outside);
    ok &= check(!outside.hit_camera, "point outside square pixel should not hit camera");

    CameraConfig collector_cfg;
    collector_cfg.enabled = true;
    collector_cfg.mode = "csv";
    collector_cfg.collector = "bezier";
    collector_cfg.collector_material = "true_reflect";
    collector_cfg.collector_exit_size_m = 0.0130;
    collector_cfg.collector_height_m = 0.0297;
    auto collector = buildLightCollector(collector_cfg, camera);
    ElectronicsResponse ideal_electronics;

    auto normal = makeHit(0.0, 0.0, {0.0, 0.0, 1.0});
    applyCameraResponse(camera, collector.get(), plane, sipm, ideal_electronics, normal);
    ok &= check(normal.hit_camera, "normal center ray should reach the SiPM");
    ok &= check(normal.collector_enabled, "normal center ray should use collector");
    ok &= check(normal.collector_reflections == 0,
                "normal center ray should pass through without wall reflections");
    ok &= check(std::abs(normal.collector_intensity - 1.0) < 1e-12,
                "normal center ray should keep full collector weight");

    auto angled = makeHit(0.010, 0.0, {0.35, 0.0, 1.0});
    applyCameraResponse(camera, collector.get(), plane, sipm, ideal_electronics, angled);
    ok &= check(angled.hit_camera, "angled ray should still reach the SiPM");
    ok &= check(angled.collector_reflections > 0,
                "angled ray should reflect inside the collector");
    ok &= check(angled.collector_intensity < normal.collector_intensity,
                "true_reflect material should reduce angled-ray collector weight");

    std::cout << "collector normal ray: reflections="
              << normal.collector_reflections
              << " weight=" << normal.collector_intensity << "\n";
    std::cout << "collector angled ray: reflections="
              << angled.collector_reflections
              << " weight=" << angled.collector_intensity << "\n";

    int tried = 0;
    int accepted = 0;
    for (double u = -0.010; u <= 0.010; u += 0.0025) {
        for (double v = -0.010; v <= 0.010; v += 0.0025) {
            for (double du : {-0.25, 0.0, 0.25}) {
                for (double dv : {-0.25, 0.0, 0.25}) {
                    auto h = makeHit(u, v, {du, dv, 1.0});
                    applyCameraResponse(camera, collector.get(), plane, sipm,
                                        ideal_electronics, h);
                    ++tried;
                    if (h.hit_camera && h.accepted && h.collector_enabled &&
                        h.collector_intensity > 0.0 &&
                        h.collector_intensity <= 1.0 + 1e-12) {
                        ++accepted;
                    }
                }
            }
        }
    }
    ok &= check(tried > 0, "collector scan should try rays");
    ok &= check(accepted > 0, "bezier collector should accept at least one sampled ray");

    std::cout << "collector sampled rays=" << tried << " accepted=" << accepted << "\n";
    return ok ? 0 : 1;
}
