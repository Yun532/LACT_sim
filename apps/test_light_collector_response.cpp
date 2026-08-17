#include "app/OpticalSimCommon.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
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
    {
        Cone::TrueReflectMaterial material;
        const double grazing_reflectivity =
            material.get_reflect_intensity(90.0);
        if (!check(std::isfinite(grazing_reflectivity) &&
                       grazing_reflectivity >= 0.0 &&
                       grazing_reflectivity <= 1.0,
                   "true-reflect grazing probability must stay in [0,1]")) {
            return 1;
        }
    }

    {
        const auto [delta, first, second] =
            MathUtils::quadratic_equation_root(0.0, 0.0, 1.0);
        if (!(delta < 0.0) || std::isfinite(first) ||
            std::isfinite(second)) {
            std::cerr << "degenerate quadratic equation was not rejected\n";
            return 1;
        }
        MathUtils::DMatrix<2, 3> original{1, 2, 3, 4, 5, 6};
        const MathUtils::DMatrix<2, 3> copied(original);
        if (copied(1, 2) != 6.0) {
            std::cerr << "non-square matrix copy changed its dimensions\n";
            return 1;
        }
    }

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
    ok &= check(direct.collector_path_length_m == 0.0 &&
                    std::abs(direct.time_ns - 12.0) < 1e-12,
                "camera response without a collector must not add optical path");

    auto backside = makeHit(0.0, 0.0, {0.0, 0.0, -1.0});
    applyCameraResponse(camera, nullptr, plane, sipm, electronics, backside);
    ok &= check(!backside.hit_camera,
                "ray crossing the camera plane from behind must be rejected");
    ok &= check(!backside.accepted && backside.pixel_id < 0,
                "camera backside rejection must happen before pixel response");

    auto outside = makeHit(0.02, 0.0, {0.0, 0.0, 1.0});
    applyCameraResponse(camera, nullptr, plane, sipm, electronics, outside);
    ok &= check(!outside.hit_camera, "point outside square pixel should not hit camera");

    ::lact::electronics::MicrocellConfig tiled;
    tiled.enabled = true;
    tiled.layout = "s17351_tiled_2x4";
    tiled.sensor_size_x_m = 0.0134;
    tiled.sensor_size_y_m = 0.0134;
    const double channel_fraction =
        ::lact::electronics::interChannelActiveFraction(tiled);
    auto missing_collector_coordinates =
        makeHit(0.0, 0.0, {0.0, 0.0, 1.0});
    try {
        applyCameraResponse(camera, nullptr, plane, sipm, electronics,
                            missing_collector_coordinates, 0.299792458,
                            &tiled, 1.0 / channel_fraction);
        std::cerr << "explicit S17351 silently accepted missing collector "
                     "exit coordinates\n";
        ok = false;
    } catch (...) {
    }
    auto conditioned = makeHit(0.0, 0.0, {0.0, 0.0, 1.0});
    conditioned.collector_exit_x_m = -0.0066875;
    conditioned.collector_exit_y_m = -0.0066875;
    applyCameraResponse(camera, nullptr, plane, sipm, electronics,
                        conditioned, 0.299792458, &tiled,
                        1.0 / channel_fraction);
    ok &= check(conditioned.accepted,
                "channel-area-conditioned hit should remain accepted");
    ok &= check(std::abs(conditioned.relative_efficiency -
                         1.0 / channel_fraction) < 1.0e-12,
                "expectation-mode channel PDE correction is incorrect");

    auto explicit_gap = makeHit(0.0, 0.0, {0.0, 0.0, 1.0});
    explicit_gap.collector_exit_x_m = 0.0;
    explicit_gap.collector_exit_y_m = 0.0;
    explicit_gap.collector_exit_z_m = 1.0e-9;
    applyCameraResponse(camera, nullptr, plane, sipm, electronics,
                        explicit_gap, 0.299792458, &tiled,
                        1.0 / channel_fraction);
    ok &= check(!explicit_gap.accepted &&
                    explicit_gap.sipm_channel_gap_rejected,
                "explicit inter-channel gap must reject before PDE weighting");

    CameraConfig collector_cfg;
    collector_cfg.enabled = true;
    collector_cfg.mode = "csv";
    collector_cfg.collector = "bezier";
    collector_cfg.collector_material = "true_reflect";
    collector_cfg.collector_exit_size_m = 0.0130;
    collector_cfg.collector_height_m = 0.0297;
    auto collector = buildLightCollector(collector_cfg, camera);
    ElectronicsResponse ideal_electronics;

    CameraGeometry hex_camera;
    CameraPixel hex_pixel;
    hex_pixel.id = 1;
    hex_pixel.center = {0.0, 0.0, 0.0};
    hex_pixel.size = 0.0244;
    hex_pixel.shape = PixelShape::Hexagonal;
    hex_camera.addPixel(hex_pixel);
    try {
        (void)buildLightCollector(collector_cfg, hex_camera);
        ok = false;
        std::cerr << "square-cone collector accepted a hexagonal pixel\n";
    } catch (...) {
    }

    CameraConfig overlapping_cfg;
    overlapping_cfg.enabled = true;
    overlapping_cfg.mode = "square_grid";
    overlapping_cfg.pixel_shape = "square";
    overlapping_cfg.pixel_size_m = 0.02;
    overlapping_cfg.pixel_pitch_m = 0.01;
    overlapping_cfg.radius_m = 0.03;
    try {
        (void)buildCameraGeometry(overlapping_cfg);
        ok = false;
        std::cerr << "camera geometry accepted overlapping pixel interiors\n";
    } catch (...) {
    }

    const auto duplicate_camera_path =
        std::filesystem::temp_directory_path() /
        "lact_duplicate_camera_ids.csv";
    {
        std::ofstream output(duplicate_camera_path);
        output << "id,x_m,y_m,shape,size_m\n"
               << "3,0,0,square,0.01\n"
               << "3,0.02,0,square,0.01\n";
    }
    try {
        (void)readCameraCsv(duplicate_camera_path.string());
        ok = false;
        std::cerr << "camera CSV accepted duplicate pixel ids\n";
    } catch (...) {
    }
    std::filesystem::remove(duplicate_camera_path);

    auto normal = makeHit(0.0, 0.0, {0.0, 0.0, 1.0});
    applyCameraResponse(camera, collector.get(), plane, sipm, ideal_electronics, normal);
    ok &= check(normal.hit_camera, "normal center ray should reach the SiPM");
    ok &= check(normal.collector_enabled, "normal center ray should use collector");
    ok &= check(normal.collector_reflections == 0,
                "normal center ray should pass through without wall reflections");
    ok &= check(std::abs(normal.collector_intensity - 1.0) < 1e-12,
                "normal center ray should keep full collector weight");
    ok &= check(normal.collector_path_length_m > 0.0,
                "normal center ray should report its collector path length");
    ok &= check(std::abs(normal.time_ns -
                        (12.0 + normal.collector_path_length_m / 0.299792458)) <
                    1e-12,
                "collector path length should be included in photon time");
    ok &= check(std::abs(normal.collector_time_delay_ns -
                        normal.collector_path_length_m / 0.299792458) < 1e-12,
                "collector should expose the same time delay added to the hit");
    ok &= check(!normal.collector_reflection_limit_reached,
                "normal center ray must not reach the collector reflection limit");

    auto angled = makeHit(0.010, 0.0, {0.35, 0.0, 1.0});
    applyCameraResponse(camera, collector.get(), plane, sipm, ideal_electronics, angled);
    ok &= check(angled.hit_camera, "angled ray should still reach the SiPM");
    ok &= check(angled.collector_reflections > 0,
                "angled ray should reflect inside the collector");
    ok &= check(angled.collector_intensity < normal.collector_intensity,
                "true_reflect material should reduce angled-ray collector weight");
    ok &= check(angled.collector_path_length_m > 0.0 &&
                    angled.time_ns > 12.0,
                "angled collector ray should report positive path and delay");
    ok &= check(!angled.collector_reflection_limit_reached,
                "ordinary angled ray must not reach the reflection limit");

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
