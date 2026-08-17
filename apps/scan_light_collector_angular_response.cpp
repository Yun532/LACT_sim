#include "app/OpticalSimCommon.hpp"

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>

using namespace lact;

namespace {

constexpr double PI = 3.14159265358979323846;

struct Args {
    std::string output = "run_logs/collector_angular_response/collector_angular_response.csv";
    int photons_per_angle = 20000;
    double angle_step_deg = 1.0;
    double max_angle_deg = 90.0;
    std::uint64_t seed = 1229;
};

void usage(const char* argv0)
{
    std::cerr
        << "Usage: " << argv0 << " [options]\n"
        << "Options:\n"
        << "  --output PATH              CSV output path\n"
        << "  --photons-per-angle N      number of entrance photons per angle\n"
        << "  --angle-step-deg DEG       angular scan step\n"
        << "  --max-angle-deg DEG        maximum incidence angle\n"
        << "  --seed N                   random seed\n";
}

Args parseArgs(int argc, char** argv)
{
    Args args;
    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        auto need_value = [&](const std::string& name) -> std::string {
            if (i + 1 >= argc) {
                throw std::runtime_error(name + " requires a value");
            }
            return argv[++i];
        };
        if (key == "--output") {
            args.output = need_value(key);
        } else if (key == "--photons-per-angle") {
            args.photons_per_angle = std::stoi(need_value(key));
        } else if (key == "--angle-step-deg") {
            args.angle_step_deg = std::stod(need_value(key));
        } else if (key == "--max-angle-deg") {
            args.max_angle_deg = std::stod(need_value(key));
        } else if (key == "--seed") {
            args.seed = static_cast<std::uint64_t>(std::stoull(need_value(key)));
        } else if (key == "--help" || key == "-h") {
            usage(argv[0]);
            std::exit(0);
        } else {
            throw std::runtime_error("unknown argument: " + key);
        }
    }
    if (args.photons_per_angle <= 0) {
        throw std::runtime_error("--photons-per-angle must be positive");
    }
    if (!(args.angle_step_deg > 0.0)) {
        throw std::runtime_error("--angle-step-deg must be positive");
    }
    if (args.max_angle_deg < 0.0 || args.max_angle_deg > 90.0) {
        throw std::runtime_error("--max-angle-deg must be in [0, 90]");
    }
    return args;
}

CameraGeometry onePixelCamera()
{
    CameraGeometry camera;
    CameraPixel p;
    p.id = 0;
    p.center = {0.0, 0.0, 0.0};
    p.size = 0.0244;
    p.shape = PixelShape::Square;
    camera.addPixel(p);
    return camera;
}

OutputPlane collectorEntrancePlane()
{
    OutputPlane plane;
    plane.point = {0.0, 0.0, -8.0};
    plane.normal = {0.0, 0.0, -1.0};
    plane.u_axis = {1.0, 0.0, 0.0};
    plane.v_axis = {0.0, 1.0, 0.0};
    return plane;
}

OpticalSurfaceHit makeEntranceHit(double u_m, double v_m, double theta_rad, double phi_rad)
{
    OpticalSurfaceHit hit;
    hit.hit_surface = true;
    hit.u_m = u_m;
    hit.v_m = v_m;
    const double s = std::sin(theta_rad);
    const double c_raw = std::cos(theta_rad);
    const double c = std::abs(c_raw) < 1e-12 ? 0.0 : c_raw;
    hit.out_dir = Vec3{s * std::cos(phi_rad), s * std::sin(phi_rad), c}.normalized();
    hit.weight = 1.0;
    hit.relative_efficiency = 1.0;
    hit.wavelength_nm = 400.0;
    return hit;
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const Args args = parseArgs(argc, argv);

        const auto camera = onePixelCamera();
        const auto plane = collectorEntrancePlane();

        CameraConfig collector_cfg;
        collector_cfg.enabled = true;
        collector_cfg.mode = "csv";
        collector_cfg.collector = "bezier";
        collector_cfg.collector_material = "true_reflect";
        collector_cfg.collector_exit_size_m = 0.0130;
        collector_cfg.collector_height_m = 0.0297;
        auto collector = buildLightCollector(collector_cfg, camera);

        SipmConfig sipm;
        sipm.size_m = 0.0130;
        ElectronicsResponse ideal_electronics;

        const std::filesystem::path output_path(args.output);
        if (output_path.has_parent_path()) {
            std::filesystem::create_directories(output_path.parent_path());
        }
        std::ofstream out(args.output);
        if (!out) {
            throw std::runtime_error("failed to open output: " + args.output);
        }
        out << std::setprecision(10);
        out << "angle_deg,input_photons,hit_sipm_photons,geometric_acceptance,"
            << "weighted_photons,weighted_acceptance,mean_collector_weight,"
            << "mean_reflections,mean_reflections_all\n";

        std::mt19937_64 rng(args.seed);
        std::uniform_real_distribution<double> unit(0.0, 1.0);
        const double half_entrance_m = 0.5 * camera.pixels().front().size;

        for (double angle = 0.0; angle <= args.max_angle_deg + 1e-9;
             angle += args.angle_step_deg) {
            const double theta = angle * PI / 180.0;
            int hit_count = 0;
            double weighted_sum = 0.0;
            double reflection_sum_hit = 0.0;
            double reflection_sum_all = 0.0;

            for (int i = 0; i < args.photons_per_angle; ++i) {
                const double u = (2.0 * unit(rng) - 1.0) * half_entrance_m;
                const double v = (2.0 * unit(rng) - 1.0) * half_entrance_m;
                const double phi = 2.0 * PI * unit(rng);

                auto hit = makeEntranceHit(u, v, theta, phi);
                applyCameraResponse(camera, collector.get(), plane, sipm,
                                    ideal_electronics, hit);

                reflection_sum_all += hit.collector_reflections;
                if (hit.hit_camera && hit.accepted) {
                    ++hit_count;
                    weighted_sum += hit.collector_intensity;
                    reflection_sum_hit += hit.collector_reflections;
                }
            }

            const double n = static_cast<double>(args.photons_per_angle);
            const double hit_n = static_cast<double>(hit_count);
            const double geom_acceptance = hit_n / n;
            const double weighted_acceptance = weighted_sum / n;
            const double mean_weight = hit_count > 0 ? weighted_sum / hit_n : 0.0;
            const double mean_ref_hit = hit_count > 0 ? reflection_sum_hit / hit_n : 0.0;
            const double mean_ref_all = reflection_sum_all / n;

            out << angle << ","
                << args.photons_per_angle << ","
                << hit_count << ","
                << geom_acceptance << ","
                << weighted_sum << ","
                << weighted_acceptance << ","
                << mean_weight << ","
                << mean_ref_hit << ","
                << mean_ref_all << "\n";
        }

        std::cout << "collector angular response written: " << args.output << "\n"
                  << "photons_per_angle=" << args.photons_per_angle
                  << " angle_step_deg=" << args.angle_step_deg
                  << " material=true_reflect collector=bezier\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "scan_light_collector_angular_response error: "
                  << e.what() << "\n";
        return 1;
    }
}
