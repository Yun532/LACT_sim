#include "optics/OpticalEfficiency.hpp"
#include "optics/OpticalTracer.hpp"

#include <cmath>
#include <cstdint>
#include <iostream>

int main()
{
    constexpr std::size_t sample_count = 200000;
    constexpr double sigma_rad = 0.01;

    MirrorTile mirror;
    mirror.id = 9;
    mirror.center = {0.0, 0.0, 0.0};
    mirror.normal = {0.0, 0.0, 1.0};
    mirror.aperture_radius = 2.0;
    mirror.type = SurfaceType::Planar;
    MirrorLayout layout;
    layout.addTile(mirror);

    OutputPlane plane;
    plane.point = {0.0, 0.0, 5.0};
    plane.normal = {0.0, 0.0, -1.0};
    plane.buildLocalFrame();

    OpticalTracer tracer(0.299792458, sigma_rad, 987654321ULL);
    OpticalEfficiency efficiency;
    long double sum_x = 0.0L;
    long double sum_y = 0.0L;
    long double sum_x2 = 0.0L;
    long double sum_y2 = 0.0L;
    std::size_t accepted = 0;

    for (std::size_t i = 0; i < sample_count; ++i) {
        Photon photon;
        photon.pos = {0.0, 0.0, 1.0};
        photon.dir = {0.0, 0.0, -1.0};
        photon.random_stream_id = static_cast<std::uint64_t>(i + 1);
        const auto hit = tracer.traceToPlane(
            photon, layout, plane, efficiency);
        if (!hit.hit_surface) {
            continue;
        }
        sum_x += hit.out_dir.x;
        sum_y += hit.out_dir.y;
        sum_x2 += hit.out_dir.x * hit.out_dir.x;
        sum_y2 += hit.out_dir.y * hit.out_dir.y;
        ++accepted;
    }

    if (accepted != sample_count) {
        std::cerr << "not all scattered rays reached the output plane\n";
        return 1;
    }
    const long double inverse_count =
        1.0L / static_cast<long double>(accepted);
    const double mean_x = static_cast<double>(sum_x * inverse_count);
    const double mean_y = static_cast<double>(sum_y * inverse_count);
    const double variance_x = static_cast<double>(
        sum_x2 * inverse_count - mean_x * mean_x);
    const double variance_y = static_cast<double>(
        sum_y2 * inverse_count - mean_y * mean_y);
    const double sigma2 = sigma_rad * sigma_rad;

    std::cout << "accepted=" << accepted
              << " mean_x=" << mean_x
              << " mean_y=" << mean_y
              << " var_x/sigma2=" << variance_x / sigma2
              << " var_y/sigma2=" << variance_y / sigma2 << '\n';

    if (std::abs(mean_x) >= 8.0e-5 || std::abs(mean_y) >= 8.0e-5 ||
        std::abs(variance_x / sigma2 - 1.0) >= 0.02 ||
        std::abs(variance_y / sigma2 - 1.0) >= 0.02) {
        std::cerr << "mirror-scatter ray distribution is inconsistent with "
                     "the configured angular sigma\n";
        return 1;
    }
    return 0;
}
