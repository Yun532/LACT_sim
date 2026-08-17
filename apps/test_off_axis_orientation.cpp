#include "app/OpticalSimCommon.hpp"
#include "io/MirrorFacetCsvReader.hpp"

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
    return true;
}

Vec3 beamDirection(double theta_deg, double phi_deg)
{
    const double theta = theta_deg * DEG_TO_RAD;
    const double phi = phi_deg * DEG_TO_RAD;
    return Vec3{
        std::sin(theta) * std::cos(phi),
        std::sin(theta) * std::sin(phi),
        -std::cos(theta)
    }.normalized();
}

std::pair<double, double> centroidForBeam(const MirrorLayout& mirrors,
                                          const OutputPlane& plane,
                                          const Vec3& dir)
{
    OpticalTracer tracer;
    double sum_u = 0.0;
    double sum_v = 0.0;
    int n = 0;
    for (const auto& tile : mirrors.tiles()) {
        Photon photon;
        photon.pos = tile.center - dir * 120.0;
        photon.dir = dir;
        photon.wavelength_nm = 400.0;
        photon.weight = 1.0;

        auto hit = tracer.traceToPlane(photon, mirrors, plane, OpticalEfficiency{});
        if (hit.hit_surface) {
            sum_u += hit.u_m;
            sum_v += hit.v_m;
            ++n;
        }
    }
    if (n == 0) {
        throw std::runtime_error("no off-axis test rays reached output plane");
    }
    return {sum_u / n, sum_v / n};
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "usage: test_off_axis_orientation <mirror_1229_facets.csv>\n";
        return 2;
    }

    std::vector<MirrorFacet> facets;
    std::string error;
    if (!readMirrorFacetsCSV(argv[1], facets, &error)) {
        std::cerr << "failed to read mirror csv: " << error << "\n";
        return 1;
    }
    MirrorLayout mirrors = makeMirrorLayoutFromFacets(facets);

    OutputPlane plane;
    plane.point = {0.0, 0.0, -8.0};
    plane.normal = {0.0, 0.0, -1.0};
    // Mirror +x is viewed camera -> mirror (West at north pointing), while
    // output +u is viewed mirror -> camera and points East.  The two axes are
    // therefore opposites; +v remains mirror-local sky-up.
    plane.u_axis = {-1.0, 0.0, 0.0};
    plane.v_axis = {0.0, 1.0, 0.0};

    const auto x_plus = centroidForBeam(mirrors, plane, beamDirection(1.0, 0.0));
    const auto x_minus = centroidForBeam(mirrors, plane, beamDirection(1.0, 180.0));
    const auto y_plus = centroidForBeam(mirrors, plane, beamDirection(1.0, 90.0));
    const auto y_minus = centroidForBeam(mirrors, plane, beamDirection(1.0, 270.0));

    bool ok = true;
    constexpr double min_shift_m = 0.05;
    ok &= check(x_plus.first < -min_shift_m,
                "east-source propagation (+mirror x/West) should image at -u");
    ok &= check(x_minus.first > min_shift_m,
                "west-source propagation (-mirror x/East) should image at +u");
    ok &= check(y_plus.second > min_shift_m, "photon +y should image at +v");
    ok &= check(y_minus.second < -min_shift_m, "photon -y should image at -v");
    ok &= check(std::abs(x_plus.second) < 0.03, "x off-axis should not strongly shift v");
    ok &= check(std::abs(y_plus.first) < 0.03, "y off-axis should not strongly shift u");

    std::cout << "centroid photon +x: u=" << x_plus.first << " v=" << x_plus.second << "\n";
    std::cout << "centroid photon -x: u=" << x_minus.first << " v=" << x_minus.second << "\n";
    std::cout << "centroid photon +y: u=" << y_plus.first << " v=" << y_plus.second << "\n";
    std::cout << "centroid photon -y: u=" << y_minus.first << " v=" << y_minus.second << "\n";
    return ok ? 0 : 1;
}
