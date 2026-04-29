#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "geometry/FacetLayoutUtils.hpp"
#include "io/MirrorFacetCsvReader.hpp"
#include "optics/OpticalEfficiency.hpp"
#include "optics/OpticalTracer.hpp"
#include "optics/OutputPlane.hpp"

namespace {

bool check(bool condition, const std::string& label) {
    if (!condition) {
        std::cerr << "FAILED: " << label << "\n";
        return false;
    }
    return true;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: test_mirror_1229_layout <mirror_1229_facets.csv>\n";
        return 2;
    }

    std::vector<MirrorFacet> facets;
    std::string error;
    if (!readMirrorFacetsCSV(argv[1], facets, &error)) {
        std::cerr << "failed to read mirror CSV: " << error << "\n";
        return 1;
    }

    MirrorLayout mirrors = makeMirrorLayoutFromFacets(facets);
    OutputPlane plane;
    plane.point = {0.0, 0.0, -8.0};
    plane.normal = {0.0, 0.0, 1.0};
    plane.buildLocalFrame();

    OpticalTracer tracer;
    OpticalEfficiency eff;

    bool ok = true;
    double sum_r2 = 0.0;
    double max_r = 0.0;

    for (const auto& facet : facets) {
        Photon photon;
        photon.pos = {facet.center.x, facet.center.y, 1.0};
        photon.dir = {0.0, 0.0, -1.0};
        photon.normalizeDirection();

        OpticalSurfaceHit hit = tracer.traceToPlane(photon, mirrors, plane, eff);
        double r = hit.hit_surface ? std::sqrt(hit.u_m * hit.u_m + hit.v_m * hit.v_m) : -1.0;
        if (hit.hit_surface) {
            sum_r2 += r * r;
            if (r > max_r) {
                max_r = r;
            }
        }

        ok &= check(hit.hit_mirror, "center ray hits a mirror");
        ok &= check(hit.hit_surface, "center ray reaches output plane");
        ok &= check(hit.mirror_id == facet.id, "center ray hits its own mirror");
        ok &= check(r >= 0.0 && r < 1e-8, "center ray reaches common focus");
    }

    double rms = facets.empty() ? 0.0 : std::sqrt(sum_r2 / facets.size());
    std::cout << "mirror_1229 facets=" << facets.size() << "\n";
    std::cout << "center_ray_rms_m=" << rms << "\n";
    std::cout << "center_ray_max_m=" << max_r << "\n";

    return ok ? 0 : 1;
}
