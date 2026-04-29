#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "geometry/FacetLayoutUtils.hpp"
#include "io/MirrorFacetCsvReader.hpp"

namespace {

bool nearlyEqual(double a, double b, double tol = 1e-12) {
    return std::abs(a - b) <= tol;
}

bool check(bool condition, const std::string& label) {
    if (!condition) {
        std::cerr << "FAILED: " << label << "\n";
        return false;
    }
    return true;
}

bool writeTextFile(const std::string& path, const std::string& text) {
    std::ofstream ofs(path);
    if (!ofs) {
        return false;
    }
    ofs << text;
    return true;
}

} // namespace

int main() {
    bool ok = true;

    const std::string valid_path = "mirror_facets_valid.csv";
    ok &= check(writeTextFile(valid_path,
        "id,center_x,center_y,center_z,normal_x,normal_y,normal_z,"
        "surface_type,radius_of_curvature,aperture_shape,size1,size2,"
        "reflectivity_scale,roughness_sigma_rad,misalign_sigma_rad\n"
        "0,0,0,0,0,0,2,Spherical,10,Circular,0.22,0,0.98,0.001,0.002\n"
        "1,0.55,0,0.030342,-0.055084,0,0.998482,Spherical,10,Circular,0.22,0,1,0,0\n"),
        "write valid CSV fixture");

    std::vector<MirrorFacet> facets;
    std::string error;
    ok &= check(readMirrorFacetsCSV(valid_path, facets, &error), "read valid CSV: " + error);
    ok &= check(facets.size() == 2, "valid CSV facet count");

    if (facets.size() == 2) {
        ok &= check(facets[0].id == 0, "first imported id");
        ok &= check(nearlyEqual(facets[0].normal.norm(), 1.0), "import normal is normalized");
        ok &= check(facets[0].surface_type == SurfaceType::Spherical, "import surface type");
        ok &= check(facets[0].aperture_shape == ApertureShape::Circular, "import aperture shape");
        ok &= check(nearlyEqual(facets[0].radius_of_curvature, 10.0), "import curvature radius");
        ok &= check(nearlyEqual(facets[0].size1, 0.22), "import size1");
        ok &= check(nearlyEqual(facets[0].reflectivity_scale, 0.98), "import reflectivity scale");
        ok &= check(nearlyEqual(facets[0].roughness_sigma_rad, 0.001), "import roughness");
        ok &= check(nearlyEqual(facets[0].misalign_sigma_rad, 0.002), "import misalignment");

        MirrorLayout layout = makeMirrorLayoutFromFacets(facets);
        ok &= check(layout.size() == facets.size(), "imported facets convert to MirrorLayout");
        ok &= check(layout.tiles()[0].id == facets[0].id, "converted tile preserves id");
        ok &= check(nearlyEqual(layout.tiles()[0].normal.norm(), 1.0), "converted tile normal");
    }

    const std::string duplicate_path = "mirror_facets_duplicate.csv";
    ok &= check(writeTextFile(duplicate_path,
        "id,center_x,center_y,center_z,normal_x,normal_y,normal_z,"
        "surface_type,radius_of_curvature,aperture_shape,size1\n"
        "0,0,0,0,0,0,1,Spherical,10,Circular,0.22\n"
        "0,1,0,0,0,0,1,Spherical,10,Circular,0.22\n"),
        "write duplicate CSV fixture");

    std::vector<MirrorFacet> bad_facets;
    error.clear();
    ok &= check(!readMirrorFacetsCSV(duplicate_path, bad_facets, &error),
                "duplicate ids are rejected");
    ok &= check(error.find("duplicate id") != std::string::npos,
                "duplicate id error message");

    const std::string polynomial_path = "mirror_facets_polynomial.csv";
    ok &= check(writeTextFile(polynomial_path,
        "id,center_x,center_y,center_z,normal_x,normal_y,normal_z,"
        "surface_type,radius_of_curvature,aperture_shape,size1\n"
        "0,0,0,0,0,0,1,Polynomial,10,Circular,0.22\n"),
        "write polynomial CSV fixture");

    error.clear();
    ok &= check(!readMirrorFacetsCSV(polynomial_path, bad_facets, &error),
                "unsupported polynomial surface is rejected");
    ok &= check(error.find("Polynomial") != std::string::npos,
                "polynomial error message");

    if (ok) {
        std::cout << "Mirror facet CSV checks passed\n";
        return 0;
    }
    return 1;
}
