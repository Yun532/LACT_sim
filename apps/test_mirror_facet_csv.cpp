#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "app/OpticalSimCommon.hpp"
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
    std::vector<MirrorFacet> bad_facets;
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
        ok &= check(facets[0].has_reflectivity_scale, "record explicit reflectivity scale");
        ok &= check(nearlyEqual(facets[0].roughness_sigma_rad, 0.001), "import roughness");
        ok &= check(facets[0].has_roughness_sigma_rad, "record explicit roughness");
        ok &= check(nearlyEqual(facets[0].misalign_sigma_rad, 0.002), "import misalignment");
        ok &= check(facets[0].has_misalign_sigma_rad, "record explicit misalignment");

        MirrorLayout layout = makeMirrorLayoutFromFacets(facets);
        ok &= check(layout.size() == facets.size(), "imported facets convert to MirrorLayout");
        ok &= check(layout.tiles()[0].id == facets[0].id, "converted tile preserves id");
        ok &= check(nearlyEqual(layout.tiles()[0].normal.norm(), 1.0), "converted tile normal");
    }

    const std::string optional_path = "mirror_facets_optional.csv";
    ok &= check(writeTextFile(optional_path,
        "id,center_x,center_y,center_z,normal_x,normal_y,normal_z,"
        "radius_of_curvature\n"
        "0,0,0,0,0,0,1,0\n"),
        "write optional-column CSV fixture");
    MirrorFacetCsvDefaults defaults;
    defaults.radius_of_curvature = 16.4;
    defaults.aperture_shape = ApertureShape::Hexagon;
    defaults.size1 = 0.8;
    facets.clear();
    error.clear();
    ok &= check(readMirrorFacetsCSV(optional_path, facets, &error, &defaults),
                "read CSV with optional optical fields: " + error);
    if (facets.size() == 1) {
        ok &= check(nearlyEqual(facets[0].radius_of_curvature, 16.4),
                    "zero curvature placeholder uses mirror default");
        ok &= check(facets[0].aperture_shape == ApertureShape::Hexagon,
                    "missing aperture shape uses mirror default");
        ok &= check(nearlyEqual(facets[0].size1, 0.8),
                    "missing aperture size uses mirror default");
        ok &= check(!facets[0].has_reflectivity_scale,
                    "missing reflectivity scale remains unresolved");
    }

    const std::string negative_radius_path = "mirror_facets_negative_radius.csv";
    ok &= check(writeTextFile(negative_radius_path,
        "id,center_x,center_y,center_z,normal_x,normal_y,normal_z,"
        "radius_of_curvature\n"
        "0,0,0,0,0,0,1,-1\n"),
        "write negative-radius CSV fixture");
    error.clear();
    ok &= check(!readMirrorFacetsCSV(
                    negative_radius_path, bad_facets, &error, &defaults),
                "negative curvature is rejected");
    ok &= check(error.find("radius_of_curvature must be finite and >= 0") !=
                    std::string::npos,
                "negative curvature error message");

    const std::string series_path = "mirror_facets_series_optional.csv";
    ok &= check(writeTextFile(series_path,
        "elevation_deg,id,center_x,center_y,center_z,normal_x,normal_y,normal_z,"
        "radius_of_curvature,reflectivity_scale,roughness_sigma_rad,misalign_sigma_rad\n"
        "40,0,0,0,0,0,0,1,16,0.9,0.001,0\n"
        "40,1,1,0,0,0,0,1,0,1,0,0\n"
        "60,0,0,0,0,0,0,1,16,0.8,0.003,0\n"
        "60,1,1,0,0,0,0,1,0,1,0,0\n"),
        "write elevation-series optional-column fixture");
    std::map<std::string, std::string> series_cfg{
        {"mirror.mode", "elevation_series"},
        {"mirror.series_format", "facet_csv"},
        {"mirror.series_csv_path", series_path},
        {"mirror.series_elevation_deg", "50"},
        {"mirror.radius_of_curvature", "17"},
        {"mirror.aperture_shape", "Hexagon"},
        {"mirror.size1", "0.8"}};
    const auto series_facets = lact::buildFacetsFromConfig(series_cfg);
    ok &= check(series_facets.size() == 2, "elevation-series facet count");
    if (series_facets.size() == 2) {
        ok &= check(nearlyEqual(series_facets[0].radius_of_curvature, 16.0),
                    "series explicit curvature is interpolated");
        ok &= check(nearlyEqual(series_facets[1].radius_of_curvature, 17.0),
                    "series zero curvature uses mirror default");
        ok &= check(nearlyEqual(series_facets[0].reflectivity_scale, 0.85),
                    "series reflectivity scale is interpolated");
        ok &= check(nearlyEqual(series_facets[0].roughness_sigma_rad, 0.002),
                    "series roughness is interpolated");
    }

    const std::string duplicate_path = "mirror_facets_duplicate.csv";
    ok &= check(writeTextFile(duplicate_path,
        "id,center_x,center_y,center_z,normal_x,normal_y,normal_z,"
        "surface_type,radius_of_curvature,aperture_shape,size1\n"
        "0,0,0,0,0,0,1,Spherical,10,Circular,0.22\n"
        "0,1,0,0,0,0,1,Spherical,10,Circular,0.22\n"),
        "write duplicate CSV fixture");

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
