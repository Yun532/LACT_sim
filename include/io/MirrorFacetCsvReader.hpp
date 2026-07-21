#pragma once
#include <string>
#include <vector>
#include "geometry/MirrorFacet.hpp"

struct MirrorFacetCsvDefaults {
    SurfaceType surface_type = SurfaceType::Spherical;
    double radius_of_curvature = 16.0;
    ApertureShape aperture_shape = ApertureShape::Circular;
    double size1 = 0.25;
    double size2 = 0.0;
    double aperture_rotation_rad = 0.0;
};

bool readMirrorFacetsCSV(const std::string& path,
                         std::vector<MirrorFacet>& facets,
                         std::string* error = nullptr,
                         const MirrorFacetCsvDefaults* defaults = nullptr);
