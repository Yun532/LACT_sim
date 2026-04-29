#pragma once
#include <cmath>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include "geometry/MirrorFacet.hpp"

inline bool validateMirrorFacets(const std::vector<MirrorFacet>& facets,
                                 std::string* error = nullptr)
{
    std::set<int> ids;

    for (const auto& facet : facets) {
        auto fail = [&](const std::string& msg) {
            if (error) {
                std::ostringstream oss;
                oss << "facet id " << facet.id << ": " << msg;
                *error = oss.str();
            }
            return false;
        };

        if (facet.id < 0) {
            return fail("id must be non-negative");
        }
        if (!ids.insert(facet.id).second) {
            return fail("duplicate id");
        }

        if (!std::isfinite(facet.center.x) ||
            !std::isfinite(facet.center.y) ||
            !std::isfinite(facet.center.z)) {
            return fail("center must be finite");
        }

        if (!std::isfinite(facet.normal.x) ||
            !std::isfinite(facet.normal.y) ||
            !std::isfinite(facet.normal.z) ||
            facet.normal.norm() <= 0.0) {
            return fail("normal must be finite and non-zero");
        }

        if (!std::isfinite(facet.size1) || facet.size1 <= 0.0) {
            return fail("size1 must be positive");
        }
        if (!std::isfinite(facet.size2)) {
            return fail("size2 must be finite");
        }
        if (!std::isfinite(facet.aperture_rotation_rad)) {
            return fail("aperture_rotation_rad must be finite");
        }

        if (facet.surface_type == SurfaceType::Spherical ||
            facet.surface_type == SurfaceType::Parabolic) {
            if (!std::isfinite(facet.radius_of_curvature) ||
                facet.radius_of_curvature <= 0.0) {
                return fail("curved surfaces require positive radius_of_curvature");
            }
        }

        if (facet.surface_type == SurfaceType::Polynomial) {
            return fail("Polynomial surface is declared but not supported by current raytrace");
        }

        if (!std::isfinite(facet.reflectivity_scale) || facet.reflectivity_scale < 0.0) {
            return fail("reflectivity_scale must be finite and non-negative");
        }
        if (!std::isfinite(facet.roughness_sigma_rad) || facet.roughness_sigma_rad < 0.0) {
            return fail("roughness_sigma_rad must be finite and non-negative");
        }
        if (!std::isfinite(facet.misalign_sigma_rad) || facet.misalign_sigma_rad < 0.0) {
            return fail("misalign_sigma_rad must be finite and non-negative");
        }
    }

    if (error) {
        error->clear();
    }
    return true;
}
