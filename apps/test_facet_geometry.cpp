#include <cmath>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "geometry/FacetFactory.hpp"
#include "geometry/MirrorFacetValidation.hpp"
#include "optics/Reflection.hpp"

namespace {

bool nearlyEqual(double a, double b, double tol = 1e-9) {
    return std::abs(a - b) <= tol;
}

bool check(bool condition, const std::string& label) {
    if (!condition) {
        std::cerr << "FAILED: " << label << "\n";
        return false;
    }
    return true;
}

FacetGridConfig defaultGrid() {
    FacetGridConfig grid;
    grid.facet_spacing = 0.55;
    grid.facet_radius = 0.22;
    grid.aperture_shape = ApertureShape::Circular;
    grid.surface_type = SurfaceType::Spherical;
    return grid;
}

DishPrescription defaultDish(DishType type) {
    DishPrescription dish;
    dish.type = type;
    dish.telescope_focal_length = 5.0;
    dish.dish_shape_length = 5.0;
    dish.dish_radius = 2.0;
    dish.vertex = {0.0, 0.0, 0.0};
    dish.optical_axis = {0.0, 0.0, 1.0};
    return dish;
}

bool checkCommonLayout(const std::vector<MirrorFacet>& facets,
                       const DishPrescription& dish,
                       const FacetGridConfig& grid)
{
    bool ok = true;
    std::string error;
    ok &= check(validateMirrorFacets(facets, &error), "generated facets validate: " + error);
    ok &= check(!facets.empty(), "generated facets are non-empty");

    std::set<int> ids;
    Vec3 focus = dish.vertex + Vec3{0.0, 0.0, dish.telescope_focal_length};
    Vec3 incoming{0.0, 0.0, -1.0};

    for (std::size_t i = 0; i < facets.size(); ++i) {
        const auto& facet = facets[i];
        ids.insert(facet.id);

        ok &= check(facet.id == static_cast<int>(i), "facet ids are contiguous");
        ok &= check(nearlyEqual(facet.size1, grid.facet_radius), "facet radius is copied");
        ok &= check(facet.aperture_shape == grid.aperture_shape, "aperture shape is copied");
        ok &= check(facet.surface_type == grid.surface_type, "surface type is copied");

        double dx = facet.center.x - dish.vertex.x;
        double dy = facet.center.y - dish.vertex.y;
        double r = std::sqrt(dx * dx + dy * dy);
        ok &= check(r <= dish.dish_radius + 1e-12, "facet center lies inside dish radius");
        ok &= check(nearlyEqual(facet.normal.norm(), 1.0), "facet normal is unit length");

        Vec3 reflected = reflectDirection(incoming, facet.normal.normalized());
        Vec3 to_focus = (focus - facet.center).normalized();
        ok &= check((reflected - to_focus).norm() < 1e-9,
                    "facet center ray reflects to focus");
    }

    ok &= check(ids.size() == facets.size(), "facet ids are unique");
    return ok;
}

bool checkDaviesCottonLayout() {
    DishPrescription dish = defaultDish(DishType::DaviesCotton);
    FacetGridConfig grid = defaultGrid();
    auto facets = FacetFactory::buildFacets(dish, grid);

    bool ok = checkCommonLayout(facets, dish, grid);
    ok &= check(facets.size() == 45, "Davies-Cotton fixture facet count");

    Vec3 focus = dish.vertex + Vec3{0.0, 0.0, dish.telescope_focal_length};
    for (const auto& facet : facets) {
        Vec3 rel = facet.center - focus;
        ok &= check(nearlyEqual(rel.norm(), dish.dish_shape_length, 1e-9),
                    "Davies-Cotton center lies on parent sphere");
        ok &= check(facet.center.z <= focus.z, "Davies-Cotton center is below focus");
        ok &= check(nearlyEqual(facet.radius_of_curvature,
                                2.0 * dish.telescope_focal_length),
                    "Davies-Cotton facet curvature radius");
    }

    return ok;
}

bool checkParabolicLayout() {
    DishPrescription dish = defaultDish(DishType::Parabolic);
    FacetGridConfig grid = defaultGrid();
    auto facets = FacetFactory::buildFacets(dish, grid);

    bool ok = checkCommonLayout(facets, dish, grid);
    ok &= check(facets.size() == 45, "parabolic fixture facet count");

    Vec3 focus = dish.vertex + Vec3{0.0, 0.0, dish.telescope_focal_length};
    for (const auto& facet : facets) {
        double dx = facet.center.x - dish.vertex.x;
        double dy = facet.center.y - dish.vertex.y;
        double expected_z = dish.vertex.z + (dx * dx + dy * dy) / (4.0 * dish.dish_shape_length);
        ok &= check(nearlyEqual(facet.center.z, expected_z, 1e-9),
                    "parabolic center lies on parent paraboloid");

        double expected_r = 2.0 * (focus - facet.center).norm();
        ok &= check(nearlyEqual(facet.radius_of_curvature, expected_r, 1e-9),
                    "parabolic facet curvature radius");
    }

    return ok;
}

} // namespace

int main() {
    bool ok = true;
    ok &= checkDaviesCottonLayout();
    ok &= checkParabolicLayout();

    if (ok) {
        std::cout << "Facet geometry checks passed\n";
        return 0;
    }
    return 1;
}
