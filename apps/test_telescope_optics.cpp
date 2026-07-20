#include "app/TelescopeOpticsCache.hpp"

#include <cmath>
#include <iostream>

using namespace lact;

namespace {

MirrorFacet makeFacet()
{
    MirrorFacet facet;
    facet.id = 1;
    facet.center = {1.0, 0.0, -8.0};
    facet.normal = {0.0, 0.0, 1.0};
    facet.radius_of_curvature = 16.0;
    facet.reflectivity_scale = 0.8;
    facet.roughness_sigma_rad = 0.002;
    facet.misalign_sigma_rad = 0.001;
    return facet;
}

} // namespace

int main()
{
    ErrorConfig error;
    error.random_seed = 1234;
    error.facet_radial_position_sigma_m = 0.01;
    error.facet_normal_sigma_deg = 0.1;
    error.radius_of_curvature_sigma_m = 0.05;
    error.reflectivity_scale_sigma = 0.02;

    TelescopeOpticsCache cache({makeFacet()}, error);
    const auto& first = cache.layoutFor(5);
    const auto& repeat = cache.layoutFor(5);
    const auto& other = cache.layoutFor(6);
    if (&first != &repeat || first.size() != 1 || other.size() != 1) {
        std::cerr << "telescope optics cache is not stable\n";
        return 1;
    }
    const auto& a = first.tiles().front();
    const auto& b = other.tiles().front();
    if (a.center.x == b.center.x && a.center.y == b.center.y &&
        a.center.z == b.center.z && a.normal.x == b.normal.x &&
        a.normal.y == b.normal.y && a.normal.z == b.normal.z) {
        std::cerr << "different telescopes reused one optical realization\n";
        return 1;
    }
    if (a.roughness_sigma_rad != 0.002) {
        std::cerr << "facet roughness was not preserved\n";
        return 1;
    }

    auto zero_reflectivity = makeFacet();
    zero_reflectivity.reflectivity_scale = 0.0;
    std::vector<MirrorFacet> zero_facets{zero_reflectivity};
    ErrorConfig reflectivity_error;
    reflectivity_error.random_seed = 42;
    reflectivity_error.reflectivity_scale_sigma = 0.2;
    applyFacetErrors(zero_facets, reflectivity_error);
    if (zero_facets.front().reflectivity_scale != 0.0) {
        std::cerr << "reflectivity error replaced rather than scaled base value\n";
        return 1;
    }

    auto misaligned = makeFacet();
    const Vec3 original_normal = misaligned.normal;
    std::vector<MirrorFacet> misaligned_facets{misaligned};
    ErrorConfig facet_only_error;
    facet_only_error.random_seed = 43;
    applyFacetErrors(misaligned_facets, facet_only_error);
    const Vec3 changed_normal = misaligned_facets.front().normal;
    if (changed_normal.x == original_normal.x &&
        changed_normal.y == original_normal.y &&
        changed_normal.z == original_normal.z) {
        std::cerr << "per-facet misalignment was ignored\n";
        return 1;
    }

    if (telescopeOpticsSeed(error.random_seed, 5) ==
        telescopeOpticsSeed(error.random_seed, 6)) {
        std::cerr << "telescope optical seeds collided\n";
        return 1;
    }
    return 0;
}
