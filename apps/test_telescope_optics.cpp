#include "app/TelescopeOpticsCache.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

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

    auto fallback_reflectivity = makeFacet();
    fallback_reflectivity.reflectivity_scale = 1.0;
    fallback_reflectivity.has_reflectivity_scale = false;
    auto legacy_reflectivity = makeFacet();
    legacy_reflectivity.id = 2;
    legacy_reflectivity.reflectivity_scale = 0.8;
    legacy_reflectivity.has_reflectivity_scale = true;
    std::vector<MirrorFacet> reflectivity_facets{
        fallback_reflectivity, legacy_reflectivity};
    std::map<std::string, std::string> efficiency_cfg{
        {"efficiency.mirror_reflectivity_scale", "0.91"}};
    applyFacetEfficiencyScales(reflectivity_facets, efficiency_cfg);
    if (std::abs(reflectivity_facets[0].reflectivity_scale - 0.91) > 1e-12 ||
        std::abs(reflectivity_facets[1].reflectivity_scale - 0.8) > 1e-12) {
        std::cerr << "uniform reflectivity fallback overrode a legacy per-facet value\n";
        return 1;
    }

    const std::string scale_csv = "facet_reflectivity_scales.csv";
    {
        std::ofstream ofs(scale_csv);
        ofs << "id,reflectivity_scale\n1,0.7\n2,0\n";
    }
    efficiency_cfg["efficiency.mirror_reflectivity_scale_csv"] = scale_csv;
    applyFacetEfficiencyScales(reflectivity_facets, efficiency_cfg);
    if (std::abs(reflectivity_facets[0].reflectivity_scale - 0.7) > 1e-12 ||
        reflectivity_facets[1].reflectivity_scale != 0.0) {
        std::cerr << "per-facet efficiency CSV was not authoritative\n";
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

    auto active_misalignment = makeFacet();
    active_misalignment.has_misalign_sigma_rad = true;
    auto explicit_zero_misalignment = makeFacet();
    explicit_zero_misalignment.id = 2;
    explicit_zero_misalignment.misalign_sigma_rad = 0.0;
    explicit_zero_misalignment.has_misalign_sigma_rad = true;
    const Vec3 explicit_zero_normal = explicit_zero_misalignment.normal;
    std::vector<MirrorFacet> mixed_misalignment{
        active_misalignment, explicit_zero_misalignment};
    ErrorConfig mixed_error;
    mixed_error.random_seed = 44;
    mixed_error.facet_normal_sigma_deg = 1.0;
    applyFacetErrors(mixed_misalignment, mixed_error);
    const auto& zero_after = mixed_misalignment[1].normal;
    if (zero_after.x != explicit_zero_normal.x ||
        zero_after.y != explicit_zero_normal.y ||
        zero_after.z != explicit_zero_normal.z) {
        std::cerr << "global normal sigma was added to an active per-facet column\n";
        return 1;
    }

    auto zero_column = makeFacet();
    zero_column.misalign_sigma_rad = 0.0;
    zero_column.has_misalign_sigma_rad = true;
    const Vec3 zero_column_normal = zero_column.normal;
    std::vector<MirrorFacet> zero_column_facets{zero_column};
    ErrorConfig fallback_normal_error;
    fallback_normal_error.random_seed = 45;
    fallback_normal_error.facet_normal_sigma_deg = 0.1;
    applyFacetErrors(zero_column_facets, fallback_normal_error);
    const auto& fallback_after = zero_column_facets[0].normal;
    if (fallback_after.x == zero_column_normal.x &&
        fallback_after.y == zero_column_normal.y &&
        fallback_after.z == zero_column_normal.z) {
        std::cerr << "all-zero per-facet normal sigma did not use global fallback\n";
        return 1;
    }

    ErrorConfig scatter_error;
    scatter_error.reflect_direction_sigma_deg = 0.2;
    std::vector<MirrorFacet> zero_roughness{makeFacet()};
    zero_roughness[0].roughness_sigma_rad = 0.0;
    zero_roughness[0].has_roughness_sigma_rad = true;
    if (effectiveReflectDirectionSigmaRad(zero_roughness, scatter_error) <= 0.0) {
        std::cerr << "all-zero roughness did not use global fallback\n";
        return 1;
    }
    zero_roughness[0].roughness_sigma_rad = 0.002;
    if (effectiveReflectDirectionSigmaRad(zero_roughness, scatter_error) != 0.0) {
        std::cerr << "global reflection sigma was added to active per-facet roughness\n";
        return 1;
    }

    if (telescopeOpticsSeed(error.random_seed, 5) ==
        telescopeOpticsSeed(error.random_seed, 6)) {
        std::cerr << "telescope optical seeds collided\n";
        return 1;
    }
    return 0;
}
