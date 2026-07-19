#include "app/TelescopeOpticsCache.hpp"

#include "geometry/FacetLayoutUtils.hpp"

#include <utility>

namespace lact {

std::uint64_t telescopeOpticsSeed(std::uint64_t base_seed, int telescope_id)
{
    const auto value = static_cast<std::uint64_t>(
        static_cast<std::int64_t>(telescope_id) + 10007);
    base_seed ^= value + 0x9e3779b97f4a7c15ULL +
                 (base_seed << 6U) + (base_seed >> 2U);
    return base_seed;
}

TelescopeOpticsCache::TelescopeOpticsCache(
    std::vector<MirrorFacet> nominal_facets,
    ErrorConfig error_config)
    : nominal_facets_(std::move(nominal_facets)),
      error_config_(std::move(error_config))
{
}

const MirrorLayout& TelescopeOpticsCache::layoutFor(int telescope_id)
{
    const auto found = layouts_.find(telescope_id);
    if (found != layouts_.end()) {
        return found->second;
    }

    std::vector<MirrorFacet> facets = nominal_facets_;
    ErrorConfig telescope_error = error_config_;
    telescope_error.random_seed = telescopeOpticsSeed(
        error_config_.random_seed, telescope_id);
    applyFacetErrors(facets, telescope_error);
    auto inserted = layouts_.emplace(
        telescope_id, makeMirrorLayoutFromFacets(facets));
    return inserted.first->second;
}

} // namespace lact
