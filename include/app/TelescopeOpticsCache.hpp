#pragma once

#include <cstdint>
#include <map>
#include <vector>

#include "app/OpticalSimCommon.hpp"
#include "geometry/MirrorLayout.hpp"

namespace lact {

std::uint64_t telescopeOpticsSeed(std::uint64_t base_seed, int telescope_id);

class TelescopeOpticsCache {
public:
    TelescopeOpticsCache(std::vector<MirrorFacet> nominal_facets,
                         ErrorConfig error_config);

    const MirrorLayout& layoutFor(int telescope_id);

private:
    std::vector<MirrorFacet> nominal_facets_;
    ErrorConfig error_config_;
    std::map<int, MirrorLayout> layouts_;
};

} // namespace lact
