#pragma once
#include <vector>
#include "geometry/MirrorFacet.hpp"
#include "geometry/MirrorLayout.hpp"

// 把 MirrorFacet 列表转换成光追模块使用的 MirrorLayout
inline MirrorLayout makeMirrorLayoutFromFacets(const std::vector<MirrorFacet>& facets)
{
    MirrorLayout layout;
    for (const auto& f : facets) {
        layout.addTile(f.toMirrorTile());
    }
    return layout;
}
