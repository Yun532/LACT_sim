#pragma once
#include <string>
#include <vector>
#include "geometry/MirrorFacet.hpp"

bool readMirrorFacetsCSV(const std::string& path,
                         std::vector<MirrorFacet>& facets,
                         std::string* error = nullptr);
