#pragma once
#include <vector>
#include "geometry/MirrorTile.hpp"

class MirrorLayout {
public:
    void addTile(const MirrorTile& tile) {
        tiles_.push_back(tile);
    }

    const std::vector<MirrorTile>& tiles() const {
        return tiles_;
    }

    std::vector<MirrorTile>& tiles() {
        return tiles_;
    }

    bool empty() const {
        return tiles_.empty();
    }

    std::size_t size() const {
        return tiles_.size();
    }

private:
    std::vector<MirrorTile> tiles_;
};
