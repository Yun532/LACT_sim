#pragma once
#include <vector>
#include <limits>
#include <cmath>
#include <cstddef>
#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include "geometry/CameraPixel.hpp"

class CameraGeometry {
public:
    double focal_surface_z = 16.0;   // m

    void addPixel(const CameraPixel& p) {
        pixels_.push_back(p);
        spatial_index_valid_ = false;
    }

    const std::vector<CameraPixel>& pixels() const {
        return pixels_;
    }

    std::vector<CameraPixel>& pixels() {
        spatial_index_valid_ = false;
        return pixels_;
    }

    bool empty() const {
        return pixels_.empty();
    }

    std::size_t size() const {
        return pixels_.size();
    }

    static bool contains(const CameraPixel& p, double x, double y) {
        double dx = x - p.center.x;
        double dy = y - p.center.y;

        if (p.shape == PixelShape::Circular) {
            double r = 0.5 * p.size;
            return dx * dx + dy * dy <= r * r + 1e-14;
        }
        if (p.shape == PixelShape::Square) {
            double half = 0.5 * p.size;
            return std::abs(dx) <= half + 1e-14 && std::abs(dy) <= half + 1e-14;
        }
        if (p.shape == PixelShape::Hexagonal) {
            // size is flat-to-flat, so half the size is the apothem.
            double apothem = 0.5 * p.size;
            double q1 = std::abs(dx);
            double q2 = std::abs(0.5 * dx + 0.8660254037844386 * dy);
            double q3 = std::abs(-0.5 * dx + 0.8660254037844386 * dy);
            return std::max(q1, std::max(q2, q3)) <= apothem + 1e-14;
        }
        return false;
    }

    const CameraPixel* findContainingPixelPtr(double x, double y) const {
        if (pixels_.empty()) return nullptr;

        ensureSpatialIndex();
        if (x < index_origin_x_ || x > index_max_x_ ||
            y < index_origin_y_ || y > index_max_y_) {
            return nullptr;
        }
        const auto it = spatial_index_.find(cellKey(cellCoordinate(x, index_origin_x_),
                                                    cellCoordinate(y, index_origin_y_)));
        if (it == spatial_index_.end()) return nullptr;

        double best_r2 = std::numeric_limits<double>::max();
        const CameraPixel* best = nullptr;

        for (const std::size_t index : it->second) {
            const auto& p = pixels_[index];
            if (!contains(p, x, y)) {
                continue;
            }
            double dx = x - p.center.x;
            double dy = y - p.center.y;
            double r2 = dx * dx + dy * dy;
            if (r2 < best_r2) {
                best_r2 = r2;
                best = &p;
            }
        }
        return best;
    }

    int findContainingPixel(double x, double y) const {
        const auto* pixel = findContainingPixelPtr(x, y);
        return pixel ? pixel->id : -1;
    }

    int findNearestPixel(double x, double y) const {
        if (pixels_.empty()) return -1;

        double best_r2 = std::numeric_limits<double>::max();
        int best_id = -1;

        for (const auto& p : pixels_) {
            double dx = x - p.center.x;
            double dy = y - p.center.y;
            double r2 = dx * dx + dy * dy;
            if (r2 < best_r2) {
                best_r2 = r2;
                best_id = p.id;
            }
        }
        return best_id;
    }

private:
    static std::int64_t cellKey(int ix, int iy) {
        const std::uint64_t ux = static_cast<std::uint32_t>(ix);
        const std::uint64_t uy = static_cast<std::uint32_t>(iy);
        return static_cast<std::int64_t>((ux << 32U) | uy);
    }

    int cellCoordinate(double value, double origin) const {
        return static_cast<int>(std::floor((value - origin) / index_cell_size_));
    }

    void ensureSpatialIndex() const {
        if (spatial_index_valid_) return;
        spatial_index_.clear();
        if (pixels_.empty()) {
            spatial_index_valid_ = true;
            return;
        }

        index_cell_size_ = 0.0;
        index_origin_x_ = std::numeric_limits<double>::max();
        index_origin_y_ = std::numeric_limits<double>::max();
        index_max_x_ = std::numeric_limits<double>::lowest();
        index_max_y_ = std::numeric_limits<double>::lowest();
        for (const auto& pixel : pixels_) {
            index_cell_size_ = std::max(index_cell_size_, pixel.size);
            const double radius = pixel.shape == PixelShape::Hexagonal
                                      ? pixel.size / std::sqrt(3.0)
                                      : 0.5 * pixel.size;
            index_origin_x_ = std::min(index_origin_x_, pixel.center.x - radius);
            index_origin_y_ = std::min(index_origin_y_, pixel.center.y - radius);
            index_max_x_ = std::max(index_max_x_, pixel.center.x + radius);
            index_max_y_ = std::max(index_max_y_, pixel.center.y + radius);
        }
        if (!(index_cell_size_ > 0.0) || !std::isfinite(index_cell_size_)) {
            index_cell_size_ = 1.0;
        }

        for (std::size_t index = 0; index < pixels_.size(); ++index) {
            const auto& pixel = pixels_[index];
            const double radius = pixel.shape == PixelShape::Hexagonal
                                      ? pixel.size / std::sqrt(3.0)
                                      : 0.5 * pixel.size;
            const int ix0 = cellCoordinate(pixel.center.x - radius, index_origin_x_);
            const int ix1 = cellCoordinate(pixel.center.x + radius, index_origin_x_);
            const int iy0 = cellCoordinate(pixel.center.y - radius, index_origin_y_);
            const int iy1 = cellCoordinate(pixel.center.y + radius, index_origin_y_);
            for (int ix = ix0; ix <= ix1; ++ix) {
                for (int iy = iy0; iy <= iy1; ++iy) {
                    spatial_index_[cellKey(ix, iy)].push_back(index);
                }
            }
        }
        spatial_index_valid_ = true;
    }

    std::vector<CameraPixel> pixels_;
    mutable bool spatial_index_valid_ = false;
    mutable double index_cell_size_ = 1.0;
    mutable double index_origin_x_ = 0.0;
    mutable double index_origin_y_ = 0.0;
    mutable double index_max_x_ = 0.0;
    mutable double index_max_y_ = 0.0;
    mutable std::unordered_map<std::int64_t, std::vector<std::size_t>> spatial_index_;
};
