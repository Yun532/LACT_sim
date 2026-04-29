#pragma once
#include <vector>
#include <limits>
#include <cmath>
#include <cstddef>
#include <algorithm>
#include "geometry/CameraPixel.hpp"

class CameraGeometry {
public:
    double focal_surface_z = 16.0;   // m

    void addPixel(const CameraPixel& p) {
        pixels_.push_back(p);
    }

    const std::vector<CameraPixel>& pixels() const {
        return pixels_;
    }

    std::vector<CameraPixel>& pixels() {
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

    int findContainingPixel(double x, double y) const {
        if (pixels_.empty()) return -1;

        double best_r2 = std::numeric_limits<double>::max();
        int best_id = -1;

        for (const auto& p : pixels_) {
            if (!contains(p, x, y)) {
                continue;
            }
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
    std::vector<CameraPixel> pixels_;
};
