#include <cmath>
#include <iostream>
#include "geometry/CameraGeometry.hpp"

namespace {

int bruteForceContainingPixel(const CameraGeometry& camera, double x, double y) {
    const CameraPixel* best = nullptr;
    double best_r2 = std::numeric_limits<double>::max();
    for (const auto& pixel : camera.pixels()) {
        if (!CameraGeometry::contains(pixel, x, y)) continue;
        const double dx = x - pixel.center.x;
        const double dy = y - pixel.center.y;
        const double r2 = dx * dx + dy * dy;
        if (r2 < best_r2) {
            best_r2 = r2;
            best = &pixel;
        }
    }
    return best ? best->id : -1;
}

} // namespace

int main() {
    CameraGeometry camera;

    CameraPixel center;
    center.id = 10;
    center.center = {0.0, 0.0, 0.0};
    center.shape = PixelShape::Hexagonal;
    center.size = 0.10;
    camera.addPixel(center);

    CameraPixel right;
    right.id = 11;
    right.center = {0.12, 0.0, 0.0};
    right.shape = PixelShape::Circular;
    right.size = 0.08;
    camera.addPixel(right);

    if (camera.findContainingPixel(0.0, 0.0) != 10) {
        std::cerr << "expected center hex pixel\n";
        return 1;
    }
    if (camera.findContainingPixel(0.12, 0.0) != 11) {
        std::cerr << "expected right circular pixel\n";
        return 1;
    }
    if (camera.findContainingPixel(0.30, 0.0) != -1) {
        std::cerr << "point outside camera should not map to a pixel\n";
        return 1;
    }
    if (camera.findNearestPixel(0.30, 0.0) != 11) {
        std::cerr << "nearest-pixel fallback changed unexpectedly\n";
        return 1;
    }

    for (int ix = -80; ix <= 160; ++ix) {
        for (int iy = -90; iy <= 90; ++iy) {
            const double x = static_cast<double>(ix) * 0.002;
            const double y = static_cast<double>(iy) * 0.002;
            const int indexed = camera.findContainingPixel(x, y);
            const int brute = bruteForceContainingPixel(camera, x, y);
            if (indexed != brute) {
                std::cerr << "spatial index disagrees with brute force at "
                          << x << ", " << y << ": " << indexed
                          << " vs " << brute << "\n";
                return 1;
            }
        }
    }

    CameraPixel upper;
    upper.id = 12;
    upper.center = {0.0, 0.12, 0.0};
    upper.shape = PixelShape::Square;
    upper.size = 0.08;
    camera.addPixel(upper);
    if (camera.findContainingPixel(0.0, 0.12) != 12) {
        std::cerr << "spatial index was not invalidated after addPixel\n";
        return 1;
    }
    if (camera.findContainingPixel(1.0e100, -1.0e100) != -1) {
        std::cerr << "far out-of-camera lookup should return no pixel\n";
        return 1;
    }

    std::cout << "camera geometry checks passed\n";
    return 0;
}
