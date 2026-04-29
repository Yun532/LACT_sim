#include <cmath>
#include <iostream>
#include "geometry/CameraGeometry.hpp"

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

    std::cout << "camera geometry checks passed\n";
    return 0;
}
