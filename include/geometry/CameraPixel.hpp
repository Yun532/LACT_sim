#pragma once
#include "core/Vec3.hpp"

enum class PixelShape {
    Circular,
    Hexagonal,
    Square
};

struct CameraPixel {
    int id = -1;
    Vec3 center;             // typically on the focal plane
    PixelShape shape = PixelShape::Circular;
    double size = 0.03;      // characteristic size in m
    double depth = 0.0;      // reserved for later
};
