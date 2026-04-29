#pragma once
#include "core/Vec3.hpp"

struct Photon {
    Vec3 pos;                 // m
    Vec3 dir{0.0, 0.0, -1.0}; // unit vector
    double wavelength_nm = 400.0;
    double time_ns = 0.0;
    double weight = 1.0;

    void normalizeDirection() {
        dir = dir.normalized();
    }
};
