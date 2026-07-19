#pragma once
#include <cstdint>
#include "core/Vec3.hpp"

struct Photon {
    Vec3 pos;                 // m
    Vec3 dir{0.0, 0.0, -1.0}; // unit vector
    double wavelength_nm = 400.0;
    double time_ns = 0.0;
    double weight = 1.0;
    // Stable per-represented-photon identity for stochastic optical effects.
    std::uint64_t random_stream_id = 0;
    // Negative-wavelength EventIO bunches are CEFFIC photoelectron bunches:
    // their wavelength-dependent detector efficiency was already applied.
    bool optical_efficiency_preapplied = false;

    void normalizeDirection() {
        dir = dir.normalized();
    }
};
