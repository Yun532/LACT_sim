#pragma once
#include "core/Vec3.hpp"

struct OpticalHit {
    bool hit_mirror = false;
    bool hit_camera = false;
    bool accepted = false;

    int mirror_id = -1;
    int pixel_id = -1;

    Vec3 mirror_point;
    Vec3 camera_point;

    double xcam = 0.0;
    double ycam = 0.0;
    double sxcam = 0.0;
    double sycam = 0.0;

    double time_to_camera_ns = 0.0;
    double relative_efficiency = 1.0;
};
