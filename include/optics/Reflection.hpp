#pragma once
#include "core/Vec3.hpp"

inline Vec3 reflectDirection(const Vec3& in_dir, const Vec3& normal_unit) {
    return (in_dir - normal_unit * (2.0 * in_dir.dot(normal_unit))).normalized();
}
