#pragma once
#include <limits>
#include "geometry/ApertureShape.hpp"
#include "core/Vec3.hpp"

enum class SurfaceType {
    Spherical,
    Parabolic,
    Polynomial,
    Planar
};

struct MirrorTile {
    int id = -1;

    // 镜片中心位置
    Vec3 center;

    // 镜片中心处法向
    Vec3 normal{0.0, 0.0, 1.0};

    // 单块镜片有效半径（目前仍按圆形镜片处理）
    double aperture_radius = 0.5;   // m, bounding radius for quick rejection
    ApertureShape aperture_shape = ApertureShape::Circular;
    double aperture_size1 = 0.5;    // circular: radius; hexagon: flat-to-flat; square: side
    double aperture_size2 = 0.0;
    double aperture_rotation_rad = 0.0;
    Vec3 aperture_u_axis{0.0, 0.0, 0.0};
    Vec3 aperture_v_axis{0.0, 0.0, 0.0};

    // 曲率半径 R
    // 对球面镜片，近轴焦距 f = R / 2
    // 对平面镜可设为 inf
    double radius_of_curvature = std::numeric_limits<double>::infinity();

    double reflectivity_scale = 1.0;
    double roughness_sigma_rad = 0.0;
    double misalign_sigma_rad = 0.0;

    SurfaceType type = SurfaceType::Planar;

    // 仅作为调试/打印辅助，不建议在光追核心里依赖它
    double focal_length() const {
        return std::isfinite(radius_of_curvature) ? 0.5 * radius_of_curvature
                                                  : std::numeric_limits<double>::infinity();
    }
};
