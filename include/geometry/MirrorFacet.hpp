#pragma once
#include <limits>
#include <cmath>
#include "core/Vec3.hpp"
#include "geometry/ApertureShape.hpp"
#include "geometry/MirrorTile.hpp"

// 单块镜片的理想几何参数
struct MirrorFacet {
    int id = -1;

    // 镜片中心（光学表面参考点）
    Vec3 center;

    // 镜片中心处的表面法向（单位向量）
    Vec3 normal{0.0, 0.0, 1.0};

    // 镜片实际表面类型
    SurfaceType surface_type = SurfaceType::Spherical;

    // 曲率半径 R
    double radius_of_curvature = std::numeric_limits<double>::infinity();

    // 镜片边界形状
    ApertureShape aperture_shape = ApertureShape::Circular;

    // 尺寸参数
    // 对圆形镜片：size1 = 半径
    // 对六边形镜片：size1 = 相对平边距离 flat-to-flat
    // 对方形镜片：size1 = 边长
    double size1 = 0.25;   // m
    double size2 = 0.0;    // 预留
    double aperture_rotation_rad = 0.0;
    Vec3 aperture_u_axis{0.0, 0.0, 0.0};
    Vec3 aperture_v_axis{0.0, 0.0, 0.0};

    // 预留给后续真实误差模型
    double reflectivity_scale = 1.0;
    double roughness_sigma_rad = 0.0;
    double misalign_sigma_rad = 0.0;

    // 调试辅助：返回近轴焦距
    double focalLength() const {
        return std::isfinite(radius_of_curvature) ? 0.5 * radius_of_curvature
                                                  : std::numeric_limits<double>::infinity();
    }

    double apertureBoundingRadius() const {
        if (aperture_shape == ApertureShape::Hexagon) {
            return size1 / std::sqrt(3.0);
        }
        if (aperture_shape == ApertureShape::Square) {
            return size1 / std::sqrt(2.0);
        }
        return size1;
    }

    void apertureFrame(Vec3& u_axis, Vec3& v_axis) const {
        if (aperture_u_axis.norm2() > 0.0 && aperture_v_axis.norm2() > 0.0) {
            u_axis = aperture_u_axis.normalized();
            v_axis = aperture_v_axis.normalized();
            return;
        }
        Vec3 n = normal.normalized();
        Vec3 ref = (std::abs(n.z) < 0.9) ? Vec3{0.0, 0.0, 1.0}
                                         : Vec3{0.0, 1.0, 0.0};
        u_axis = ref.cross(n).normalized();
        v_axis = n.cross(u_axis).normalized();
    }

    // 转成当前追迹器兼容的 MirrorTile
    MirrorTile toMirrorTile() const {
        MirrorTile t;
        t.id = id;
        t.center = center;
        t.normal = normal;
        t.aperture_radius = apertureBoundingRadius();
        t.aperture_shape = aperture_shape;
        t.aperture_size1 = size1;
        t.aperture_size2 = size2;
        t.aperture_rotation_rad = aperture_rotation_rad;
        apertureFrame(t.aperture_u_axis, t.aperture_v_axis);
        t.radius_of_curvature = radius_of_curvature;
        t.reflectivity_scale = reflectivity_scale;
        t.roughness_sigma_rad = roughness_sigma_rad;
        t.misalign_sigma_rad = misalign_sigma_rad;
        t.type = surface_type;
        return t;
    }
};
