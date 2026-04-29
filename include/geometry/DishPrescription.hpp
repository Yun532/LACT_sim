#pragma once
#include "core/Vec3.hpp"
#include "geometry/DishType.hpp"

// 主镜整体设计参数
// 注意：这是“整台主镜”的设计，不是单块镜片参数
struct DishPrescription {
    // 主镜母面类型：DC 或 Parabolic
    DishType type = DishType::Parabolic;

    // 整台系统的设计焦距
    double telescope_focal_length = 5.0;   // m

    // 母面参数：
    // - 对 DC：镜片中心所在球面的半径
    // - 对 Parabolic：母面抛物线的焦距
    double dish_shape_length = 5.0;        // m

    // 主镜总口径半径
    double dish_radius = 2.0;              // m

    // 主镜顶点位置
    Vec3 vertex{0.0, 0.0, 0.0};

    // 主镜光轴方向
    // 第一版先默认沿 +z
    Vec3 optical_axis{0.0, 0.0, 1.0};
};
