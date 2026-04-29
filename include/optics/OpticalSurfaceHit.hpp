#pragma once
#include "core/Vec3.hpp"

// 光学部分的输出：
// 只记录光子打到“输出参考平面”时的信息，不涉及任何像素布局。
struct OpticalSurfaceHit {
    bool hit_mirror = false;
    bool hit_surface = false;

    int mirror_id = -1;

    // 镜面命中点
    Vec3 mirror_point;

    // 输出平面命中点（全局三维坐标）
    Vec3 surface_point;

    // 光子打到输出平面时的传播方向
    Vec3 out_dir;

    // 在输出平面局部坐标系中的二维坐标
    double u_m = 0.0;
    double v_m = 0.0;

    // 可选相机几何映射：
    // camera_enabled=false 表示这次运行只使用白板输出。
    bool camera_enabled = false;
    bool hit_camera = false;
    bool accepted = false;
    int pixel_id = -1;
    double camera_x_m = 0.0;
    double camera_y_m = 0.0;

    // 可选光收集器诊断量：
    // 坐标为单个 pixel/collector 的局部坐标，单位 m。
    bool collector_enabled = false;
    bool hit_collector = false;
    int collector_reflections = 0;
    double collector_intensity = 1.0;
    double collector_exit_x_m = 0.0;
    double collector_exit_y_m = 0.0;
    double collector_exit_z_m = 0.0;
    double collector_dir_u = 0.0;
    double collector_dir_v = 0.0;
    double collector_dir_w = 0.0;

    // 时间、波长、权重
    double time_ns = 0.0;
    double wavelength_nm = 0.0;
    double weight = 1.0;

    // 到该平面为止的相对光学效率
    double relative_efficiency = 1.0;
};
