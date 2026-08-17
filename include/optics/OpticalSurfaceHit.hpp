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

    // 可选：进入光追前的望远镜本地输入光子。该快照在将光子放入
    // 全局镜面/焦平面几何之前保存，因此不随原始 source.coordinate_frame 改变。
    bool has_input_photon = false;
    Vec3 input_pos_local;
    Vec3 input_dir_local;

    // 在输出平面局部坐标系中的二维坐标
    double u_m = 0.0;
    double v_m = 0.0;

    // 遮挡诊断标记：
    // 默认正式模拟会丢弃被遮挡光子；只有 obstruction.mark_only=true 时，
    // 被遮挡光子仍会写出，用这些字段在图里单独上色。
    bool obstruction_blocked = false;
    bool obstruction_blocked_incoming = false;
    bool obstruction_blocked_reflected = false;

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
    double collector_path_length_m = 0.0;
    double collector_time_delay_ns = 0.0;
    bool collector_reflection_limit_reached = false;
    double collector_exit_x_m = 0.0;
    double collector_exit_y_m = 0.0;
    double collector_exit_z_m = 0.0;
    double collector_dir_u = 0.0;
    double collector_dir_v = 0.0;
    double collector_dir_w = 0.0;

    // Optional S17351 channel geometry diagnostics, populated before PDE.
    bool sipm_geometry_enabled = false;
    bool sipm_gap_rejected = false;
    bool sipm_channel_gap_rejected = false;
    int sipm_grid_column = -1;
    int sipm_grid_row = -1;
    int sipm_channel_id = -1;
    int sipm_microcell_id = -1;

    // 时间、波长、权重
    double time_ns = 0.0;
    double wavelength_nm = 0.0;
    double weight = 1.0;

    // 到该平面为止的相对光学效率
    double relative_efficiency = 1.0;
};
