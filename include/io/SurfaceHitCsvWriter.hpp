#pragma once
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include "optics/OpticalSurfaceHit.hpp"

// 把输出平面上的光子结果保存成 CSV
// 这样 Python 可以直接读，不需要先做 pybind11。
inline bool writeSurfaceHitsCSV(const std::string& path,
                                const std::vector<OpticalSurfaceHit>& hits,
                                bool include_input_photon = false)
{
    const std::filesystem::path out_path(path);
    if (out_path.has_parent_path()) {
        std::filesystem::create_directories(out_path.parent_path());
    }
    std::ofstream ofs(path);
    if (!ofs) return false;

    ofs << std::setprecision(10);

    ofs << "hit_mirror,hit_surface,mirror_id,"
        << "mirror_x,mirror_y,mirror_z,"
        << "surface_x,surface_y,surface_z,"
        << "dir_x,dir_y,dir_z,"
        << "u_m,v_m,"
        << "obstruction_blocked,obstruction_blocked_incoming,obstruction_blocked_reflected,"
        << "camera_enabled,hit_camera,accepted,pixel_id,camera_x_m,camera_y_m,"
        << "collector_enabled,hit_collector,collector_reflections,collector_intensity,"
        << "collector_path_length_m,collector_time_delay_ns,"
        << "collector_reflection_limit_reached,"
        << "collector_exit_x_m,collector_exit_y_m,collector_exit_z_m,"
        << "collector_dir_u,collector_dir_v,collector_dir_w,"
        << "time_ns,wavelength_nm,weight,relative_efficiency";
    if (include_input_photon) {
        ofs << ",input_local_x_m,input_local_y_m,input_local_z_m,"
            << "input_local_dir_x,input_local_dir_y,input_local_dir_z";
    }
    ofs << "\n";

    for (const auto& h : hits) {
        ofs << (h.hit_mirror ? 1 : 0) << ","
            << (h.hit_surface ? 1 : 0) << ","
            << h.mirror_id << ","
            << h.mirror_point.x << "," << h.mirror_point.y << "," << h.mirror_point.z << ","
            << h.surface_point.x << "," << h.surface_point.y << "," << h.surface_point.z << ","
            << h.out_dir.x << "," << h.out_dir.y << "," << h.out_dir.z << ","
            << h.u_m << "," << h.v_m << ","
            << (h.obstruction_blocked ? 1 : 0) << ","
            << (h.obstruction_blocked_incoming ? 1 : 0) << ","
            << (h.obstruction_blocked_reflected ? 1 : 0) << ","
            << (h.camera_enabled ? 1 : 0) << ","
            << (h.hit_camera ? 1 : 0) << ","
            << (h.accepted ? 1 : 0) << ","
            << h.pixel_id << ","
            << h.camera_x_m << "," << h.camera_y_m << ","
            << (h.collector_enabled ? 1 : 0) << ","
            << (h.hit_collector ? 1 : 0) << ","
            << h.collector_reflections << ","
            << h.collector_intensity << ","
            << h.collector_path_length_m << ","
            << h.collector_time_delay_ns << ","
            << (h.collector_reflection_limit_reached ? 1 : 0) << ","
            << h.collector_exit_x_m << "," << h.collector_exit_y_m << ","
            << h.collector_exit_z_m << ","
            << h.collector_dir_u << "," << h.collector_dir_v << ","
            << h.collector_dir_w << ","
            << h.time_ns << ","
            << h.wavelength_nm << ","
            << h.weight << ","
            << h.relative_efficiency;
        if (include_input_photon) {
            ofs << "," << h.input_pos_local.x
                << "," << h.input_pos_local.y
                << "," << h.input_pos_local.z
                << "," << h.input_dir_local.x
                << "," << h.input_dir_local.y
                << "," << h.input_dir_local.z;
        }
        ofs << "\n";
    }

    return true;
}
