#pragma once
#include <cmath>
#include "core/Vec3.hpp"

// 输出参考平面
// 光学追迹的结果先打到这个平面上，不引入像素几何。
// 这样后面可以独立保存完整光子文件，再做相机映射。
struct OutputPlane {
    // 平面上一点
    Vec3 point{0.0, 0.0, 0.0};

    // 平面法向（单位向量）
    Vec3 normal{0.0, 0.0, 1.0};

    // 平面内的局部坐标基底
    // 用于把三维命中点投影成二维 (u, v)，方便画 spot 图
    Vec3 u_axis{1.0, 0.0, 0.0};
    Vec3 v_axis{0.0, 1.0, 0.0};

    // 根据 normal 自动构造平面内两个正交基底
    void buildLocalFrame() {
        normal = normal.normalized();

        // 选一个不和 normal 平行的参考向量
        Vec3 ref = (std::abs(normal.z) < 0.9) ? Vec3{0.0, 0.0, 1.0}
                                              : Vec3{0.0, 1.0, 0.0};

        // 构造平面内正交基
        u_axis = ref.cross(normal).normalized();
        v_axis = normal.cross(u_axis).normalized();
    }
};
