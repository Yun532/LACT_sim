#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include "core/Photon.hpp"

#include "geometry/DishPrescription.hpp"
#include "geometry/FacetFactory.hpp"
#include "geometry/FacetLayoutUtils.hpp"

#include "optics/OutputPlane.hpp"
#include "optics/OpticalTracer.hpp"
#include "optics/OpticalEfficiency.hpp"

// ============================================================
// 这个程序的目的：
// 对当前 DC 主镜的“每一块镜片中心”发射一条轴上平行入射光，
// 看它在输出平面上最终落到哪里。
// 
// 光子构造方式：
//   对于每个 facet：
//   - 光子初始 x,y 取 facet.center.x, facet.center.y
//   - 光子从上方 z = z_source 平行向下入射
//   - dir = (0,0,-1)
//
// 这样这条光线会穿过该镜片中心所在的竖直线，
// 用来检查“镜片中心光线是否被正确反射到焦面附近”。
// ============================================================

int main() {
    // ========================================================
    // 1) 构造 DC 主镜设计
    // ========================================================
    DishPrescription dish;
    dish.type = DishType::DaviesCotton;
    dish.telescope_focal_length = 5.0;   // 系统焦距
    dish.dish_shape_length = 5.0;        // DC 母面球半径，先取 = f
    dish.dish_radius = 2.0;              // 主镜半径
    dish.vertex = {0.0, 0.0, 0.0};
    dish.optical_axis = {0.0, 0.0, 1.0};

    // ========================================================
    // 2) 构造镜片网格
    // ========================================================
    FacetGridConfig grid;
    grid.facet_spacing = 0.55;
    grid.facet_radius = 0.22;
    grid.aperture_shape = ApertureShape::Circular;
    grid.surface_type = SurfaceType::Spherical;

    auto facets = FacetFactory::buildFacets(dish, grid);
    MirrorLayout mirrors = makeMirrorLayoutFromFacets(facets);

    // ========================================================
    // 3) 定义输出平面
    //    这里先放在 z = f_tel
    // ========================================================
    OutputPlane plane;
    plane.point = {0.0, 0.0, dish.vertex.z + dish.telescope_focal_length};
    plane.normal = {0.0, 0.0, 1.0};
    plane.buildLocalFrame();

    // ========================================================
    // 4) 光学追迹器
    // ========================================================
    OpticalTracer tracer;
    OpticalEfficiency eff;

    // ========================================================
    // 5) 对每个镜片中心打一条平行光
    // ========================================================
    const double z_source = 30.0;  // 光子起始平面高度
    const Vec3 dir_in{0.0, 0.0, -1.0};

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Number of facets = " << facets.size() << "\n\n";

    // 统计量
    int n_total = 0;
    int n_hit_mirror = 0;
    int n_hit_surface = 0;
    int n_failed = 0;

    double sum_r2 = 0.0;
    double max_r = 0.0;
    const double focus_tolerance_m = 1e-9;

    for (const auto& facet : facets) {
        ++n_total;

        Photon p;
        // 关键点：x,y 取镜片中心，z 放在高处，方向平行向下
        p.pos = {facet.center.x, facet.center.y, z_source};
        p.dir = dir_in;
        p.wavelength_nm = 400.0;
        p.time_ns = 0.0;
        p.weight = 1.0;
        p.normalizeDirection();

        auto hit = tracer.traceToPlane(p, mirrors, plane, eff);

        if (hit.hit_mirror) ++n_hit_mirror;
        if (hit.hit_surface) {
            ++n_hit_surface;
            double r = std::sqrt(hit.u_m * hit.u_m + hit.v_m * hit.v_m);
            sum_r2 += r * r;
            if (r > max_r) max_r = r;
        }

        double r = hit.hit_surface ? std::sqrt(hit.u_m * hit.u_m + hit.v_m * hit.v_m) : -1.0;
        bool ray_ok = hit.hit_mirror
                   && hit.hit_surface
                   && hit.mirror_id == facet.id
                   && r >= 0.0
                   && r <= focus_tolerance_m;
        if (!ray_ok) {
            ++n_failed;
        }

        std::cout
            << "facet_id=" << facet.id
            << "  facet_center=("
            << facet.center.x << ", "
            << facet.center.y << ", "
            << facet.center.z << ")"
            << "  normal=("
            << facet.normal.x << ", "
            << facet.normal.y << ", "
            << facet.normal.z << ")"
            << "  Rcurv=" << facet.radius_of_curvature
            << "  mirror_id=" << hit.mirror_id
            << "  mirror_point=("
            << hit.mirror_point.x << ", "
            << hit.mirror_point.y << ", "
            << hit.mirror_point.z << ")"
            << "  out_dir=("
            << hit.out_dir.x << ", "
            << hit.out_dir.y << ", "
            << hit.out_dir.z << ")"
            << "  hit_mirror=" << hit.hit_mirror
            << "  hit_surface=" << hit.hit_surface
            << "  u=" << hit.u_m
            << "  v=" << hit.v_m
            << "  r=" << r
            << "  ok=" << ray_ok
            << "\n";
    }

    std::cout << "\n================ Summary ================\n";
    std::cout << "Total center rays   = " << n_total << "\n";
    std::cout << "Hit mirror          = " << n_hit_mirror << "\n";
    std::cout << "Hit output plane    = " << n_hit_surface << "\n";
    std::cout << "Failed checks       = " << n_failed << "\n";

    if (n_hit_surface > 0) {
        double rms = std::sqrt(sum_r2 / n_hit_surface);
        std::cout << "Center-ray RMS [m]  = " << rms << "\n";
        std::cout << "Center-ray max r[m] = " << max_r << "\n";
    } else {
        std::cout << "Center-ray RMS [m]  = nan\n";
        std::cout << "Center-ray max r[m] = nan\n";
    }

    return n_failed == 0 ? 0 : 1;
}
