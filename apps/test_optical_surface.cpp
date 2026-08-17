#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include "io/SyntheticPhotonSource.hpp"
#include "io/SurfaceHitCsvWriter.hpp"

#include "geometry/DishPrescription.hpp"
#include "geometry/FacetFactory.hpp"
#include "geometry/FacetLayoutUtils.hpp"

#include "optics/OutputPlane.hpp"
#include "optics/OpticalSurfaceHit.hpp"
#include "optics/OpticalTracer.hpp"
#include "optics/OpticalEfficiency.hpp"

int main() {
    // =========================================
    // 1) 构造主镜整体设计
    // =========================================
    DishPrescription dish;
    dish.type = DishType::DaviesCotton;       // 可改成 DishType::DaviesCotton
    dish.telescope_focal_length = 5.0;     // 系统设计焦距
    dish.dish_shape_length = 5.0;          // 母面参数
    dish.dish_radius = 2.0;                // 主镜半径
    dish.vertex = {0.0, 0.0, 0.0};
    dish.optical_axis = {0.0, 0.0, 1.0};

    // =========================================
    // 2) 构造镜片网格配置
    // =========================================
    FacetGridConfig grid;
    grid.facet_spacing = 0.55;
    grid.facet_radius = 0.22;
    grid.aperture_shape = ApertureShape::Circular;
    grid.surface_type = SurfaceType::Spherical;

    // 生成镜片列表
    auto facets = FacetFactory::buildFacets(dish, grid);

    // 转成当前追迹器兼容的 MirrorLayout
    MirrorLayout mirrors = makeMirrorLayoutFromFacets(facets);

    // =========================================
    // 3) 设置输出参考平面
    //    这里先取 z = 设计焦距 的平面
    // =========================================
    OutputPlane plane;
    plane.point = {0.0, 0.0, dish.vertex.z + dish.telescope_focal_length};
    plane.normal = {0.0, 0.0, 1.0};
    plane.buildLocalFrame();

    // =========================================
    // 4) 构造输入光子
    // =========================================
    SyntheticPhotonConfig src_cfg;
    src_cfg.mode = SyntheticMode::ParallelBeam;
    src_cfg.n_bunches = 200000;
    src_cfg.multiplicity = 1.0;
    src_cfg.wavelength_nm = 400.0;
    src_cfg.time_ns = 0.0;
    src_cfg.source_plane_z = 30.0;
    src_cfg.beam_radius_m = 1.5;
    src_cfg.beam_direction = {0.0, 0.0, -1.0};

    SyntheticPhotonSource source(src_cfg);

    // =========================================
    // 5) 光学追迹
    // =========================================
    OpticalEfficiency eff;
    OpticalTracer tracer;

    std::vector<OpticalSurfaceHit> hits;
    hits.reserve(src_cfg.n_bunches);

    PhotonBunch bunch;
    int n_total = 0;
    int n_hit_mirror = 0;
    int n_hit_surface = 0;

    double sum_w = 0.0;
    double sum_r2 = 0.0;

    while (source.next(bunch)) {
        ++n_total;

        Photon p = bunch.photon;
        p.normalizeDirection();

        // 把 bunch 的 multiplicity 合并进 photon.weight
        p.weight *= bunch.multiplicity;

        OpticalSurfaceHit hit = tracer.traceToPlane(p, mirrors, plane, eff);

        if (hit.hit_mirror) ++n_hit_mirror;
        if (hit.hit_surface) {
            ++n_hit_surface;
            hits.push_back(hit);

            double w = hit.weight * hit.relative_efficiency;
            double r2 = hit.u_m * hit.u_m + hit.v_m * hit.v_m;

            sum_w += w;
            sum_r2 += w * r2;
        }
    }

    // =========================================
    // 6) 保存结果
    // =========================================
    const std::string out_csv = "surface_hits.csv";
    bool ok = writeSurfaceHitsCSV(out_csv, hits);

    // =========================================
    // 7) 打印统计
    // =========================================
    std::cout << std::fixed << std::setprecision(6);

    std::cout << "Number of facets      = " << mirrors.size() << "\n";
    std::cout << "Total photons         = " << n_total << "\n";
    std::cout << "Hit mirror            = " << n_hit_mirror << "\n";
    std::cout << "Hit output plane      = " << n_hit_surface << "\n";

    if (sum_w > 0.0) {
        double rms = std::sqrt(sum_r2 / sum_w);
        std::cout << "Weighted spot RMS [m] = " << rms << "\n";
    } else {
        std::cout << "Weighted spot RMS [m] = nan\n";
    }

    std::cout << "CSV written           = " << (ok ? "yes" : "no") << "\n";
    std::cout << "Output file           = " << out_csv << "\n";

    // 只示范打印前几条
    int nshow = std::min<int>(5, hits.size());
    for (int i = 0; i < nshow; ++i) {
        const auto& h = hits[i];
        std::cout
            << "hit[" << i << "]"
            << " mirror_id=" << h.mirror_id
            << " u=" << h.u_m
            << " v=" << h.v_m
            << " time_ns=" << h.time_ns
            << " eff=" << h.relative_efficiency
            << "\n";
    }

    return 0;
}
