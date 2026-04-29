#pragma once
#include <vector>
#include "geometry/MirrorFacet.hpp"
#include "geometry/DishPrescription.hpp"

// 用于生成镜片列表的简化配置
struct FacetGridConfig {
    // 镜片中心的近似间距
    double facet_spacing = 0.55;   // m

    // 单块镜片边界尺寸
    double facet_radius = 0.22;    // m

    // 镜片边界形状
    ApertureShape aperture_shape = ApertureShape::Circular;

    // 镜片表面类型
    SurfaceType surface_type = SurfaceType::Spherical;
};

// 负责根据主镜整体设计，生成理想镜片参数
class FacetFactory {
public:
    static std::vector<MirrorFacet> buildFacets(const DishPrescription& dish,
                                                const FacetGridConfig& grid);

private:
    static MirrorFacet makeDaviesCottonFacet(int id, double x, double y,
                                             const DishPrescription& dish,
                                             const FacetGridConfig& grid);

    static MirrorFacet makeParabolicFacet(int id, double x, double y,
                                          const DishPrescription& dish,
                                          const FacetGridConfig& grid);

    static Vec3 designNormalToFocus(const Vec3& facet_center, const Vec3& focus);
};
