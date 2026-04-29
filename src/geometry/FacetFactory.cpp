#include "geometry/FacetFactory.hpp"
#include <cmath>

std::vector<MirrorFacet> FacetFactory::buildFacets(const DishPrescription& dish,
                                                   const FacetGridConfig& grid)
{
    std::vector<MirrorFacet> facets;

    double R = dish.dish_radius;
    double step = grid.facet_spacing;
    int nmax = static_cast<int>(std::ceil(R / step));

    int id = 0;

    // 这里先用规则方格生成镜片中心，再用圆口径裁切。
    for (int iy = -nmax; iy <= nmax; ++iy) {
        for (int ix = -nmax; ix <= nmax; ++ix) {
            double x = ix * step;
            double y = iy * step;

            double r = std::sqrt(x * x + y * y);
            if (r > R) continue;

            MirrorFacet facet;
            if (dish.type == DishType::DaviesCotton) {
                facet = makeDaviesCottonFacet(id, x, y, dish, grid);
            } else {
                facet = makeParabolicFacet(id, x, y, dish, grid);
            }

            facets.push_back(facet);
            ++id;
        }
    }

    return facets;
}

Vec3 FacetFactory::designNormalToFocus(const Vec3& facet_center, const Vec3& focus)
{
    // 假设轴上平行光从 +z 往 -z 入射。
    // 从镜片看回“来光方向”是 +z；再和“去焦点方向”取角平分线，
    // 得到镜片中心理想法向。
    Vec3 to_source{0.0, 0.0, 1.0};
    Vec3 to_focus = (focus - facet_center).normalized();
    return (to_source + to_focus).normalized();
}

MirrorFacet FacetFactory::makeDaviesCottonFacet(int id, double x, double y,
                                                const DishPrescription& dish,
                                                const FacetGridConfig& grid)
{
    MirrorFacet facet;
    facet.id = id;
    facet.aperture_shape = grid.aperture_shape;
    facet.surface_type = grid.surface_type;
    facet.size1 = grid.facet_radius;

    // DC 结构：镜片中心落在球面母面上
    double R_dish = dish.dish_shape_length;
    Vec3 focus = dish.vertex + Vec3{0.0, 0.0, dish.telescope_focal_length};

    double r2 = x * x + y * y;
    double z = dish.vertex.z;

    if (r2 < R_dish * R_dish) {
        z = focus.z - std::sqrt(R_dish * R_dish - r2);
    }

    facet.center = {dish.vertex.x + x, dish.vertex.y + y, z};

    // 镜片中心法向：让中心光线反射到系统焦点
    facet.normal = designNormalToFocus(facet.center, focus);

    // DC 的标准近似：单块球面镜片焦距取系统焦距
    // 所以曲率半径 R = 2 f_tel
    facet.radius_of_curvature = 2.0 * dish.telescope_focal_length;

    return facet;
}

MirrorFacet FacetFactory::makeParabolicFacet(int id, double x, double y,
                                             const DishPrescription& dish,
                                             const FacetGridConfig& grid)
{
    MirrorFacet facet;
    facet.id = id;
    facet.aperture_shape = grid.aperture_shape;
    facet.surface_type = grid.surface_type;
    facet.size1 = grid.facet_radius;

    // 抛物面母面：z = (x^2 + y^2) / (4 f_dish)
    double f_dish = dish.dish_shape_length;
    double z = dish.vertex.z + (x * x + y * y) / (4.0 * f_dish);

    facet.center = {dish.vertex.x + x, dish.vertex.y + y, z};

    Vec3 focus = dish.vertex + Vec3{0.0, 0.0, dish.telescope_focal_length};

    // 镜片中心法向：让中心光线反射到系统焦点
    facet.normal = designNormalToFocus(facet.center, focus);

    // 抛物面母面上的 spherical facet：
    // 第一版先取镜片中心到焦点距离作为镜片近轴焦距，
    // 所以曲率半径 R = 2 * distance(center, focus)
    double fi = (focus - facet.center).norm();
    facet.radius_of_curvature = 2.0 * fi;

    return facet;
}
