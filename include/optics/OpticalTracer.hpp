#pragma once
#include <cstdint>
#include <optional>
#include <utility>
#include "core/Photon.hpp"
#include "geometry/MirrorLayout.hpp"
#include "optics/OutputPlane.hpp"
#include "optics/OpticalSurfaceHit.hpp"
#include "optics/OpticalEfficiency.hpp"

class OpticalTracer {
public:
    OpticalTracer(double speed_of_light_m_per_ns = 0.299792458,
                  double reflect_direction_sigma_rad = 0.0,
                  std::uint64_t random_seed = 123456789ULL);

    // 只追到某个输出平面，不做像素映射
    OpticalSurfaceHit traceToPlane(const Photon& photon,
                                   const MirrorLayout& mirrors,
                                   const OutputPlane& plane,
                                   const OpticalEfficiency& eff) const;

    // 用于处理平移后的 2D EventIO 记录平面。
    // 镜面交点在完整直线 p + t*dir 上寻找，因此 t 可以为正
    // （从记录平面继续向前传播）或为负（记录平面已经在镜面之后）。
    // 反射仍然使用原始入射方向 photon.dir，时间计算中保留 t 的符号。
    OpticalSurfaceHit traceBackprojectedToPlane(const Photon& photon,
                                                const MirrorLayout& mirrors,
                                                const OutputPlane& plane,
                                                const OpticalEfficiency& eff) const;

private:
    struct MirrorIntersection {
        double t = 0.0;
        Vec3 point;
        Vec3 normal;
    };

    static std::optional<MirrorIntersection>
    intersectMirror(const Vec3& p0, const Vec3& d, const MirrorTile& tile,
                    bool allow_negative_t = false);

    static std::optional<MirrorIntersection>
    intersectPlaneDisk(const Vec3& p0, const Vec3& d, const MirrorTile& tile,
                       bool allow_negative_t = false);

    static std::optional<MirrorIntersection>
    intersectSphericalFacet(const Vec3& p0, const Vec3& d, const MirrorTile& tile,
                            bool allow_negative_t = false);

    static std::optional<MirrorIntersection>
    intersectParaboloid(const Vec3& p0, const Vec3& d, const MirrorTile& tile,
                        bool allow_negative_t = false);

    static std::optional<std::pair<double, Vec3>>
    intersectOutputPlane(const Vec3& p0, const Vec3& d, const OutputPlane& plane);

    double speed_of_light_m_per_ns_ = 0.299792458;
    double reflect_direction_sigma_rad_ = 0.0;
    std::uint64_t random_seed_ = 123456789ULL;
};
