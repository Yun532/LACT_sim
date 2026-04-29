#pragma once
#include <optional>
#include <utility>
#include "core/Photon.hpp"
#include "core/OpticalHit.hpp"
#include "geometry/MirrorLayout.hpp"
#include "geometry/CameraGeometry.hpp"
#include "optics/OpticalEfficiency.hpp"

class RayTracer {
public:
    OpticalHit trace(const Photon& photon,
                     const MirrorLayout& mirrors,
                     const CameraGeometry& camera,
                     const OpticalEfficiency& eff) const;

private:
    struct MirrorIntersection {
        double t = 0.0;
        Vec3 point;
        Vec3 normal;
    };

    static std::optional<MirrorIntersection>
    intersectMirror(const Vec3& p0, const Vec3& d, const MirrorTile& tile);

    static std::optional<MirrorIntersection>
    intersectPlaneDisk(const Vec3& p0, const Vec3& d, const MirrorTile& tile);

    static std::optional<MirrorIntersection>
    intersectSphericalFacet(const Vec3& p0, const Vec3& d, const MirrorTile& tile);

    static std::optional<MirrorIntersection>
    intersectParaboloid(const Vec3& p0, const Vec3& d, const MirrorTile& tile);

    static std::optional<std::pair<double, Vec3>>
    intersectPlaneZ(const Vec3& p0, const Vec3& d, double zplane);
};
