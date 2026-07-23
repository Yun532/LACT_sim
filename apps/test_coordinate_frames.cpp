#include "app/OpticalSimCommon.hpp"

#include <cmath>
#include <iostream>

using namespace lact;

namespace {

bool near(double a, double b, double eps = 1e-12)
{
    return std::abs(a - b) <= eps;
}

bool nearVec(const Vec3& a, const Vec3& b, double eps = 1e-12)
{
    return near(a.x, b.x, eps) && near(a.y, b.y, eps) && near(a.z, b.z, eps);
}

bool check(bool condition, const char* message)
{
    if (!condition) {
        std::cerr << message << '\n';
        return false;
    }
    return true;
}

} // namespace

int main()
{
    bool ok = true;
    TelescopeConfig telescope;
    telescope.position_m = {100.0, 200.0, 4.4};
    telescope.pointing_az_deg = 0.0;
    telescope.pointing_el_deg = 70.0;

    const double el = 70.0 * DEG_TO_RAD;
    const TelescopeFrame frame = buildCorsikaNwuTelescopeFrame(telescope);
    ok &= check(nearVec(frame.x_axis, {-std::sin(el), 0.0, std::cos(el)}),
                "CORSIKA local +x must point toward increasing elevation");
    ok &= check(nearVec(frame.y_axis, {0.0, -1.0, 0.0}),
                "CORSIKA local +y must point toward increasing azimuth/East at az=0");
    ok &= check(nearVec(frame.z_axis, {std::cos(el), 0.0, std::sin(el)}),
                "CORSIKA local +z must be the boresight");
    ok &= check(nearVec(frame.x_axis.cross(frame.y_axis), frame.z_axis),
                "CORSIKA optical frame must be right-handed");

    PhotonBunch local;
    local.photon.pos = {1.0, -2.0, 3.0};
    local.photon.dir = Vec3{0.02, -0.01, -1.0}.normalized();
    const Vec3 nwu_relative_pos = frame.rotateVector(local.photon.pos);
    const Vec3 nwu_direction = frame.rotateVector(local.photon.dir).normalized();

    PhotonBunch relative = local;
    relative.photon.pos = nwu_relative_pos;
    relative.photon.dir = nwu_direction;
    const PhotonBunch relative_local = transformBunchToTelescopeLocal(
        relative, telescope, "corsika_nwu_relative");
    ok &= check(nearVec(relative_local.photon.pos, local.photon.pos),
                "relative NWU positions must rotate without subtracting telescope position");
    ok &= check(nearVec(relative_local.photon.dir, local.photon.dir),
                "relative NWU direction round trip failed");

    PhotonBunch global = relative;
    global.photon.pos = telescope.position_m + nwu_relative_pos;
    const PhotonBunch global_local = transformBunchToTelescopeLocal(
        global, telescope, "corsika_nwu_global");
    ok &= check(nearVec(global_local.photon.pos, local.photon.pos),
                "global NWU positions must subtract telescope.position_m exactly once");

    // ENU east-start is the same physical frame after:
    //   East = -West, North = North, az_ENU = 90 deg - az_NWU.
    TelescopeConfig enu_telescope = telescope;
    enu_telescope.position_m = {-telescope.position_m.y,
                                 telescope.position_m.x,
                                 telescope.position_m.z};
    enu_telescope.pointing_az_deg = 90.0 - telescope.pointing_az_deg;
    const TelescopeFrame enu_frame = buildEnuEastTelescopeFrame(enu_telescope);
    const Vec3 enu_relative_pos{-nwu_relative_pos.y,
                                nwu_relative_pos.x,
                                nwu_relative_pos.z};
    const Vec3 enu_direction{-nwu_direction.y,
                             nwu_direction.x,
                             nwu_direction.z};
    ok &= check(nearVec(enu_frame.x_axis,
                        {-frame.x_axis.y, frame.x_axis.x, frame.x_axis.z}),
                "ENU and NWU elevation axes must describe the same physical direction");
    ok &= check(nearVec(enu_frame.y_axis,
                        {-frame.y_axis.y, frame.y_axis.x, frame.y_axis.z}),
                "ENU and NWU camera-y axes must describe the same physical direction");
    ok &= check(nearVec(enu_frame.z_axis,
                        {-frame.z_axis.y, frame.z_axis.x, frame.z_axis.z}),
                "ENU and NWU boresights must describe the same physical direction");

    PhotonBunch enu_relative = local;
    enu_relative.photon.pos = enu_relative_pos;
    enu_relative.photon.dir = enu_direction;
    const PhotonBunch enu_relative_local = transformBunchToTelescopeLocal(
        enu_relative, enu_telescope, "enu_east_relative");
    ok &= check(nearVec(enu_relative_local.photon.pos, local.photon.pos),
                "relative ENU positions must rotate without subtracting telescope position");
    ok &= check(nearVec(enu_relative_local.photon.dir, local.photon.dir),
                "relative ENU direction round trip failed");

    PhotonBunch enu_global = enu_relative;
    enu_global.photon.pos = enu_telescope.position_m + enu_relative_pos;
    const PhotonBunch enu_global_local = transformBunchToTelescopeLocal(
        enu_global, enu_telescope, "enu_east_global");
    ok &= check(nearVec(enu_global_local.photon.pos, local.photon.pos),
                "global ENU positions must subtract telescope.position_m exactly once");

    TelescopeConfig east_pointing;
    east_pointing.pointing_az_deg = 0.0;
    east_pointing.pointing_el_deg = 0.0;
    const TelescopeFrame east_frame = buildEnuEastTelescopeFrame(east_pointing);
    ok &= check(nearVec(east_frame.z_axis, {1.0, 0.0, 0.0}),
                "ENU east-start az=0 must point East");
    east_pointing.pointing_az_deg = 90.0;
    const TelescopeFrame north_frame = buildEnuEastTelescopeFrame(east_pointing);
    ok &= check(nearVec(north_frame.z_axis, {0.0, 1.0, 0.0}),
                "ENU east-start az=90 must point North");

    const TelescopeFrame existing_trace_frame = buildTelescopeFrame(telescope);
    const Vec3 reconstructed_world = sourceDirectionInWorld(
        local, telescope, "telescope_local");
    ok &= check(nearVec(reconstructed_world,
                        existing_trace_frame.rotateVector(local.photon.dir).normalized()),
                "local atmosphere direction must preserve the existing telescope frame");
    ok &= check(near(sourceDirectionInWorld(relative, telescope,
                                             "corsika_nwu_relative").z,
                     nwu_direction.z),
                "NWU atmosphere direction must be used directly");

    // A horizontally symmetric atmosphere must depend on zenith angle, not on
    // camera azimuth around the boresight. Convert two NWU directions with the
    // same vertical component through local coordinates and recover that z.
    const double horizontal = std::sqrt(1.0 - 0.75 * 0.75);
    for (const Vec3& world_dir : {Vec3{horizontal, 0.0, -0.75},
                                  Vec3{0.0, horizontal, -0.75}}) {
        PhotonBunch local_direction;
        local_direction.photon.dir =
            existing_trace_frame.rotateVectorToLocal(world_dir).normalized();
        ok &= check(near(sourceDirectionInWorld(local_direction, telescope,
                                                "telescope_local").z,
                         -0.75),
                    "atmosphere zenith cosine must be invariant under world azimuth");
    }

    const SourceRuntimeConfig default_csv = buildSourceRuntimeConfig({
        {"source.mode", "PhotonCsv"},
        {"source.csv_path", "photons.csv"},
    });
    ok &= check(default_csv.coordinate_frame == "telescope_local",
                "PhotonCsv must default to the existing telescope-local frame");
    const TelescopeConfig default_telescope = buildTelescopeConfig({});
    ok &= check(nearVec(default_telescope.position_m, {0.0, 0.0, 0.0}),
                "unset telescope.position_m must default to the array origin");

    ok &= check(normalizeSourceCoordinateFrame("corsika_iact") ==
                    "corsika_nwu_relative",
                "legacy corsika_iact alias must remain supported");
    ok &= check(normalizeSourceCoordinateFrame("local") == "telescope_local",
                "legacy local alias must remain supported");
    ok &= check(normalizeSourceCoordinateFrame("enu_relative") ==
                    "enu_east_relative",
                "ENU relative alias must normalize to the explicit east-start name");
    ok &= check(normalizeSourceCoordinateFrame("enu_global") ==
                    "enu_east_global",
                "ENU global alias must normalize to the explicit east-start name");

    return ok ? 0 : 1;
}
