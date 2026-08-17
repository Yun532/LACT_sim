#include "app/ArrayTimingCorrection.hpp"

#include <cmath>
#include <string>
#include <stdexcept>

namespace lact {
namespace {

constexpr double kDegToRad = 3.14159265358979323846 / 180.0;

void requireFinite(double value, const char* name)
{
    if (!std::isfinite(value)) {
        throw std::runtime_error(std::string(name) + " must be finite");
    }
}

} // namespace

double airRefractiveIndexAtAltitudeKm(double altitude_km)
{
    requireFinite(altitude_km, "array timing observation altitude");
    return 1.0 + 0.0002814 * std::exp(
        -0.0947982 * altitude_km -
        0.00134614 * altitude_km * altitude_km);
}

double airLightSpeedAtAltitudeM(double altitude_m)
{
    requireFinite(altitude_m, "array timing observation altitude");
    constexpr double vacuum_speed_m_per_ns = 0.299792458;
    return vacuum_speed_m_per_ns /
           airRefractiveIndexAtAltitudeKm(altitude_m * 1.0e-3);
}

Vec3 corsikaNwuViewingDirection(double azimuth_north_to_east_deg,
                                double altitude_deg)
{
    requireFinite(azimuth_north_to_east_deg, "array timing azimuth");
    requireFinite(altitude_deg, "array timing altitude");
    if (altitude_deg < -90.0 || altitude_deg > 90.0) {
        throw std::runtime_error(
            "array timing altitude must be within [-90, 90] degrees");
    }

    const double azimuth = azimuth_north_to_east_deg * kDegToRad;
    const double altitude = altitude_deg * kDegToRad;
    const double cos_altitude = std::cos(altitude);
    return Vec3{
        cos_altitude * std::cos(azimuth),
        -cos_altitude * std::sin(azimuth),
        std::sin(altitude),
    }.normalized();
}

double planeWavefrontGeometricDelayNs(
    const Vec3& telescope_position_nwu_m,
    const Vec3& viewing_direction_nwu,
    double wavefront_speed_m_per_ns)
{
    requireFinite(telescope_position_nwu_m.x, "array timing telescope x");
    requireFinite(telescope_position_nwu_m.y, "array timing telescope y");
    requireFinite(telescope_position_nwu_m.z, "array timing telescope z");
    requireFinite(wavefront_speed_m_per_ns, "array wavefront speed");
    if (wavefront_speed_m_per_ns <= 0.0) {
        throw std::runtime_error("array wavefront speed must be > 0");
    }
    if (!std::isfinite(viewing_direction_nwu.x) ||
        !std::isfinite(viewing_direction_nwu.y) ||
        !std::isfinite(viewing_direction_nwu.z) ||
        viewing_direction_nwu.norm() <= 0.0) {
        throw std::runtime_error(
            "array timing viewing direction must be finite and non-zero");
    }

    const Vec3 direction = viewing_direction_nwu.normalized();
    return telescope_position_nwu_m.dot(direction) /
           wavefront_speed_m_per_ns;
}

ArrayTimingCorrectionResult correctArrayTriggerTimePlaneWave(
    double raw_trigger_time_ns,
    const Vec3& telescope_position_nwu_m,
    double azimuth_north_to_east_deg,
    double altitude_deg,
    double wavefront_speed_m_per_ns)
{
    requireFinite(raw_trigger_time_ns, "raw array trigger time");
    const Vec3 direction = corsikaNwuViewingDirection(
        azimuth_north_to_east_deg, altitude_deg);
    ArrayTimingCorrectionResult result;
    result.geometric_delay_ns = planeWavefrontGeometricDelayNs(
        telescope_position_nwu_m, direction, wavefront_speed_m_per_ns);
    result.corrected_trigger_time_ns =
        raw_trigger_time_ns + result.geometric_delay_ns;
    return result;
}

void applyPlaneWavefrontTimingCorrection(
    std::vector<TelescopeTriggerTime>& telescope_triggers,
    const std::map<int, Vec3>& telescope_positions_nwu_m,
    double azimuth_north_to_east_deg,
    double altitude_deg,
    double wavefront_speed_m_per_ns)
{
    for (auto& trigger : telescope_triggers) {
        const auto position = telescope_positions_nwu_m.find(
            trigger.telescope_id);
        if (position == telescope_positions_nwu_m.end()) {
            throw std::runtime_error(
                "array plane-wave timing has no position for telescope " +
                std::to_string(trigger.telescope_id));
        }
        const auto correction = correctArrayTriggerTimePlaneWave(
            trigger.trigger_time_ns,
            position->second,
            azimuth_north_to_east_deg,
            altitude_deg,
            wavefront_speed_m_per_ns);
        trigger.geometric_delay_ns = correction.geometric_delay_ns;
        trigger.coincidence_time_ns = correction.corrected_trigger_time_ns;
    }
}

} // namespace lact
