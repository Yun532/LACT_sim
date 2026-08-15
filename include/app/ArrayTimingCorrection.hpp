#pragma once

#include <map>
#include <vector>

#include "app/TriggerResponse.hpp"
#include "core/Vec3.hpp"

namespace lact {

struct ArrayTimingCorrectionResult {
    double geometric_delay_ns = 0.0;
    double corrected_trigger_time_ns = 0.0;
};

// Same observation-level refractivity approximation used by sim_telarray.
double airRefractiveIndexAtAltitudeKm(double altitude_km);
double airLightSpeedAtAltitudeM(double altitude_m);

// Unit vector pointing from the array toward the shower source in the
// CORSIKA magnetic North-West-Up frame. Azimuth increases North-to-East,
// while the stored array y coordinate is positive toward West.
Vec3 corsikaNwuViewingDirection(double azimuth_north_to_east_deg,
                                double altitude_deg);

double showerAxisImpactParameterM(
    const Vec3& telescope_position_nwu_m,
    const Vec3& shower_core_nwu_m,
    const Vec3& viewing_direction_nwu);

// Plane-wave correction used before applying an array coincidence window.
// A positive projection means that the wavefront reaches this telescope
// earlier, so the projection/c term is added to the raw local trigger time.
double planeWavefrontGeometricDelayNs(
    const Vec3& telescope_position_nwu_m,
    const Vec3& viewing_direction_nwu,
    double wavefront_speed_m_per_ns);

ArrayTimingCorrectionResult correctArrayTriggerTimePlaneWave(
    double raw_trigger_time_ns,
    const Vec3& telescope_position_nwu_m,
    double azimuth_north_to_east_deg,
    double altitude_deg,
    double wavefront_speed_m_per_ns);

void applyPlaneWavefrontTimingCorrection(
    std::vector<TelescopeTriggerTime>& telescope_triggers,
    const std::map<int, Vec3>& telescope_positions_nwu_m,
    double azimuth_north_to_east_deg,
    double altitude_deg,
    double wavefront_speed_m_per_ns);

} // namespace lact
