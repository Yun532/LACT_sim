#include "app/ArrayTimingCorrection.hpp"
#include "io/EventIOArrayTiming.hpp"

#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>

using namespace lact;

namespace {

bool check(bool condition, const char* message)
{
    if (!condition) {
        std::cerr << message << '\n';
    }
    return condition;
}

bool near(double a, double b, double tolerance = 1.0e-12)
{
    return std::abs(a - b) <= tolerance;
}

} // namespace

int main()
{
    bool ok = true;
    constexpr double c = 0.299792458;
    ok &= check(near(airRefractiveIndexAtAltitudeKm(4.4),
                     1.00018065764828, 1.0e-14),
                "observation-level refractive index differs from sim_telarray");
    ok &= check(near(airLightSpeedAtAltitudeM(4400.0),
                     0.299738307982179, 1.0e-14),
                "observation-level air light speed differs from sim_telarray");

    const Vec3 north = corsikaNwuViewingDirection(0.0, 0.0);
    ok &= check(near(north.x, 1.0) && near(north.y, 0.0) && near(north.z, 0.0),
                "azimuth 0 must point North (+x)");
    const Vec3 east = corsikaNwuViewingDirection(90.0, 0.0);
    ok &= check(near(east.x, 0.0, 1.0e-10) && near(east.y, -1.0) &&
                    near(east.z, 0.0),
                "azimuth 90 must point East (-West array y)");
    const Vec3 zenith = corsikaNwuViewingDirection(37.0, 90.0);
    ok &= check(near(zenith.x, 0.0, 1.0e-10) &&
                    near(zenith.y, 0.0, 1.0e-10) && near(zenith.z, 1.0),
                "altitude 90 must point Up (+z)");

    const double delay = planeWavefrontGeometricDelayNs(
        Vec3{30.0, 0.0, 0.0}, north, c);
    ok &= check(near(delay, 30.0 / c),
                "northward plane-wave delay has wrong magnitude or sign");

    const double reference_raw_time = 1000.0;
    const auto displaced = correctArrayTriggerTimePlaneWave(
        reference_raw_time - delay, Vec3{30.0, 0.0, 0.0},
        0.0, 0.0, c);
    ok &= check(near(displaced.corrected_trigger_time_ns, reference_raw_time),
                "plane-wave correction did not align telescope times");

    std::vector<TelescopeTriggerTime> triggers{
        {1, reference_raw_time},
        {2, reference_raw_time - delay},
    };
    applyPlaneWavefrontTimingCorrection(
        triggers, {{1, Vec3{0.0, 0.0, 0.0}},
                   {2, Vec3{30.0, 0.0, 0.0}}},
        0.0, 0.0, c);
    ok &= check(near(triggers[0].coincidence_time_ns,
                     triggers[1].coincidence_time_ns),
                "vector plane-wave correction did not align telescope times");

    EventIOMetadata metadata;
    EventIOEventHeader event;
    event.shower_event_id = 12;
    event.altitude_deg = 0.0;
    event.azimuth_north_to_east_deg = 0.0;
    metadata.events.push_back(event);
    metadata.observation_altitude_m = 4400.0;
    metadata.output_event_identity[1203] = {12, 3};
    metadata.telescopes.push_back({1, 0.0, 0.0, 0.0, 0.0});
    metadata.telescopes.push_back({2, 30.0, 0.0, 0.0, 0.0});
    TriggerConfig eventio_trigger;
    eventio_trigger.array_time_correction = "plane_wave";
    eventio_trigger.array_wavefront_speed_m_per_ns = c;
    TelescopeConfig telescope;
    std::vector<TelescopeTriggerTime> eventio_triggers{
        {1, reference_raw_time},
        {2, reference_raw_time - delay},
    };
    applyEventIOArrayTimingCorrection(
        eventio_triggers, 1203, "event_array100", eventio_trigger,
        telescope, metadata);
    ok &= check(near(eventio_triggers[0].coincidence_time_ns,
                     eventio_triggers[1].coincidence_time_ns),
                "EventIO timing context did not align telescope times");
    const auto all_delays = eventIOArrayGeometricDelaysNs(
        {1, 2}, 1203, "event_array100", eventio_trigger, telescope, metadata);
    ok &= check(near(all_delays.at(1), 0.0) &&
                    near(all_delays.at(2), delay),
                "EventIO geometric delays were not available for all telescopes");
    eventio_trigger.array_time_correction = "none";
    const auto disabled_delays = eventIOArrayGeometricDelaysNs(
        {1, 2}, 1203, "event_array100", eventio_trigger, telescope, metadata);
    ok &= check(near(disabled_delays.at(1), 0.0) &&
                    near(disabled_delays.at(2), 0.0),
                "disabled timing correction must serialize finite zero delays");
    eventio_trigger.array_time_correction = "plane_wave";
    eventio_trigger.array_wavefront_speed_m_per_ns = 0.0;
    ok &= check(near(resolveEventIOArrayWavefrontSpeedMPerNs(
                         eventio_trigger, metadata),
                     airLightSpeedAtAltitudeM(4400.0)),
                "EventIO automatic wavefront speed did not use observation altitude");

    bool rejected = false;
    try {
        (void)corsikaNwuViewingDirection(0.0, 100.0);
    } catch (const std::runtime_error&) {
        rejected = true;
    }
    ok &= check(rejected, "invalid altitude was not rejected");
    return ok ? 0 : 1;
}
