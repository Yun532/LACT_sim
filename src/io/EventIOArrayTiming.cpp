#include "io/EventIOArrayTiming.hpp"

#include "app/ArrayTimingCorrection.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

namespace lact {
namespace {

int showerEventForOutputEvent(int output_event_id,
                              const std::string& event_id_mode,
                              const EventIOMetadata& metadata)
{
    const auto identity = metadata.output_event_identity.find(output_event_id);
    if (identity != metadata.output_event_identity.end()) {
        return identity->second.first;
    }
    const std::string mode = lowerCopy(trim(event_id_mode));
    if (mode == "event_array100" || mode == "runid") {
        return output_event_id / 100;
    }
    return output_event_id;
}

std::map<int, Vec3> telescopePositions(const EventIOMetadata& metadata,
                                      const TelescopeConfig& telescope_cfg)
{
    std::map<int, Vec3> positions;
    for (const auto& telescope : metadata.telescopes) {
        positions[telescope.telescope_id] = {
            telescope.x_m, telescope.y_m, telescope.z_m};
    }
    if (positions.empty()) {
        positions[telescope_cfg.id] = telescope_cfg.position_m;
    }
    return positions;
}

} // namespace

double resolveEventIOArrayWavefrontSpeedMPerNs(
    const TriggerConfig& trigger_cfg,
    const EventIOMetadata& metadata)
{
    if (trigger_cfg.array_wavefront_speed_m_per_ns > 0.0) {
        return trigger_cfg.array_wavefront_speed_m_per_ns;
    }
    if (!std::isfinite(metadata.observation_altitude_m)) {
        throw std::runtime_error(
            "automatic array wavefront speed requires EventIO observation altitude");
    }
    return airLightSpeedAtAltitudeM(metadata.observation_altitude_m);
}

void applyEventIOArrayTimingCorrection(
    std::vector<TelescopeTriggerTime>& telescope_triggers,
    int output_event_id,
    const std::string& event_id_mode,
    const TriggerConfig& trigger_cfg,
    const TelescopeConfig& telescope_cfg,
    const EventIOMetadata& metadata)
{
    if (trigger_cfg.array_time_correction == "none" ||
        telescope_triggers.empty()) {
        return;
    }
    if (trigger_cfg.array_time_correction != "plane_wave") {
        throw std::runtime_error("unsupported array timing correction mode: " +
                                 trigger_cfg.array_time_correction);
    }

    const int shower_event = showerEventForOutputEvent(
        output_event_id, event_id_mode, metadata);
    const auto event = std::find_if(
        metadata.events.begin(), metadata.events.end(),
        [shower_event](const EventIOEventHeader& candidate) {
            return candidate.shower_event_id == shower_event;
        });
    if (event == metadata.events.end()) {
        throw std::runtime_error(
            "array plane-wave timing has no EventIO shower metadata for event " +
            std::to_string(output_event_id));
    }
    const double altitude_deg = std::isfinite(event->altitude_deg)
        ? event->altitude_deg
        : 90.0 - event->theta_deg;
    if (!std::isfinite(altitude_deg) ||
        !std::isfinite(event->azimuth_north_to_east_deg)) {
        throw std::runtime_error(
            "array plane-wave timing requires EventIO shower altitude and azimuth");
    }

    applyPlaneWavefrontTimingCorrection(
        telescope_triggers,
        telescopePositions(metadata, telescope_cfg),
        event->azimuth_north_to_east_deg,
        altitude_deg,
        resolveEventIOArrayWavefrontSpeedMPerNs(trigger_cfg, metadata));
}

} // namespace lact
