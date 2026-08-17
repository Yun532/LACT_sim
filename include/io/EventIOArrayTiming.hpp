#pragma once

#include <map>
#include <string>
#include <vector>

#include "app/OpticalSimCommon.hpp"
#include "app/TriggerResponse.hpp"
#include "io/EventIOPhotonSource.hpp"

namespace lact {

double resolveEventIOArrayWavefrontSpeedMPerNs(
    const TriggerConfig& trigger_cfg,
    const EventIOMetadata& metadata);

// Return the configured plane-wave delay for every requested telescope.
// This quantity depends only on event geometry, never on whether the camera
// triggered.  In `none` mode every requested telescope receives 0 ns.
std::map<int, double> eventIOArrayGeometricDelaysNs(
    const std::vector<int>& telescope_ids,
    int output_event_id,
    const std::string& event_id_mode,
    const TriggerConfig& trigger_cfg,
    const TelescopeConfig& telescope_cfg,
    const EventIOMetadata& metadata);

// Apply the configured EventIO array timing model in-place. The event identity,
// shower direction, and NWU telescope positions are resolved once here so ROOT
// and HDF5 cannot drift into separate timing conventions.
void applyEventIOArrayTimingCorrection(
    std::vector<TelescopeTriggerTime>& telescope_triggers,
    int output_event_id,
    const std::string& event_id_mode,
    const TriggerConfig& trigger_cfg,
    const TelescopeConfig& telescope_cfg,
    const EventIOMetadata& metadata);

} // namespace lact
