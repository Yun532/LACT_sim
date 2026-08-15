#include "io/CorsikaTraceEventMetadata.hpp"

#include <algorithm>

#include "app/OpticalSimCommon.hpp"

namespace lact {

int showerEventFromOutputEvent(int event_id, const std::string& event_id_mode)
{
    const std::string mode = lowerCopy(trim(event_id_mode));
    if (mode == "event_array100" || mode == "runid") {
        return event_id / 100;
    }
    return event_id;
}

int arrayIdFromOutputEvent(int event_id, const std::string& event_id_mode)
{
    const std::string mode = lowerCopy(trim(event_id_mode));
    if (mode == "event_array100" || mode == "runid") {
        return event_id % 100;
    }
    return 0;
}

double mirrorFrontReferenceZ(const MirrorLayout& mirrors)
{
    double z = -std::numeric_limits<double>::infinity();
    for (const auto& tile : mirrors.tiles()) {
        z = std::max(z, tile.center.z);
    }
    return std::isfinite(z) ? z : -16.0;
}

bool shouldBackprojectEventIO2d(const SourceRuntimeConfig& source_runtime_cfg)
{
    const std::string mode = lowerCopy(trim(source_runtime_cfg.eventio_2d_plane_mode));
    if (mode == "forward") {
        return false;
    }
    // A 2D EventIO bunch supplies an anchor on the unperturbed photon line,
    // not a creation point. In both auto and explicit backproject modes the
    // mirror may therefore lie at either sign of t relative to that anchor.
    return true;
}


OutputEventMetadata outputEventMetadata(int event_id,
                                        const std::string& event_id_mode,
                                        const EventIOMetadata& metadata)
{
    OutputEventMetadata out;
    out.event_id = event_id;
    const auto identity = metadata.output_event_identity.find(event_id);
    if (identity != metadata.output_event_identity.end()) {
        out.shower_event = identity->second.first;
        out.array_id = identity->second.second;
    } else {
        out.shower_event = showerEventFromOutputEvent(event_id, event_id_mode);
        out.array_id = arrayIdFromOutputEvent(event_id, event_id_mode);
    }

    auto event_it = std::find_if(
        metadata.events.begin(), metadata.events.end(),
        [&out](const EventIOEventHeader& event) {
            return event.shower_event_id == out.shower_event;
        });
    if (event_it == metadata.events.end()) {
        return out;
    }

    out.found = true;
    out.energy_gev = event_it->energy_gev;
    out.core_x_north_m = event_it->core_x_m;
    out.core_y_west_m = event_it->core_y_m;
    out.azimuth_north_to_east_deg = event_it->azimuth_north_to_east_deg;

    if (auto offsets = metadata.arrayOffsetsForShower(out.shower_event)) {
        out.array_time_offset_ns = offsets->time_offset_ns;
        const std::size_t offset_index = static_cast<std::size_t>(out.array_id);
        if (out.array_id >= 0 && offset_index < offsets->x_m.size() &&
            offset_index < offsets->y_m.size()) {
            // MC_TELOFF stores the offset of the detector array with respect
            // to the shower core.  The core position in the telescope/input
            // array frame is therefore the opposite vector.
            out.core_x_north_m = -offsets->x_m[offset_index];
            out.core_y_west_m = -offsets->y_m[offset_index];
            out.used_array_offset = true;
        }
        if (out.array_id >= 0 && offset_index < offsets->weight.size()) {
            out.area_weight_m2 = offsets->weight[offset_index];
            out.has_explicit_area_weight = offsets->has_explicit_weights;
        }
    }
    return out;
}

} // namespace lact
