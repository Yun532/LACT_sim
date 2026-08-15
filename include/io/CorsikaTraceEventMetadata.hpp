#pragma once

// CORSIKA event identity and truth resolved from EventIO metadata. Shared by
// the tracing app and the HDF5 writer so both derive event ids the same way.

#include <string>

#include "app/OpticalSimCommon.hpp"
#include "io/EventIOPhotonSource.hpp"

namespace lact {

int showerEventFromOutputEvent(int event_id, const std::string& event_id_mode);
int arrayIdFromOutputEvent(int event_id, const std::string& event_id_mode);

// Front-most mirror z in the telescope-local frame, used as the fallback
// EventIO input plane.
double mirrorFrontReferenceZ(const MirrorLayout& mirrors);

// A 2D EventIO record is an anchor on the photon line, so the mirror may sit
// at either sign of the ray parameter.
bool shouldBackprojectEventIO2d(const SourceRuntimeConfig& source_runtime_cfg);

struct OutputEventMetadata {
    int event_id = 0;
    int shower_event = 0;
    int array_id = 0;
    double energy_gev = 0.0;
    double core_x_north_m = 0.0;
    double core_y_west_m = 0.0;
    double azimuth_north_to_east_deg = 0.0;
    double array_time_offset_ns = 0.0;
    double area_weight_m2 = 0.0;
    bool has_explicit_area_weight = false;
    bool found = false;
    bool used_array_offset = false;
};

OutputEventMetadata outputEventMetadata(int event_id,
                                        const std::string& event_id_mode,
                                        const EventIOMetadata& metadata);

} // namespace lact
