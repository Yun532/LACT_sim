#pragma once

// Per-run configuration flags, accumulators and row types shared by the
// tracing app and its output/report modules.

#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "core/PhotonBunch.hpp"
#include "io/CorsikaTraceOutputTypes.hpp"
#include "optics/OpticalSurfaceHit.hpp"

namespace lact {

struct CollectorDebugConfig {
    bool photon_output = false;
    std::string photon_csv = "collector_debug_photons.csv";
    std::uint64_t max_photons = 100000;
};

struct WavelengthRange {
    double min_nm = 0.0;
    double max_nm = 0.0;
};

struct ProfileConfig {
    bool enabled = false;
};

struct AtmosphereHistogramConfig {
    bool enabled = false;
    std::string csv_path = "atmosphere_height_histogram.csv";
    double min_altitude_km = 4.4;
    double max_altitude_km = 100.0;
    double bin_width_km = 1.0;
};

struct AtmosphereHistogramBin {
    double low_km = 0.0;
    double high_km = 0.0;
    std::uint64_t bunches = 0;
    double before_weight = 0.0;
    double after_weight = 0.0;
    double theory_weight = 0.0;
};

struct ProfileStats {
    double eventio_stream_s = 0.0;
    double transform_s = 0.0;
    double trace_to_plane_s = 0.0;
    double obstruction_s = 0.0;
    double camera_response_s = 0.0;
    double whiteboard_accumulate_s = 0.0;
    double camera_accumulate_s = 0.0;
    double hdf5_write_s = 0.0;
};

struct TelescopeEventAccumulator {
    int telescope_id = 0;
    std::set<int> output_events;
    std::uint64_t input_bunches = 0;
    double input_photons = 0.0;
    std::uint64_t blocked_by_obstruction = 0;
    std::uint64_t blocked_incoming = 0;
    std::uint64_t blocked_reflected = 0;
    std::uint64_t hit_mirror_before_obstruction = 0;
    std::uint64_t hit_output_before_obstruction = 0;
    std::uint64_t hit_mirror = 0;
    std::uint64_t hit_output_plane = 0;
    std::uint64_t hit_camera = 0;
    std::uint64_t accepted_camera = 0;
    std::uint64_t lost_between_pixels = 0;
    std::set<int> unique_pixels;
    double weighted_signal = 0.0;
    double weighted_time_sum = 0.0;
    double weighted_time2_sum = 0.0;
    double first_cherenkov_time_ns = std::numeric_limits<double>::infinity();
};

struct CollectorDebugPhotonRow {
    int event_id = 0;
    int telescope_id = 0;
    int pixel_id = -1;
    int accepted = 0;
    int collector_reflections = 0;
    int collector_reflection_limit_reached = 0;
    double time_ns = 0.0;
    double collector_path_length_m = 0.0;
    double collector_time_delay_ns = 0.0;
    double wavelength_nm = 0.0;
    double photon_weight = 0.0;
    double relative_efficiency = 0.0;
    double signal_weight = 0.0;
    double exit_x_m = 0.0;
    double exit_y_m = 0.0;
    double exit_z_m = 0.0;
    double dir_u = 0.0;
    double dir_v = 0.0;
    double dir_w = 0.0;
    int sipm_gap_rejected = 0;
    int sipm_channel_gap_rejected = 0;
    int sipm_grid_column = -1;
    int sipm_grid_row = -1;
    int sipm_channel_id = -1;
    int sipm_microcell_id = -1;
};

} // namespace lact
