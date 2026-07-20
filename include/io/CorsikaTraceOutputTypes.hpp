#pragma once

#include <cstdint>
#include <limits>
#include <set>
#include <string>
#include <tuple>

namespace lact {

struct CorsikaTraceOutputConfig {
    std::string hits_csv = "corsika_whiteboard_hits.csv";
    std::string pixel_csv = "corsika_pixel_image.csv";
    std::string summary_csv = "corsika_trace_summary.csv";
    std::string hdf5_path = "corsika_trace.h5";
    std::string lact_root_path = "lact_events.root";
    std::string lact_profile = "image_pe";
    std::string format = "hdf5";
    std::string hdf5_storage = "dense";
    std::string hdf5_waveform_storage = "sparse";
    bool hdf5_write_components = false;
    bool hdf5_write_waveforms = true;
    bool hdf5_write_parent_component_waveforms = false;
    bool whiteboard_emitter_info = false;
    bool lact_root_write_components = false;
    bool lact_root_enabled = false;
    bool save_only_triggered = false;
    bool write_pixel_time_stats = false;
    double lact_root_auto_flush_mb = 200.0;
    int lact_root_flush_events = 20;
};

struct WaveformOutputConfig {
    bool enabled = false;
    std::string source = "none";
    std::string time_reference = "absolute";
    double time_bin_width_ns = 1.0;
    double time_window_start_ns = 0.0;
    double time_window_end_ns = 100.0;
};

struct TraceSummary {
    int event_id = 0;
    int telescope_id = 0;
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

using SummaryKey = std::pair<int, int>;
using WaveformKey = std::tuple<int, int, int, int>;

struct WaveformPixelAccumulator {
    int event_id = 0;
    int telescope_id = 0;
    int pixel_id = -1;
    int time_bin = -1;
    std::uint64_t photon_count = 0;
    double pe = 0.0;
};

struct RawWaveformHit {
    int event_id = 0;
    int telescope_id = 0;
    int pixel_id = -1;
    double time_ns = 0.0;
    std::uint64_t photon_count = 0;
    double pe = 0.0;
};

struct RawParentComponentWaveformHit {
    int event_id = 0;
    int telescope_id = 0;
    int pixel_id = -1;
    int component_id = 5;
    double time_ns = 0.0;
    std::uint64_t photon_count = 0;
    double pe = 0.0;
};

} // namespace lact
