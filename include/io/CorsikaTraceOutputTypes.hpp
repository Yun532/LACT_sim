#pragma once

#include <cmath>
#include <cstdint>
#include <limits>
#include <set>
#include <string>
#include <tuple>

#include "electronics/DetectorPipeline.hpp"

namespace lact {

struct CorsikaTraceOutputConfig {
    std::string hits_csv = "corsika_whiteboard_hits.csv";
    std::string pixel_csv = "corsika_pixel_image.csv";
    std::string summary_csv = "corsika_trace_summary.csv";
    std::string waveform_csv = "corsika_waveforms.csv";
    std::string trigger_csv = "corsika_triggers.csv";
    std::string primary_pe_csv = "corsika_primary_pe.csv";
    std::string fired_pe_csv = "corsika_fired_pe.csv";
    std::string mirror_diagnostic_csv;
    std::string hdf5_path = "corsika_trace.h5";
    std::string lact_root_path = "lact_events.root";
    std::string lact_profile = "image_pe";
    std::string format = "hdf5";
    std::string hdf5_storage = "dense";
    std::string hdf5_waveform_storage = "sparse";
    bool hdf5_write_components = false;
    bool hdf5_write_waveforms = true;
    bool lact_root_write_components = false;
    bool lact_root_enabled = false;
    bool save_only_triggered = false;
    bool write_pixel_time_stats = false;
    bool whiteboard_emitter_info = false;
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
    double cherenkov_pe = 0.0;
    double nsb_pe = 0.0;
    double dark_pe = 0.0;
};

struct RawWaveformHit {
    int event_id = 0;
    int telescope_id = 0;
    int pixel_id = -1;
    double time_ns = 0.0;
    std::uint64_t photon_count = 0;
    double pe = 0.0;
    double sensor_x_m = 0.0;
    double sensor_y_m = 0.0;
    double wavelength_nm = 0.0;
    electronics::HitOrigin origin = electronics::HitOrigin::Cherenkov;
};

struct WhiteboardHdf5Row {
    std::int64_t event_id;
    std::int32_t telescope_id;
    std::int64_t photon_index;
    std::int32_t mirror_id;
    float surface_x_m;
    float surface_y_m;
    float surface_z_m;
    float u_m;
    float v_m;
    float dir_x;
    float dir_y;
    float dir_z;
    float time_ns;
    float wavelength_nm;
    float weight;
    float relative_efficiency;
    float signal_weight;
    std::uint8_t has_emitter;
    float emitter_mass_gev;
    float emitter_charge;
    float emitter_energy_gev;
    float emitter_time_ns;
};

inline bool outputWantsCsv(const CorsikaTraceOutputConfig& cfg)
{
    return cfg.format == "csv" || cfg.format == "both";
}

inline bool outputWantsHdf5(const CorsikaTraceOutputConfig& cfg)
{
    return cfg.format == "hdf5" || cfg.format == "h5" || cfg.format == "both";
}

inline bool outputWantsLactRoot(const CorsikaTraceOutputConfig& cfg)
{
    return cfg.lact_root_enabled;
}

// Waveform time-bin geometry. Inline so both the app and the HDF5
// writer resolve the same definition.
inline std::size_t waveformBinCount(const WaveformOutputConfig& cfg)
{
    if (!cfg.enabled) {
        return 0;
    }
    const double span = cfg.time_window_end_ns - cfg.time_window_start_ns;
    return static_cast<std::size_t>(std::llround(span / cfg.time_bin_width_ns));
}

inline int waveformBinForTime(const WaveformOutputConfig& cfg, double time_ns)
{
    if (!cfg.enabled ||
        time_ns < cfg.time_window_start_ns ||
        time_ns >= cfg.time_window_end_ns) {
        return -1;
    }
    const auto bin = static_cast<int>(
        std::floor((time_ns - cfg.time_window_start_ns) / cfg.time_bin_width_ns));
    const auto n_bins = static_cast<int>(waveformBinCount(cfg));
    return bin >= 0 && bin < n_bins ? bin : -1;
}

inline bool waveformUsesImageMeanReference(const WaveformOutputConfig& cfg)
{
    return cfg.enabled && cfg.time_reference == "image_mean";
}

inline bool waveformUsesImageReference(const WaveformOutputConfig& cfg)
{
    return cfg.enabled &&
        (cfg.time_reference == "image_mean" || cfg.time_reference == "image_first");
}

} // namespace lact
