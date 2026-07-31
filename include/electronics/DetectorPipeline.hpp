#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace lact::electronics {

enum class HitOrigin {
    Cherenkov,
    Nsb,
    Dark,
};

const char* hitOriginName(HitOrigin origin);
HitOrigin parseHitOrigin(const std::string& text);

// One detected primary p.e. after the optical chain and PDE.  In the
// stochastic production path primary_pe is exactly one.  Aggregated values
// are accepted by the standalone interface, but explicit microcell tracking
// expands only integral values because a fractional p.e. has no cell identity.
struct PrimaryPeHit {
    int event_id = 0;
    int telescope_id = 0;
    int pixel_id = -1;
    double time_ns = 0.0;
    double sensor_x_m = 0.0;
    double sensor_y_m = 0.0;
    double wavelength_nm = 0.0;
    double primary_pe = 1.0;
    HitOrigin origin = HitOrigin::Cherenkov;
};

struct MicrocellConfig {
    bool enabled = false;
    std::string model = "explicit_no_recovery";
    std::string layout = "uniform_interleaved";
    double sensor_size_x_m = 0.0130;
    double sensor_size_y_m = 0.0130;
    int grid_columns = 528;
    int grid_rows = 512;
    int channels_per_pixel = 8;
    int microcells_per_channel = 33792;
};

struct SinglePeConfig {
    bool enabled = false;
    std::string model = "analytic";
    std::string csv_path;
    // pe_charge: template area is normalized to one p.e.; samples contain
    // integrated p.e. charge.  mV: template amplitudes are preserved and
    // samples contain the average voltage over each sample.
    std::string unit = "pe_charge";
    double rise_ns = 1.0;
    double fall_ns = 6.0;
    double support_ns = 60.0;
    double amplitude_scale = 1.0;
};

struct SamplingConfig {
    double width_ns = 4.0;
    double start_ns = -20.0;
    double end_ns = 100.0;
};

struct CameraTriggerConfig {
    bool enabled = false;
    std::string mode = "pe_count";
    double pixel_threshold_pe = 5.0;
    double pixel_threshold_mv = 5.0;
    int multiplicity = 3;
    double coincidence_window_ns = 20.0;
};

struct DetectorPipelineConfig {
    bool enabled = false;
    MicrocellConfig microcell;
    SinglePeConfig single_pe;
    SamplingConfig sampling;
    CameraTriggerConfig camera_trigger;
    // Serialization controls. The pipeline still keeps the internal hit
    // sequence long enough to evaluate downstream stages and triggering.
    bool save_primary_sequence = false;
    bool save_fired_sequence = false;
    bool save_microcell_decisions = false;
    bool save_channel_waveforms = false;
};

struct MicrocellDecision {
    std::size_t primary_index = 0;
    int event_id = 0;
    int telescope_id = 0;
    int pixel_id = -1;
    double time_ns = 0.0;
    double sensor_x_m = 0.0;
    double sensor_y_m = 0.0;
    int grid_column = -1;
    int grid_row = -1;
    int channel_id = -1;
    int microcell_id = -1;
    bool fired = false;
    HitOrigin origin = HitOrigin::Cherenkov;
};

struct FiredCellHit {
    int event_id = 0;
    int telescope_id = 0;
    int pixel_id = -1;
    double time_ns = 0.0;
    int channel_id = -1;
    int microcell_id = -1;
    double fired_pe = 1.0;
    HitOrigin origin = HitOrigin::Cherenkov;
};

struct PixelSummary {
    double primary_cherenkov_pe = 0.0;
    double primary_nsb_pe = 0.0;
    double primary_dark_pe = 0.0;
    double fired_cherenkov_pe = 0.0;
    double fired_nsb_pe = 0.0;
    double fired_dark_pe = 0.0;
    double saturation_lost_pe = 0.0;
};

struct CameraTriggerDecision {
    bool triggered = false;
    int max_pixels_above_threshold = 0;
    int first_trigger_bin = -1;
    double trigger_time_ns = 0.0;
    std::vector<int> pixels_above_threshold;
};

struct DetectorPipelineResult {
    int event_id = 0;
    int telescope_id = 0;
    std::size_t n_pixels = 0;
    std::size_t n_samples = 0;
    std::string sample_unit = "none";
    std::vector<double> time_edges_ns;
    std::vector<double> time_centers_ns;
    std::vector<PrimaryPeHit> primary_hits;
    std::vector<MicrocellDecision> microcell_decisions;
    std::vector<FiredCellHit> fired_hits;
    std::vector<PixelSummary> pixels;
    // Row-major [pixel][sample].
    std::vector<double> waveform;
    // Optional row-major [pixel][channel][sample].
    std::vector<double> channel_waveform;
    CameraTriggerDecision camera_trigger;
};

void validateDetectorPipelineConfig(const DetectorPipelineConfig& config);

DetectorPipelineResult runDetectorPipeline(
    const DetectorPipelineConfig& config,
    std::size_t n_pixels,
    const std::vector<PrimaryPeHit>& primary_hits);

std::vector<PrimaryPeHit> generateUniformNsbPrimaryHits(
    int event_id,
    int telescope_id,
    std::size_t n_pixels,
    double rate_pe_per_ns_per_pixel,
    double start_ns,
    double end_ns,
    double sensor_size_x_m,
    double sensor_size_y_m,
    std::uint64_t seed);

} // namespace lact::electronics
