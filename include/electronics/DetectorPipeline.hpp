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
    // False for NSB warm-up/cool-down hits generated outside the stored
    // sampling window. Such hits participate in microcell occupancy and
    // waveform convolution, but not in the window-integrated image truth.
    bool count_in_integrated_image = true;
};

struct MicrocellConfig {
    bool enabled = false;
    bool saturation_enabled = true;
    // If enabled, a previously fired cell recharges exponentially.  The
    // emitted charge-equivalent p.e. fraction is
    // 1 - exp(-delta_t / recovery_time_ns).  The default time constant is a
    // provisional value and must be replaced when device calibration exists.
    bool recovery_enabled = false;
    double recovery_time_ns = 10.0;
    std::string model = "explicit_no_recovery";
    std::string layout = "uniform_interleaved";
    double sensor_size_x_m = 0.0130;
    double sensor_size_y_m = 0.0130;
    int grid_columns = 528;
    int grid_rows = 512;
    int channels_per_pixel = 8;
    int microcells_per_channel = 33792;
    // S17351 physical layout.  The 8 channels are 2 x 4 tiles, each with
    // a 6.6 x 3.2 mm^2 active envelope and a 0.2 mm inter-channel gap.
    int channel_columns = 2;
    int channel_rows = 4;
    double channel_size_x_m = 0.0066;
    double channel_size_y_m = 0.0032;
    double channel_gap_x_m = 0.0002;
    double channel_gap_y_m = 0.0002;
    int microcell_columns_per_channel = 264;
    int microcell_rows_per_channel = 128;
    // The tiled layout contains the specified 0.2 mm inter-channel gaps.
    // Datasheet PDE is already effective, so no sub-microcell active region
    // or fill-factor rejection is modelled.
    // If true, the supplied PDE is also assumed to have been averaged over
    // the full sensor envelope including the inter-channel gaps.
    bool pde_includes_inter_channel_gaps = true;
};

struct MicrocellAddress {
    bool inside_sensor = false;
    bool inside_channel = false;
    bool channel_gap = false;
    int grid_column = -1;
    int grid_row = -1;
    int channel_id = -1;
    int microcell_id = -1;
};

struct SinglePeConfig {
    struct ChargeFluctuationConfig {
        bool enabled = false;
        std::string model = "empirical";
        std::string csv_path;
        std::string csv_column = "charge_factor_mean_one";
        std::vector<double> empirical_samples;
    };

    struct TimeJitterConfig {
        bool enabled = false;
        std::string model = "gaussian";
        double sigma_ns = 0.0;
    };

    bool enabled = false;
    std::string model = "analytic";
    std::string csv_path;
    // measured_csv only: preserve the CSV time coordinate when peak is used;
    // shift the first CSV point to zero when onset is used.
    std::string template_time_reference = "onset";
    // pe_charge: template area is normalized to one p.e.; samples contain
    // integrated p.e. charge.  mV: template amplitudes are preserved and
    // samples contain the average voltage over each sample.
    std::string unit = "pe_charge";
    double rise_ns = 1.0;
    double fall_ns = 6.0;
    double support_ns = 60.0;
    double amplitude_scale = 1.0;
    ChargeFluctuationConfig charge_fluctuation;
    TimeJitterConfig time_jitter;
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
    std::uint64_t random_seed = 20260810ULL;
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
    bool gap_rejected = false;
    bool saturation_rejected = false;
    HitOrigin origin = HitOrigin::Cherenkov;
};

struct FiredCellHit {
    int event_id = 0;
    int telescope_id = 0;
    int pixel_id = -1;
    double time_ns = 0.0;
    int channel_id = -1;
    int microcell_id = -1;
    // Charge-equivalent avalanche amplitude.  It is one for a fully charged
    // cell and can be between zero and one during exponential recovery.
    double fired_pe = 1.0;
    HitOrigin origin = HitOrigin::Cherenkov;
    // These stochastic electronics quantities affect the waveform only.
    // fired_pe remains the integer avalanche truth used by saturation studies.
    double charge_factor = 1.0;
    double time_jitter_ns = 0.0;
    bool count_in_integrated_image = true;
};

struct PrimaryHitGenerationWindow {
    double start_ns = 0.0;
    double end_ns = 0.0;
};

struct PixelSummary {
    double primary_cherenkov_pe = 0.0;
    double primary_nsb_pe = 0.0;
    double primary_dark_pe = 0.0;
    double fired_cherenkov_pe = 0.0;
    double fired_nsb_pe = 0.0;
    double fired_dark_pe = 0.0;
    double gap_lost_pe = 0.0;
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
    double single_pe_area_mv_ns = 0.0;
    bool charge_fluctuation_enabled = false;
    bool time_jitter_enabled = false;
    std::string template_time_reference = "none";
    // The common high-resolution single-p.e. pulse is file-level calibration
    // metadata. Writers serialize it once, never once per event or pixel.
    std::vector<double> reference_pulse_time_ns;
    std::vector<double> reference_pulse_amplitude;
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

MicrocellAddress mapMicrocellPosition(
    const MicrocellConfig& config,
    double sensor_x_m,
    double sensor_y_m);

double interChannelActiveFraction(const MicrocellConfig& config);

// Return the primary-hit interval whose single-p.e. pulses can overlap the
// stored sampling window. With single-p.e. shaping disabled this is exactly
// the sampling window. A finite six-sigma guard is added for enabled Gaussian
// timing jitter.
PrimaryHitGenerationWindow waveformContributingPrimaryWindow(
    const DetectorPipelineConfig& config);

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
    const MicrocellConfig& microcell,
    std::uint64_t seed);

} // namespace lact::electronics
