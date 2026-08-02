#include "electronics/DetectorPipeline.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

double waveformSum(const lact::electronics::DetectorPipelineResult& result,
                   std::size_t pixel)
{
    double sum = 0.0;
    for (std::size_t bin = 0; bin < result.n_samples; ++bin) {
        sum += result.waveform[pixel * result.n_samples + bin];
    }
    return sum;
}

} // namespace

int main()
{
    using namespace lact::electronics;
    DetectorPipelineConfig config;
    config.microcell.enabled = true;
    config.save_microcell_decisions = true;
    config.single_pe.enabled = false;
    config.sampling.start_ns = 0.0;
    config.sampling.end_ns = 40.0;
    config.sampling.width_ns = 4.0;
    config.camera_trigger.enabled = true;
    config.camera_trigger.mode = "pe_count";
    config.camera_trigger.pixel_threshold_pe = 1.0;
    config.camera_trigger.multiplicity = 2;
    config.camera_trigger.coincidence_window_ns = 8.0;

    std::vector<PrimaryPeHit> hits = {
        {1, 2, 0, 4.0, 0.0, 0.0, 420.0, 1.0,
         HitOrigin::Cherenkov},
        {1, 2, 0, 5.0, 0.0, 0.0, 420.0, 1.0,
         HitOrigin::Cherenkov},
        {1, 2, 0, 6.0, 0.001, 0.0, 420.0, 1.0,
         HitOrigin::Nsb},
        {1, 2, 1, 6.0, 0.0, 0.0, 420.0, 1.0,
         HitOrigin::Cherenkov},
    };
    auto result = runDetectorPipeline(config, 2, hits);
    require(result.primary_hits.size() == 4,
            "all primary hits must be preserved");
    require(result.microcell_decisions.size() == 4,
            "one decision is required per stochastic primary p.e.");
    require(result.fired_hits.size() == 3,
            "one repeated cell must be suppressed");
    require(result.pixels[0].saturation_lost_pe == 1.0,
            "repeated cell must be counted as saturation loss");
    require(result.pixels[0].primary_cherenkov_pe == 2.0,
            "primary Cherenkov truth must be preserved");
    require(result.pixels[0].fired_cherenkov_pe == 1.0,
            "only one Cherenkov avalanche may fire the repeated cell");
    require(result.pixels[0].fired_nsb_pe == 1.0,
            "different-position NSB hit must fire");
    require(result.camera_trigger.triggered,
            "two pixels in the coincidence window must trigger");
    require(result.camera_trigger.first_trigger_bin == 1 &&
                std::abs(result.camera_trigger.trigger_time_ns - 6.0) < 1.0e-12,
            "camera trigger time must be the causal window-end sample");

    MicrocellConfig s17351;
    s17351.layout = "s17351_tiled_2x4";
    s17351.sensor_size_x_m = 0.0134;
    s17351.sensor_size_y_m = 0.0134;
    const auto active_cell = mapMicrocellPosition(
        s17351, -0.0066875, -0.0066875);
    require(active_cell.inside_channel,
            "S17351 cell center must map to a channel");
    require(active_cell.channel_id == 7 && active_cell.microcell_id == 0,
            "S17351 front-side B-4 lower-left mapping is incorrect");
    const auto channel_gap = mapMicrocellPosition(
        s17351, 0.0, -0.0060);
    require(channel_gap.channel_gap && !channel_gap.inside_channel,
            "S17351 0.2 mm inter-channel gap must be explicit");
    const double expected_channel_fraction =
        8.0 * 0.0066 * 0.0032 / (0.0134 * 0.0134);
    require(std::abs(interChannelActiveFraction(s17351) -
                     expected_channel_fraction) < 1.0e-12,
            "S17351 channel active fraction must follow package geometry");
    const auto cell_edge = mapMicrocellPosition(
        s17351, -0.006699, -0.0066875);
    require(cell_edge.inside_channel &&
                cell_edge.microcell_id == active_cell.microcell_id,
            "all positions inside a 25 um pitch must map to that microcell");

    DetectorPipelineConfig gap_config;
    gap_config.microcell = s17351;
    gap_config.microcell.enabled = true;
    gap_config.microcell.saturation_enabled = false;
    gap_config.save_microcell_decisions = true;
    gap_config.sampling.start_ns = 0.0;
    gap_config.sampling.end_ns = 40.0;
    gap_config.sampling.width_ns = 4.0;
    std::vector<PrimaryPeHit> gap_hits = {
        {1, 2, 0, 4.0, -0.0066875, -0.0066875, 420.0, 1.0,
         HitOrigin::Cherenkov},
        {1, 2, 0, 5.0, -0.006699, -0.0066875, 420.0, 1.0,
         HitOrigin::Cherenkov},
    };
    const auto gap_result = runDetectorPipeline(gap_config, 1, gap_hits);
    require(gap_result.fired_hits.size() == 2,
            "datasheet-PDE p.e. must not be rejected inside a microcell");
    require(gap_result.pixels[0].gap_lost_pe == 0.0,
            "no sub-microcell fill-factor loss may be applied");
    require(gap_result.pixels[0].saturation_lost_pe == 0.0,
            "gap-only mode must not report saturation");

    config.single_pe.enabled = true;
    config.single_pe.model = "analytic";
    config.single_pe.unit = "pe_charge";
    config.single_pe.support_ns = 30.0;
    config.camera_trigger.enabled = false;
    result = runDetectorPipeline(config, 2, hits);
    require(std::abs(waveformSum(result, 0) - 2.0) < 1.0e-6,
            "normalized waveform must integrate to fired p.e.");
    require(std::abs(waveformSum(result, 1) - 1.0) < 1.0e-6,
            "single fired cell must integrate to one p.e.");

    config.microcell.enabled = false;
    config.single_pe.unit = "mv";
    config.single_pe.amplitude_scale = 10.0;
    config.camera_trigger.enabled = true;
    config.camera_trigger.mode = "voltage";
    config.camera_trigger.pixel_threshold_mv = 1.0;
    config.camera_trigger.multiplicity = 2;
    result = runDetectorPipeline(config, 2, hits);
    require(result.sample_unit == "mV",
            "absolute waveform mode must expose mV units");
    require(result.camera_trigger.triggered,
            "voltage threshold must operate on sampled waveform");

    std::cout << "detector pipeline tests passed\n";
    return 0;
}
