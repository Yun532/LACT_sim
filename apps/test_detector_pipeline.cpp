#include "electronics/DetectorPipeline.hpp"
#include "io/CameraElectronicsEvent.hpp"

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

int main(int argc, char** argv)
{
    using namespace lact::electronics;
    lact::WaveformOutputConfig trigger_time_axis;
    trigger_time_axis.enabled = false;
    trigger_time_axis.time_reference = "image_first";
    lact::TraceSummary trigger_time_summary;
    trigger_time_summary.first_cherenkov_time_ns = -321.5;
    require(std::abs(lact::waveformReferenceTimeNs(
                         trigger_time_axis, trigger_time_summary) + 321.5) <
                1.0e-12,
            "p.e.-count trigger must retain per-telescope image-first timing "
            "when waveform serialization is disabled");

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
    require(lact::hasFinalIntegratedImageSignal(result),
            "a fired p.e. must make the final integrated image non-empty");
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

    DetectorPipelineResult empty_result;
    empty_result.pixels.resize(2);
    empty_result.pixels[0].primary_cherenkov_pe = 3.0;
    empty_result.pixels[0].gap_lost_pe = 3.0;
    require(!lact::hasFinalIntegratedImageSignal(empty_result),
            "primary p.e. rejected before firing must leave an empty image");

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
    require(result.reference_pulse_time_ns.size() ==
                result.reference_pulse_amplitude.size() &&
                result.reference_pulse_time_ns.size() > 2,
            "single-p.e. reference pulse must be exposed as file metadata");
    double reference_area = 0.0;
    for (std::size_t i = 1; i < result.reference_pulse_time_ns.size(); ++i) {
        reference_area +=
            0.5 * (result.reference_pulse_amplitude[i - 1] +
                   result.reference_pulse_amplitude[i]) *
            (result.reference_pulse_time_ns[i] -
             result.reference_pulse_time_ns[i - 1]);
    }
    require(std::abs(reference_area - result.single_pe_area_mv_ns) < 1.0e-9,
            "serialized reference pulse must reproduce single-p.e. area");
    require(result.camera_trigger.triggered,
            "voltage threshold must operate on sampled waveform");

    DetectorPipelineConfig boundary_config;
    boundary_config.microcell.enabled = false;
    boundary_config.single_pe.enabled = true;
    boundary_config.single_pe.model = "analytic";
    boundary_config.single_pe.unit = "pe_charge";
    boundary_config.single_pe.support_ns = 30.0;
    boundary_config.sampling.start_ns = 0.0;
    boundary_config.sampling.end_ns = 40.0;
    boundary_config.sampling.width_ns = 4.0;
    const auto generation_window =
        waveformContributingPrimaryWindow(boundary_config);
    require(std::abs(generation_window.start_ns + 30.0) < 1.0e-12 &&
                std::abs(generation_window.end_ns - 40.0) < 1.0e-12,
            "causal pulse support must extend the NSB generation prehistory");
    if (argc > 1) {
        boundary_config.single_pe.model = "measured_csv";
        boundary_config.single_pe.csv_path = argv[1];
        boundary_config.single_pe.template_time_reference = "peak";
        boundary_config.sampling.start_ns = -40.0;
        boundary_config.sampling.end_ns = 220.0;
        const auto measured_window =
            waveformContributingPrimaryWindow(boundary_config);
        require(std::abs(measured_window.start_ns + 220.0) < 1.0e-9 &&
                    std::abs(measured_window.end_ns - 260.0) < 1.0e-9,
                "measured peak-aligned pulse must pad both NSB boundaries");
    }
    std::vector<PrimaryPeHit> boundary_hits = {
        {3, 4, 0, -10.0, 0.0, 0.0, 0.0, 1.0,
         HitOrigin::Nsb, false},
        {3, 4, 0, 10.0, 0.0, 0.0, 0.0, 1.0,
         HitOrigin::Nsb, true},
    };
    const auto boundary_result =
        runDetectorPipeline(boundary_config, 1, boundary_hits);
    require(boundary_result.pixels[0].primary_nsb_pe == 1.0 &&
                boundary_result.pixels[0].fired_nsb_pe == 1.0,
            "NSB padding hits must not inflate the integrated image truth");
    require(boundary_result.fired_hits.size() == 2 &&
                waveformSum(boundary_result, 0) > 1.0,
            "NSB padding hits must still contribute to the stored waveform");

    std::vector<PrimaryPeHit> image_gate_hits = {
        {5, 6, 0, -1.0, 0.0, 0.0, 0.0, 1.0,
         HitOrigin::Nsb, false},
        {5, 6, 0, 0.0, 0.001, 0.0, 0.0, 1.0,
         HitOrigin::Nsb, true},
        {5, 6, 0, 31.999, 0.002, 0.0, 0.0, 1.0,
         HitOrigin::Nsb, true},
        {5, 6, 0, 32.0, 0.003, 0.0, 0.0, 1.0,
         HitOrigin::Nsb, false},
    };
    const auto image_gate_result =
        runDetectorPipeline(boundary_config, 1, image_gate_hits);
    require(image_gate_result.pixels[0].primary_nsb_pe == 2.0 &&
                image_gate_result.pixels[0].fired_nsb_pe == 2.0,
            "electronic image truth must include exactly the configured "
            "32 ns NSB gate");
    require(image_gate_result.fired_hits.size() == 4,
            "NSB outside the image gate must remain available to the waveform");

    config.camera_trigger.enabled = false;
    config.single_pe.charge_fluctuation.enabled = true;
    config.single_pe.charge_fluctuation.empirical_samples = {0.5, 1.5};
    const auto fluctuated = runDetectorPipeline(config, 2, hits);
    const auto repeated_fluctuated = runDetectorPipeline(config, 2, hits);
    require(fluctuated.charge_fluctuation_enabled,
            "result metadata must expose enabled charge fluctuation");
    require(fluctuated.fired_hits.size() == repeated_fluctuated.fired_hits.size(),
            "fixed random seed must preserve avalanche count");
    std::vector<double> charge_factor_sum(2, 0.0);
    for (std::size_t i = 0; i < fluctuated.fired_hits.size(); ++i) {
        const auto& hit = fluctuated.fired_hits[i];
        require(hit.charge_factor == repeated_fluctuated.fired_hits[i].charge_factor,
                "fixed random seed must reproduce empirical charge draws");
        charge_factor_sum[static_cast<std::size_t>(hit.pixel_id)] +=
            hit.charge_factor * hit.fired_pe;
    }
    for (std::size_t pixel = 0; pixel < charge_factor_sum.size(); ++pixel) {
        const double waveform_area =
            waveformSum(fluctuated, pixel) * config.sampling.width_ns;
        const double expected_area =
            charge_factor_sum[pixel] * fluctuated.single_pe_area_mv_ns;
        require(std::abs(waveform_area - expected_area) < 1.0e-6,
                "mV waveform area must follow empirical avalanche charge");
    }

    std::cout << "detector pipeline tests passed\n";
    return 0;
}
