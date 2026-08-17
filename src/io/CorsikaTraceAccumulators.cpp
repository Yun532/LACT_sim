#include "io/CorsikaTraceAccumulators.hpp"

#include <algorithm>
#include <cmath>

namespace lact {

std::vector<AtmosphereHistogramBin> makeAtmosphereHistogramBins(
    const AtmosphereHistogramConfig& cfg)
{
    std::vector<AtmosphereHistogramBin> bins;
    if (!cfg.enabled) {
        return bins;
    }
    for (double low = cfg.min_altitude_km; low < cfg.max_altitude_km;
         low += cfg.bin_width_km) {
        bins.push_back({low, std::min(cfg.max_altitude_km, low + cfg.bin_width_km)});
    }
    return bins;
}

void accumulateAtmosphereHistogram(std::vector<AtmosphereHistogramBin>& bins,
                                   const AtmosphereHistogramConfig& cfg,
                                   double altitude_km,
                                   double before_weight,
                                   double after_weight,
                                   double theory_weight)
{
    if (!cfg.enabled || bins.empty() || !std::isfinite(altitude_km)) {
        return;
    }
    if (altitude_km < cfg.min_altitude_km || altitude_km >= cfg.max_altitude_km) {
        return;
    }
    std::size_t index = static_cast<std::size_t>(
        (altitude_km - cfg.min_altitude_km) / cfg.bin_width_km);
    if (index >= bins.size()) {
        index = bins.size() - 1;
    }
    auto& bin = bins[index];
    bin.bunches += 1;
    bin.before_weight += before_weight;
    bin.after_weight += after_weight;
    bin.theory_weight += theory_weight;
}

void addElapsed(ProfileStats& stats,
                double ProfileStats::*field,
                const std::chrono::steady_clock::time_point& start)
{
    stats.*field +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}



void accumulateWaveformHit(std::map<WaveformKey, WaveformPixelAccumulator>& waveform,
                           std::vector<RawWaveformHit>& raw_waveform_hits,
                           const WaveformOutputConfig& cfg,
                           bool capture_detector_hit,
                           int event_id,
                           int telescope_id,
                           const OpticalSurfaceHit& hit)
{
    if (!hit.hit_camera || hit.pixel_id < 0) {
        return;
    }
    const double pe = hit.weight * hit.relative_efficiency;
    if (capture_detector_hit) {
        raw_waveform_hits.push_back(RawWaveformHit{
            event_id,
            telescope_id,
            hit.pixel_id,
            hit.time_ns,
            1,
            pe,
            hit.collector_exit_x_m,
            hit.collector_exit_y_m,
            hit.wavelength_nm,
            electronics::HitOrigin::Cherenkov,
        });
    }
    if (!cfg.enabled) {
        return;
    }
    if (waveformUsesImageReference(cfg)) {
        if (!capture_detector_hit) {
            raw_waveform_hits.push_back(RawWaveformHit{
                event_id,
                telescope_id,
                hit.pixel_id,
                hit.time_ns,
                1,
                pe,
                hit.collector_exit_x_m,
                hit.collector_exit_y_m,
                hit.wavelength_nm,
                electronics::HitOrigin::Cherenkov,
            });
        }
        return;
    }
    const int bin = waveformBinForTime(cfg, hit.time_ns);
    if (bin < 0) {
        return;
    }
    auto& acc = waveform[{event_id, telescope_id, hit.pixel_id, bin}];
    acc.event_id = event_id;
    acc.telescope_id = telescope_id;
    acc.pixel_id = hit.pixel_id;
    acc.time_bin = bin;
    acc.photon_count += 1;
    acc.pe += pe;
}

void appendCollectorDebugPhoton(std::vector<CollectorDebugPhotonRow>& rows,
                                const CollectorDebugConfig& cfg,
                                int event_id,
                                int telescope_id,
                                const OpticalSurfaceHit& hit)
{
    if (!cfg.photon_output || !hit.collector_enabled ||
        rows.size() >= cfg.max_photons) {
        return;
    }
    rows.push_back(CollectorDebugPhotonRow{
        event_id,
        telescope_id,
        hit.pixel_id,
        hit.accepted ? 1 : 0,
        hit.collector_reflections,
        hit.collector_reflection_limit_reached ? 1 : 0,
        hit.time_ns,
        hit.collector_path_length_m,
        hit.collector_time_delay_ns,
        hit.wavelength_nm,
        hit.weight,
        hit.relative_efficiency,
        hit.weight * hit.relative_efficiency,
        hit.collector_exit_x_m,
        hit.collector_exit_y_m,
        hit.collector_exit_z_m,
        hit.collector_dir_u,
        hit.collector_dir_v,
        hit.collector_dir_w,
        hit.sipm_gap_rejected ? 1 : 0,
        hit.sipm_channel_gap_rejected ? 1 : 0,
        hit.sipm_grid_column,
        hit.sipm_grid_row,
        hit.sipm_channel_id,
        hit.sipm_microcell_id,
    });
}

} // namespace lact
