#pragma once

// Run-level accumulation: waveform binning, collector debug rows,
// the emission-altitude histogram and the profiling clock.

#include <map>
#include <vector>

#include "app/OpticalSimCommon.hpp"
#include "core/PhotonBunch.hpp"
#include "io/CorsikaTraceStats.hpp"

namespace lact {

std::vector<AtmosphereHistogramBin> makeAtmosphereHistogramBins(
    const AtmosphereHistogramConfig& cfg);

void accumulateAtmosphereHistogram(std::vector<AtmosphereHistogramBin>& bins,
                                   const AtmosphereHistogramConfig& cfg,
                                   double altitude_km,
                                   double before_weight,
                                   double after_weight,
                                   double theory_weight);

void addElapsed(ProfileStats& stats,
                double ProfileStats::*field,
                const std::chrono::steady_clock::time_point& start);

void accumulateWaveformHit(std::map<WaveformKey, WaveformPixelAccumulator>& waveform,
                           std::vector<RawWaveformHit>& raw_waveform_hits,
                           const WaveformOutputConfig& cfg,
                           bool capture_detector_hit,
                           int event_id,
                           int telescope_id,
                           const OpticalSurfaceHit& hit,
                           PhotonOrigin origin = PhotonOrigin::Cherenkov);

void appendCollectorDebugPhoton(std::vector<CollectorDebugPhotonRow>& rows,
                                const CollectorDebugConfig& cfg,
                                int event_id,
                                int telescope_id,
                                const OpticalSurfaceHit& hit);

} // namespace lact
