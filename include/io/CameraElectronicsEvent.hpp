#pragma once

#include <map>
#include <limits>
#include <vector>

#include "app/OpticalSimCommon.hpp"
#include "io/CorsikaTraceOutputTypes.hpp"

namespace lact {

// Canonical event/telescope electronics product.  All serializers consume
// this object so that ROOT, HDF5 and CSV differ only in representation, not in
// the detector response, NSB realization or camera trigger decision.
struct CameraElectronicsEvent {
    int event_id = 0;
    int telescope_id = 0;
    double reference_time_ns = 0.0;
    electronics::DetectorPipelineResult detector;

    // Absolute EventIO-time trigger fields.  The detector result above keeps
    // its natural time coordinates relative to reference_time_ns.
    double trigger_time_ns = std::numeric_limits<double>::quiet_NaN();
    double geometric_delay_ns = std::numeric_limits<double>::quiet_NaN();
    double coincidence_time_ns = std::numeric_limits<double>::quiet_NaN();
    bool array_triggered = false;
    bool selected_for_output = true;
};

using CameraElectronicsEventMap =
    std::map<SummaryKey, CameraElectronicsEvent>;

double waveformReferenceTimeNs(const WaveformOutputConfig& waveform,
                               const TraceSummary& summary);

CameraElectronicsEventMap buildCameraElectronicsEvents(
    const electronics::DetectorPipelineConfig& detector,
    const WaveformOutputConfig& waveform,
    const NsbConfig& nsb,
    const std::vector<int>& pixel_id_axis,
    const std::map<SummaryKey, TraceSummary>& summaries,
    const std::vector<RawWaveformHit>& raw_hits);

} // namespace lact
