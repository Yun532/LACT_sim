#include "io/CameraElectronicsEvent.hpp"

#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <unordered_map>

namespace lact {

double waveformReferenceTimeNs(const WaveformOutputConfig& waveform,
                               const TraceSummary& summary)
{
    if (waveform.time_reference == "image_first" &&
        std::isfinite(summary.first_cherenkov_time_ns)) {
        return summary.first_cherenkov_time_ns;
    }
    if (waveform.time_reference == "image_mean" &&
        summary.weighted_signal > 0.0) {
        return summary.weighted_time_sum / summary.weighted_signal;
    }
    return 0.0;
}

CameraElectronicsEventMap buildCameraElectronicsEvents(
    const electronics::DetectorPipelineConfig& detector,
    const WaveformOutputConfig& waveform,
    const NsbConfig& nsb,
    const std::vector<int>& pixel_id_axis,
    const std::map<SummaryKey, TraceSummary>& summaries,
    const std::vector<RawWaveformHit>& raw_hits)
{
    CameraElectronicsEventMap output;
    if (!detector.enabled || pixel_id_axis.empty()) {
        return output;
    }

    std::unordered_map<int, std::size_t> pixel_to_column;
    pixel_to_column.reserve(pixel_id_axis.size());
    for (std::size_t column = 0; column < pixel_id_axis.size(); ++column) {
        pixel_to_column.emplace(pixel_id_axis[column], column);
    }

    std::set<SummaryKey> keys;
    for (const auto& item : summaries) {
        keys.insert(item.first);
    }
    for (const auto& hit : raw_hits) {
        keys.insert({hit.event_id, hit.telescope_id});
    }

    for (const auto& key : keys) {
        const auto summary_it = summaries.find(key);
        TraceSummary empty_summary;
        empty_summary.event_id = key.first;
        empty_summary.telescope_id = key.second;
        const TraceSummary& summary =
            summary_it == summaries.end() ? empty_summary : summary_it->second;

        CameraElectronicsEvent event;
        event.event_id = key.first;
        event.telescope_id = key.second;
        event.reference_time_ns = waveformReferenceTimeNs(waveform, summary);

        std::vector<electronics::PrimaryPeHit> primary_hits;
        for (const auto& hit : raw_hits) {
            if (hit.event_id != key.first || hit.telescope_id != key.second) {
                continue;
            }
            const auto column = pixel_to_column.find(hit.pixel_id);
            if (column == pixel_to_column.end()) {
                continue;
            }
            primary_hits.push_back({
                hit.event_id,
                hit.telescope_id,
                static_cast<int>(column->second),
                hit.time_ns - event.reference_time_ns,
                hit.sensor_x_m,
                hit.sensor_y_m,
                hit.wavelength_nm,
                hit.pe,
                hit.origin,
            });
        }

        if (nsb.enabled && nsb.rate_pe_per_ns_per_pixel > 0.0) {
            const auto generation_window =
                electronics::waveformContributingPrimaryWindow(detector);
            auto nsb_hits = electronics::generateUniformNsbPrimaryHits(
                key.first,
                key.second,
                pixel_id_axis.size(),
                nsb.rate_pe_per_ns_per_pixel,
                generation_window.start_ns,
                generation_window.end_ns,
                detector.microcell,
                nsb.seed);
            for (auto& hit : nsb_hits) {
                hit.count_in_integrated_image =
                    hit.time_ns >= detector.sampling.start_ns &&
                    hit.time_ns < detector.sampling.end_ns;
            }
            primary_hits.insert(primary_hits.end(),
                                nsb_hits.begin(), nsb_hits.end());
        }

        event.detector = electronics::runDetectorPipeline(
            detector, pixel_id_axis.size(), primary_hits);
        if (event.detector.camera_trigger.triggered) {
            event.trigger_time_ns =
                event.reference_time_ns +
                event.detector.camera_trigger.trigger_time_ns;
        }
        output.emplace(key, std::move(event));
    }
    return output;
}

} // namespace lact
