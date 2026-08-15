#pragma once

// CSV serialisation of the CORSIKA trace: whiteboard hits, mirror
// diagnostics, per-pixel images, summaries and the electronics tables.

#include <fstream>
#include <string>

#include "app/OpticalSimCommon.hpp"
#include "io/CameraElectronicsEvent.hpp"
#include "io/CorsikaTraceEventMetadata.hpp"
#include "io/CorsikaTraceStats.hpp"

namespace lact {

void writeAtmosphereHistogramCsv(const AtmosphereHistogramConfig& cfg,
                                 const std::vector<AtmosphereHistogramBin>& bins);

void writeCollectorDebugCsv(const CollectorDebugConfig& cfg,
                            const std::vector<CollectorDebugPhotonRow>& rows);

void writeCorsikaWhiteboardHeader(std::ofstream& ofs, bool include_emitter_info);

void writeCorsikaWhiteboardHit(std::ofstream& ofs,
                               const PhotonBunch& bunch,
                               std::uint64_t photon_index,
                               const OpticalSurfaceHit& hit,
                               bool include_emitter_info);

void writeCorsikaMirrorDiagnosticHeader(std::ofstream& ofs);

void writeCorsikaMirrorDiagnosticHit(std::ofstream& ofs,
                                     const PhotonBunch& bunch,
                                     std::uint64_t photon_index,
                                     const OpticalSurfaceHit& hit,
                                     const char* status);

WhiteboardHdf5Row makeWhiteboardHdf5Row(const PhotonBunch& bunch,
                                        std::uint64_t photon_index,
                                        const OpticalSurfaceHit& hit);

void writeSummaryCsv(const std::string& path,
                     const std::map<SummaryKey, TraceSummary>& summaries,
                     const CameraElectronicsEventMap& electronics_events);

void ensureCsvParent(const std::string& path);

void writeElectronicsPixelCsv(
    const std::string& path,
    const std::vector<int>& pixel_id_axis,
    const CameraElectronicsEventMap& events,
    const std::map<PixelKey, PixelAccumulator>& optical_pixels);

void writeElectronicsWaveformCsv(
    const std::string& path,
    const std::vector<int>& pixel_id_axis,
    const CameraElectronicsEventMap& events);

void writeElectronicsTriggerCsv(const std::string& path,
                                const CameraElectronicsEventMap& events);

void writeElectronicsHitCsv(const std::string& primary_path,
                            const std::string& fired_path,
                            const std::vector<int>& pixel_id_axis,
                            const CameraElectronicsEventMap& events,
                            bool write_primary,
                            bool write_fired);

} // namespace lact
