#pragma once

#include <map>
#include <string>
#include <vector>

#include "app/OpticalSimCommon.hpp"
#include "io/CorsikaTraceOutputTypes.hpp"
#include "io/EventIOPhotonSource.hpp"

namespace lact {

void writeLactEventRoot(const CorsikaTraceOutputConfig& output_cfg,
                        const WaveformOutputConfig& waveform_cfg,
                        const std::string& main_config_path,
                        const std::map<std::string, std::string>& cfg,
                        const SourceRuntimeConfig& source_runtime_cfg,
                        const TelescopeConfig& telescope_cfg,
                        const EventIOMetadata& metadata,
                        const CameraGeometry& camera,
                        const std::vector<MirrorFacet>& facets,
                        const NsbConfig& nsb_cfg,
                        const TriggerConfig& trigger_cfg,
                        const std::map<SummaryKey, TraceSummary>& summaries,
                        const std::map<PixelKey, PixelAccumulator>& pixels,
                        const std::map<WaveformKey, WaveformPixelAccumulator>& waveforms,
                        const std::vector<RawWaveformHit>& raw_waveform_hits);

} // namespace lact
