#pragma once

// Native .h5 output for the CORSIKA trace. Compiled only when the build found
// an HDF5 C library; the app checks outputWantsHdf5() before calling in.

#ifdef LACT_HAS_HDF5

#include <map>
#include <string>
#include <vector>

#include "app/OpticalSimCommon.hpp"
#include "io/CameraElectronicsEvent.hpp"
#include "io/CorsikaTraceEventMetadata.hpp"
#include "io/CorsikaTraceOutputTypes.hpp"

namespace lact {

void writeNativeTraceHdf5(const CorsikaTraceOutputConfig& output_cfg,
                          const WaveformOutputConfig& waveform_cfg,
                          const std::string& main_config_path,
                          const std::map<std::string, std::string>& cfg,
                          const ComponentConfigPaths& component_paths,
                          const SourceRuntimeConfig& source_runtime_cfg,
                          const TelescopeConfig& telescope_cfg,
                          const EventIOMetadata& metadata,
                          const CameraGeometry& camera,
                          const std::vector<MirrorFacet>& facets,
                          const SipmConfig& sipm_cfg,
                          const ElectronicsConfig& electronics_cfg,
                          const electronics::DetectorPipelineConfig& detector_cfg,
                          const OpticalEfficiencyConfig& efficiency_cfg,
                          const NsbConfig& nsb_cfg,
                          const TriggerConfig& trigger_cfg,
                          const std::map<SummaryKey, TraceSummary>& summaries,
                          const std::map<PixelKey, PixelAccumulator>& pixels,
                          const std::map<WaveformKey, WaveformPixelAccumulator>& waveforms,
                          const std::vector<RawWaveformHit>& raw_waveform_hits,
                          const CameraElectronicsEventMap& electronics_events,
                          const std::vector<WhiteboardHdf5Row>& whiteboard_hits);

} // namespace lact

#endif // LACT_HAS_HDF5
