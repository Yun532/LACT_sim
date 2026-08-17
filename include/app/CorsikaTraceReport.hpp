#pragma once

// Human-readable console reporting for the CORSIKA trace: configuration
// echo, per-event and per-stream summaries, profiling and the analytic
// pure-NSB trigger estimate. Nothing here affects simulation results.

#include <map>
#include <string>
#include <vector>

#include "app/OpticalSimCommon.hpp"
#include "app/PhotonResponseSampler.hpp"
#include "io/CameraElectronicsEvent.hpp"
#include "io/CorsikaTraceEventMetadata.hpp"
#include "io/CorsikaTraceStats.hpp"

namespace lact {

std::string formatStreamSummaryLine(const TraceSummary& s,
                                    bool camera_enabled,
                                    const std::string& event_id_mode);

std::string formatTelescopeEventLine(const TelescopeEventAccumulator& s,
                                     bool camera_enabled);

void printTraceStreamSummary(const std::map<SummaryKey, TraceSummary>& summaries,
                             bool camera_enabled,
                             const std::string& event_id_mode);

void printEventSummary(const std::map<SummaryKey, TraceSummary>& summaries,
                       bool camera_enabled,
                       const std::string& event_id_mode,
                       const EventIOMetadata& metadata);

void printCorsikaOpticalConfiguration(
    const std::map<std::string, std::string>& cfg,
    const TelescopeConfig& telescope_cfg,
    const TelescopeFrame& source_adapter_frame,
    const MirrorLayout& mirrors,
    const SyntheticPhotonConfig& source_cfg,
    const SourceRuntimeConfig& source_runtime_cfg,
    const OutputPlane& plane,
    const CorsikaTraceOutputConfig& output_cfg,
    const CameraConfig& camera_cfg,
    const CameraGeometry& camera,
    const std::unique_ptr<Cone::SquareCone>& light_collector,
    const SipmConfig& sipm_cfg,
    const ElectronicsConfig& electronics_cfg,
    const electronics::DetectorPipelineConfig& detector_cfg,
    const WaveformOutputConfig& waveform_cfg,
    const CollectorDebugConfig& collector_debug_cfg,
    const NsbConfig& nsb_cfg,
    const TriggerConfig& trigger_cfg,
    const PhotonResponseConfig& response_cfg,
    const OpticalEfficiencyConfig& efficiency_cfg,
    const AtmosphereTransmissionConfig& atmosphere_cfg,
    const ErrorConfig& error_cfg,
    const ObstructionMask& obstruction,
    const PropagationConfig& propagation_cfg,
    const EventIOPhotonConfig& eventio_cfg,
    const std::string& missing_wavelength_range_source,
    double eventio_mirror_front_z_m,
    bool eventio_2d_backproject);

void printProfileStats(const ProfileStats& stats, double trace_time_s);

void printPureNsbTriggerEstimate(const NsbConfig& nsb_cfg,
                                 const WaveformOutputConfig& waveform_cfg,
                                 const TriggerConfig& trigger_cfg,
                                 std::size_t n_pixels,
                                 std::size_t n_telescopes);

bool shouldHideInputCardLine(const std::string& line);

void printEventIOMetadataSummary(const EventIOMetadata& metadata);

} // namespace lact
