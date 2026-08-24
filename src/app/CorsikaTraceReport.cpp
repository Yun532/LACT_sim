#include "app/CorsikaTraceReport.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace lact {

std::string formatStreamSummaryLine(const TraceSummary& s,
                                    bool camera_enabled,
                                    const std::string& event_id_mode)
{
    const int array_id = arrayIdFromOutputEvent(s.event_id, event_id_mode);
    const double mean = s.weighted_signal > 0.0
        ? s.weighted_time_sum / s.weighted_signal
        : 0.0;
    const double var = s.weighted_signal > 0.0
        ? std::max(0.0, s.weighted_time2_sum / s.weighted_signal - mean * mean)
        : 0.0;
    const double output_transmission = s.hit_output_before_obstruction > 0
        ? static_cast<double>(s.hit_output_plane) /
          static_cast<double>(s.hit_output_before_obstruction)
        : 0.0;
    std::ostringstream value;
    value << "array_id=" << array_id
          << " telescope=" << s.telescope_id
          << " event_id=" << s.event_id
          << " input_bunches=" << s.input_bunches
          << " input_photons=" << doubleToString(s.input_photons, 3)
          << " blocked=" << s.blocked_by_obstruction
          << " blocked_incoming=" << s.blocked_incoming
          << " blocked_reflected=" << s.blocked_reflected
          << " hit_mirror_before_obstruction=" << s.hit_mirror_before_obstruction
          << " hit_output_before_obstruction=" << s.hit_output_before_obstruction
          << " hit_mirror=" << s.hit_mirror
          << " hit_output=" << s.hit_output_plane
          << " output_transmission_after_obstruction="
          << doubleToString(output_transmission, 6);
    if (camera_enabled) {
        value << " hit_camera=" << s.hit_camera
              << " accepted=" << s.accepted_camera
              << " lost=" << s.lost_between_pixels
              << " unique_pixels=" << s.unique_pixels.size();
    }
    value << (camera_enabled ? " pe=" : " signal=")
          << doubleToString(s.weighted_signal, 3)
          << " time_mean_ns=" << doubleToString(mean, 3)
          << " time_rms_ns=" << doubleToString(std::sqrt(var), 3)
          << " time_first_ns="
          << (std::isfinite(s.first_cherenkov_time_ns)
                  ? doubleToString(s.first_cherenkov_time_ns, 3)
                  : "nan");
    return value.str();
}

std::string formatTelescopeEventLine(const TelescopeEventAccumulator& s,
                                     bool camera_enabled)
{
    const double mean = s.weighted_signal > 0.0
        ? s.weighted_time_sum / s.weighted_signal
        : 0.0;
    const double var = s.weighted_signal > 0.0
        ? std::max(0.0, s.weighted_time2_sum / s.weighted_signal - mean * mean)
        : 0.0;
    const double output_transmission = s.hit_output_before_obstruction > 0
        ? static_cast<double>(s.hit_output_plane) /
          static_cast<double>(s.hit_output_before_obstruction)
        : 0.0;
    std::ostringstream value;
    value << "telescope=" << s.telescope_id
          << " output_events=" << s.output_events.size()
          << " input_bunches=" << s.input_bunches
          << " input_photons=" << doubleToString(s.input_photons, 3)
          << " blocked=" << s.blocked_by_obstruction
          << " blocked_incoming=" << s.blocked_incoming
          << " blocked_reflected=" << s.blocked_reflected
          << " hit_mirror_before_obstruction=" << s.hit_mirror_before_obstruction
          << " hit_output_before_obstruction=" << s.hit_output_before_obstruction
          << " hit_mirror=" << s.hit_mirror
          << " hit_output=" << s.hit_output_plane
          << " output_transmission_after_obstruction="
          << doubleToString(output_transmission, 6);
    if (camera_enabled) {
        value << " hit_camera=" << s.hit_camera
              << " accepted=" << s.accepted_camera
              << " lost=" << s.lost_between_pixels
              << " unique_pixels=" << s.unique_pixels.size();
    }
    value << (camera_enabled ? " pe=" : " signal=")
          << doubleToString(s.weighted_signal, 3)
          << " time_mean_ns=" << doubleToString(mean, 3)
          << " time_rms_ns=" << doubleToString(std::sqrt(var), 3)
          << " time_first_ns="
          << (std::isfinite(s.first_cherenkov_time_ns)
                  ? doubleToString(s.first_cherenkov_time_ns, 3)
                  : "nan");
    return value.str();
}

void printTraceStreamSummary(const std::map<SummaryKey, TraceSummary>& summaries,
                             bool camera_enabled,
                             const std::string& event_id_mode)
{
    printSection("Per-stream summary");
    printField("columns",
               camera_enabled
                   ? "event_id shower_event array_id telescope input_bunches input_photons hit_output hit_camera accepted unique_pixels pe time_mean_ns time_rms_ns time_first_ns"
                   : "event_id shower_event array_id telescope input_bunches input_photons hit_output signal time_mean_ns time_rms_ns time_first_ns");

    for (const auto& kv : summaries) {
        const auto& s = kv.second;
        const int shower_event = showerEventFromOutputEvent(s.event_id, event_id_mode);
        std::ostringstream value;
        value << "shower_event=" << shower_event << " "
              << formatStreamSummaryLine(s, camera_enabled, event_id_mode);
        printField("stream", value.str());
    }
}

void printEventSummary(const std::map<SummaryKey, TraceSummary>& summaries,
                       bool camera_enabled,
                       const std::string& event_id_mode,
                       const EventIOMetadata& metadata)
{
    struct EventAggregate {
        int event_id = 0;
        int shower_event = 0;
        int array_id = 0;
        std::set<int> output_events;
        std::set<int> telescopes;
        std::uint64_t input_bunches = 0;
        double input_photons = 0.0;
        std::uint64_t blocked_by_obstruction = 0;
        std::uint64_t blocked_incoming = 0;
        std::uint64_t blocked_reflected = 0;
        std::uint64_t hit_mirror_before_obstruction = 0;
        std::uint64_t hit_output_before_obstruction = 0;
        std::uint64_t hit_output_plane = 0;
        std::uint64_t hit_camera = 0;
        std::uint64_t accepted_camera = 0;
        double weighted_signal = 0.0;
        double weighted_time_sum = 0.0;
        double weighted_time2_sum = 0.0;
        double first_cherenkov_time_ns = std::numeric_limits<double>::infinity();
    };

    std::map<int, EventAggregate> events;
    for (const auto& kv : summaries) {
        const auto& s = kv.second;
        const auto identity = outputEventMetadata(s.event_id, event_id_mode, metadata);
        auto& e = events[s.event_id];
        e.event_id = s.event_id;
        e.shower_event = identity.shower_event;
        e.array_id = identity.array_id;
        e.output_events.insert(s.event_id);
        e.telescopes.insert(s.telescope_id);
        e.input_bunches += s.input_bunches;
        e.input_photons += s.input_photons;
        e.blocked_by_obstruction += s.blocked_by_obstruction;
        e.blocked_incoming += s.blocked_incoming;
        e.blocked_reflected += s.blocked_reflected;
        e.hit_mirror_before_obstruction += s.hit_mirror_before_obstruction;
        e.hit_output_before_obstruction += s.hit_output_before_obstruction;
        e.hit_output_plane += s.hit_output_plane;
        e.hit_camera += s.hit_camera;
        e.accepted_camera += s.accepted_camera;
        e.weighted_signal += s.weighted_signal;
        e.weighted_time_sum += s.weighted_time_sum;
        e.weighted_time2_sum += s.weighted_time2_sum;
        e.first_cherenkov_time_ns =
            std::min(e.first_cherenkov_time_ns, s.first_cherenkov_time_ns);
    }

    printSection("Per-event summary");
    printField("columns",
               camera_enabled
                   ? "event_id shower_event array_id energy_gev core_N_m core_E_m core_source telescopes input_bunches input_photons blocked hit_output hit_camera accepted pe time_mean_ns time_rms_ns time_first_ns"
                   : "event_id shower_event array_id energy_gev core_N_m core_E_m core_source telescopes input_bunches input_photons blocked hit_output signal time_mean_ns time_rms_ns time_first_ns");
    for (const auto& kv : events) {
        const auto& e = kv.second;
        const OutputEventMetadata meta = outputEventMetadata(e.event_id, event_id_mode, metadata);
        const double mean = e.weighted_signal > 0.0
            ? e.weighted_time_sum / e.weighted_signal
            : 0.0;
        const double var = e.weighted_signal > 0.0
            ? std::max(0.0, e.weighted_time2_sum / e.weighted_signal - mean * mean)
            : 0.0;
        std::ostringstream value;
        value << "event_id=" << e.event_id
              << " shower_event=" << e.shower_event
              << " array_id=" << e.array_id;
        if (meta.found) {
            value << " energy_gev=" << doubleToString(meta.energy_gev, 6)
                  << " core_N_m=" << doubleToString(meta.core_x_north_m, 3)
                  << " core_E_m=" << doubleToString(-meta.core_y_west_m, 3)
                  << " core_source="
                  << (meta.used_array_offset ? "negative_MC_TELOFF_array_offset"
                                              : "shower_header");
        } else {
            value << " energy_gev=unknown core_N_m=unknown core_E_m=unknown"
                  << " core_source=missing_metadata";
        }
        value << " output_events=" << e.output_events.size()
              << " telescopes=" << e.telescopes.size()
              << " input_bunches=" << e.input_bunches
              << " input_photons=" << doubleToString(e.input_photons, 3)
              << " blocked=" << e.blocked_by_obstruction
              << " blocked_incoming=" << e.blocked_incoming
              << " blocked_reflected=" << e.blocked_reflected
              << " hit_mirror_before_obstruction=" << e.hit_mirror_before_obstruction
              << " hit_output_before_obstruction=" << e.hit_output_before_obstruction
              << " hit_output=" << e.hit_output_plane;
        if (e.hit_output_before_obstruction > 0) {
            value << " output_transmission_after_obstruction="
                  << doubleToString(static_cast<double>(e.hit_output_plane) /
                                    static_cast<double>(e.hit_output_before_obstruction), 6);
        }
        if (camera_enabled) {
            value << " hit_camera=" << e.hit_camera
                  << " accepted=" << e.accepted_camera;
        }
        value << (camera_enabled ? " pe=" : " signal=")
              << doubleToString(e.weighted_signal, 3)
              << " time_mean_ns=" << doubleToString(mean, 3)
              << " time_rms_ns=" << doubleToString(std::sqrt(var), 3)
              << " time_first_ns="
              << (std::isfinite(e.first_cherenkov_time_ns)
                      ? doubleToString(e.first_cherenkov_time_ns, 3)
                      : "nan");
        printField("event", value.str());
    }
}

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
    bool eventio_2d_backproject)
{
    const std::string mirror_mode = lowerCopy(getString(cfg, "mirror.mode", "generated"));
    const std::string mirror_csv = getString(cfg, "mirror.csv_path", "");

    printSection("Telescope");
    printField("id", intToString(static_cast<std::uint64_t>(telescope_cfg.id)));
    printField("name", telescope_cfg.name);
    printField("position_m", vec3ToString(telescope_cfg.position_m));
    printField("pointing_az_deg", doubleToString(telescope_cfg.pointing_az_deg));
    printField("pointing_el_deg", doubleToString(telescope_cfg.pointing_el_deg));
    printField("focal_length_m", doubleToString(telescope_cfg.focal_length_m));
    printField("coordinate_system", telescope_cfg.coordinate_system);
    printField("trace_geometry_frame", "telescope-local; mirror/output geometry is not globally rotated");
    printField("source_adapter_transform", "input frame -> telescope-local by projection onto the axes below");
    printField("local_x_in_input_frame", vec3ToString(source_adapter_frame.x_axis));
    printField("local_y_in_input_frame", vec3ToString(source_adapter_frame.y_axis));
    printField("local_z_boresight_in_input_frame", vec3ToString(source_adapter_frame.z_axis));

    printSection("Mirror");
    printField("mode", mirror_mode);
    if (!mirror_csv.empty()) {
        printField("csv", mirror_csv);
    }
    if (mirror_mode == "elevation_series" || mirror_mode == "series") {
        printField("series_elevation_deg",
                   doubleToString(getDouble(cfg, "mirror.series_elevation_deg",
                                            telescope_cfg.pointing_el_deg)));
        printField("series_angles_deg", getString(cfg, "mirror.series_angles_deg", ""));
        printField("series_csv_pattern", getString(cfg, "mirror.series_csv_pattern", ""));
    }
    printField("facets", intToString(mirrors.size()));

    printSection("Source");
    printField("mode", source_runtime_cfg.use_photon_csv ? "PhotonCsv" : "EventIO");
    if (source_runtime_cfg.use_photon_csv) {
        printField("csv_path", source_runtime_cfg.csv_path);
    } else {
        printField("eventio_path", source_runtime_cfg.eventio_path);
    }
    printField("event_id_mode", source_runtime_cfg.event_id_mode);
    printField("source_coordinate_frame", source_runtime_cfg.coordinate_frame);
    printField("coordinate_interpretation",
               sourceCoordinateFrameDescription(source_runtime_cfg.coordinate_frame));
    printField("rotation_center_local_m",
               vec3ToString(source_runtime_cfg.eventio_rotation_center_local_m));
    printField("eventio_reference_z_m_compat",
               doubleToString(source_runtime_cfg.eventio_reference_z_m));
    printField("eventio_2d_plane_mode", source_runtime_cfg.eventio_2d_plane_mode);
    printField("eventio_mirror_front_z_m", doubleToString(eventio_mirror_front_z_m));
    printField("eventio_2d_trace_direction",
               eventio_2d_backproject ? "signed_line_to_mirror_then_reflect"
                                      : "forward_to_mirror_then_reflect");
    const std::string normalized_frame = normalizeSourceCoordinateFrame(
        source_runtime_cfg.coordinate_frame);
    const bool position_is_applied = normalized_frame == "corsika_nwu_global" ||
                                     normalized_frame == "lact_generic_global";
    printField("eventio_telescope_position",
               position_is_applied
                   ? (source_runtime_cfg.use_eventio_telescope_position
                          ? "subtract EventIO telescope position from global photon positions"
                          : "subtract telescope.position_m from global photon positions")
                   : "metadata only; local/relative photon positions are not shifted");
    printField("filter_telescope_id",
               source_runtime_cfg.filter_telescope_id
                   ? intToString(source_runtime_cfg.selected_telescope_id)
                   : "off");
    printField("filter_event_id",
               source_runtime_cfg.filter_event_id
                   ? intToString(source_runtime_cfg.selected_event_id)
                   : "off");
    printField("filter_event_ids",
               source_runtime_cfg.selected_event_ids.empty()
                   ? "off"
                   : (intToString(source_runtime_cfg.selected_event_ids.size()) +
                      " selected"));
    printField("filter_shower_event_id",
               source_runtime_cfg.filter_shower_event_id
                   ? intToString(source_runtime_cfg.selected_shower_event_id)
                   : "off");
    printField("max_shower_events",
               source_runtime_cfg.max_shower_events > 0
                   ? intToString(source_runtime_cfg.max_shower_events)
                   : "off");
    printField("default_wavelength_nm", doubleToString(source_cfg.wavelength_nm));
    printField("missing_wavelength_model", eventio_cfg.missing_wavelength_model);
    const std::string wavelength_model = lowerCopy(trim(eventio_cfg.missing_wavelength_model));
    if (wavelength_model != "default" &&
        wavelength_model != "fixed" &&
        wavelength_model != "constant" &&
        wavelength_model != "none" &&
        wavelength_model != "off") {
        printField("missing_wavelength_range_nm",
                   doubleToString(eventio_cfg.missing_wavelength_min_nm) + " .. " +
                       doubleToString(eventio_cfg.missing_wavelength_max_nm));
        printField("missing_wavelength_range_source", missing_wavelength_range_source);
        printField("missing_wavelength_seed", intToString(eventio_cfg.missing_wavelength_seed));
    }
    printField("default_weight", doubleToString(source_cfg.photon_weight));
    printField("default_multiplicity", doubleToString(source_cfg.multiplicity));
    printField("read_emitter_info", eventio_cfg.read_emitter_info ? "true" : "false");

    printSection("Output plane");
    printField("point", vec3ToString(plane.point));
    printField("normal", vec3ToString(plane.normal));
    printField("format", output_cfg.format);
    if (outputWantsHdf5(output_cfg)) {
        printField("hdf5_path", output_cfg.hdf5_path);
        printField("hdf5_storage", output_cfg.hdf5_storage);
        printField("hdf5_write_components",
                   output_cfg.hdf5_write_components ? "true" : "false");
        printField("hdf5_write_waveforms",
                   output_cfg.hdf5_write_waveforms ? "true" : "false");
        printField("hdf5_waveform_storage",
                   output_cfg.hdf5_waveform_storage);
        printField("lact_root_write_components",
                   output_cfg.lact_root_write_components ? "true" : "false");
        printField("save_only_triggered",
                   output_cfg.save_only_triggered ? "true" : "false");
        printField("write_pixel_time_stats",
                   output_cfg.write_pixel_time_stats ? "true" : "false");
    }
    if (camera_cfg.enabled && outputWantsCsv(output_cfg)) {
        printField("pixel_csv", output_cfg.pixel_csv);
    } else if (!camera_cfg.enabled && outputWantsCsv(output_cfg)) {
        printField("hits_csv", output_cfg.hits_csv);
    }
    if (!camera_cfg.enabled) {
        printField("whiteboard_emitter_info",
                   output_cfg.whiteboard_emitter_info ? "true" : "false");
    }
    if (outputWantsCsv(output_cfg)) {
        printField("summary_csv", output_cfg.summary_csv);
    }
    if (outputWantsLactRoot(output_cfg)) {
        printField("lact_root_path", output_cfg.lact_root_path);
        printField("lact_profile", output_cfg.lact_profile);
        printField("lact_root_auto_flush_mb",
                   doubleToString(output_cfg.lact_root_auto_flush_mb));
        printField("lact_root_flush_events",
                   output_cfg.lact_root_flush_events > 0
                       ? intToString(output_cfg.lact_root_flush_events)
                       : "off");
    }

    printSection("Camera");
    printField("enabled", camera_cfg.enabled ? "true" : "false");
    if (camera_cfg.enabled) {
        printField("mode", camera_cfg.mode);
        if (!camera_cfg.csv_path.empty()) {
            printField("csv", camera_cfg.csv_path);
        }
        printField("pixels", intToString(static_cast<std::uint64_t>(camera.size())));
        if (!camera.empty()) {
            double min_size = std::numeric_limits<double>::max();
            double max_size = 0.0;
            for (const auto& pixel : camera.pixels()) {
                min_size = std::min(min_size, pixel.size);
                max_size = std::max(max_size, pixel.size);
            }
            printField("pixel_shape", pixelShapeName(camera.pixels().front().shape));
            printField("pixel_size_range_m",
                       doubleToString(min_size) + " .. " + doubleToString(max_size));
        } else {
            printField("pixel_shape", camera_cfg.pixel_shape);
            printField("pixel_size_m", doubleToString(camera_cfg.pixel_size_m));
        }
        if (lowerCopy(trim(camera_cfg.mode)) != "csv") {
            printField("pixel_pitch_m", doubleToString(camera_cfg.pixel_pitch_m));
            printField("radius_m", doubleToString(camera_cfg.radius_m));
        }
        printField("coordinates", "output-plane local u/v");
        printField("collector", light_collector ? camera_cfg.collector
                                                : "not set -> direct pixel containment");
        if (light_collector) {
            printField("collector_material", camera_cfg.collector_material);
            if (!isDisabledText(camera_cfg.collector_reflectivity_csv)) {
                printField("collector_reflectivity_csv",
                           camera_cfg.collector_reflectivity_csv);
            }
            printField("collector_entrance_m",
                       doubleToString(cameraPixelSizeForCollector(camera_cfg, camera)));
            printField("collector_exit_m", doubleToString(camera_cfg.collector_exit_size_m));
            printField("collector_height_m", doubleToString(camera_cfg.collector_height_m));
            printField("collector_debug_photon_output",
                       collector_debug_cfg.photon_output ? "true" : "false");
            if (collector_debug_cfg.photon_output) {
                printField("collector_debug_photon_csv", collector_debug_cfg.photon_csv);
                printField("collector_debug_max_photons",
                           intToString(collector_debug_cfg.max_photons));
            }
        }
    } else {
        printField("mode", camera_cfg.whiteboard
                               ? "whiteboard"
                               : "implicit_whiteboard_legacy");
    }

    printSection("SiPM");
    printField("size_m", doubleToString(sipm_cfg.size_m));
    printField("pde", factorDescription(efficiency_cfg.sipm_pde));

    printSection("Electronics");
    printField("electronics_enabled", detector_cfg.enabled ? "true" : "false");
    printField("pde_stage", "upstream sipm.pde, applied exactly once");
    printField("microcell_enabled",
               detector_cfg.microcell.enabled ? "true" : "false");
    printField("microcell_saturation_enabled",
               detector_cfg.microcell.saturation_enabled ? "true" : "false");
    printField("microcell_recovery_enabled",
               detector_cfg.microcell.recovery_enabled ? "true" : "false");
    printField("microcell_recovery_time_ns",
               doubleToString(detector_cfg.microcell.recovery_time_ns));
    printField("microcell_model", detector_cfg.microcell.model);
    printField("microcell_layout", detector_cfg.microcell.layout);
    printField("inter_channel_gap_geometry", "part of tiled layout");
    printField(
        "pde_includes_inter_channel_gaps",
        detector_cfg.microcell.pde_includes_inter_channel_gaps ? "true"
                                                               : "false");
    printField(
        "inter_channel_active_fraction",
        doubleToString(
            ::lact::electronics::interChannelActiveFraction(
                detector_cfg.microcell),
            9));
    printField(
        "microcell_grid",
        intToString(detector_cfg.microcell.grid_columns) + " x " +
            intToString(detector_cfg.microcell.grid_rows));
    printField("microcells_per_pixel",
               intToString(
                   static_cast<std::uint64_t>(
                       detector_cfg.microcell.channels_per_pixel) *
                   static_cast<std::uint64_t>(
                       detector_cfg.microcell.microcells_per_channel)));
    printField("channels_per_pixel",
               intToString(detector_cfg.microcell.channels_per_pixel));
    printField("channel_merge", "direct sum");
    printField("single_pe_enabled",
               detector_cfg.single_pe.enabled ? "true" : "false");
    printField("single_pe_model", detector_cfg.single_pe.model);
    printField("single_pe_unit", detector_cfg.single_pe.unit);
    printField("template_time_reference",
               detector_cfg.single_pe.template_time_reference);
    printField("charge_fluctuation_enabled",
               detector_cfg.single_pe.charge_fluctuation.enabled ? "true"
                                                                 : "false");
    printField("charge_fluctuation_model",
               detector_cfg.single_pe.charge_fluctuation.model);
    printField("charge_sample_count",
               intToString(detector_cfg.single_pe.charge_fluctuation
                               .empirical_samples.size()));
    printField("time_jitter_enabled",
               detector_cfg.single_pe.time_jitter.enabled ? "true" : "false");
    printField("time_jitter_sigma_ns",
               doubleToString(detector_cfg.single_pe.time_jitter.sigma_ns));
    printField("sampling_width_ns",
               doubleToString(detector_cfg.sampling.width_ns));
    printField("save_primary_sequence",
               detector_cfg.save_primary_sequence ? "true" : "false");
    printField("save_fired_sequence",
               detector_cfg.save_fired_sequence ? "true" : "false");

    printSection("Waveform");
    printField("enabled", waveform_cfg.enabled ? "true" : "false");
    printField("source", waveform_cfg.source);
    printField("time_reference", waveform_cfg.time_reference);
    printField("time_bin_width_ns", doubleToString(waveform_cfg.time_bin_width_ns));
    printField("time_window_start_ns", doubleToString(waveform_cfg.time_window_start_ns));
    printField("time_window_end_ns", doubleToString(waveform_cfg.time_window_end_ns));
    printField("model",
               waveform_cfg.source == "electronics"
                   ? "single-p.e. superposition and configured sampling"
                   : "proxy time-binned camera output");

    printSection("Profile");
    printField("enabled", getBool(cfg, "profile.enabled", false) ? "true" : "false");
    printField("atmosphere_height_histogram",
               getString(cfg, "atmosphere.height_histogram_csv",
                         getString(cfg, "atmosphere.histogram_csv", "off")));

    printSection("NSB");
    printField("enabled", nsb_cfg.enabled ? "true" : "false");
    printField("model", nsb_cfg.model);
    printField("rate_pe_per_ns_per_pixel",
               doubleToString(nsb_cfg.rate_pe_per_ns_per_pixel));
    printField("window_ns", doubleToString(nsb_cfg.window_ns));
    printField("seed", intToString(nsb_cfg.seed));
    if (nsb_cfg.model == "spectral_flux") {
        printField("spectrum_csv", nsb_cfg.spectrum_csv);
        printField("spectrum_unit", nsb_cfg.spectrum_unit);
        printField("effective_area_m2", doubleToString(nsb_cfg.effective_area_m2));
        printField("collector_mean_transmission",
                   doubleToString(nsb_cfg.collector_mean_transmission));
        printField("microcell_geometric_acceptance",
                   doubleToString(nsb_cfg.microcell_geometric_acceptance));
        printField("pixel_solid_angle_sr", doubleToString(nsb_cfg.pixel_solid_angle_sr));
        printField("computed_from_spectrum",
                   nsb_cfg.computed_from_spectrum ? "true" : "false");
        printField("spectral_integral_pe_s_sr_m2",
                   doubleToString(nsb_cfg.spectral_integral_pe_s_sr_m2));
    }

    printSection("Trigger");
    printField("enabled", trigger_cfg.enabled ? "true" : "false");
    printField("mode",
               detector_cfg.enabled
                   ? detector_cfg.camera_trigger.mode
                   : "pe_count");
    if (detector_cfg.enabled &&
        detector_cfg.camera_trigger.mode == "voltage") {
        printField(
            "pixel_threshold_mv",
            doubleToString(
                detector_cfg.camera_trigger.pixel_threshold_mv));
    } else {
        printField("pixel_threshold_pe",
                   doubleToString(trigger_cfg.pixel_threshold_pe));
    }
    printField("camera_multiplicity", intToString(trigger_cfg.camera_multiplicity));
    printField("array_enabled",
               trigger_cfg.array_enabled ? "true" : "false");
    printField("array_multiplicity", intToString(trigger_cfg.array_multiplicity));
    printField("camera_coincidence_window_ns",
               doubleToString(trigger_cfg.camera_coincidence_window_ns));
    printField("array_coincidence_window_ns",
               doubleToString(trigger_cfg.array_coincidence_window_ns));
    printField("array_time_correction", trigger_cfg.array_time_correction);
    printField("array_wavefront_speed_m/ns",
               trigger_cfg.array_wavefront_speed_m_per_ns > 0.0
                   ? doubleToString(
                         trigger_cfg.array_wavefront_speed_m_per_ns, 9)
                   : "auto (EventIO observation altitude)");

    printSection("Photon response");
    printField("mode", response_cfg.modeName());
    printField("seed", intToString(response_cfg.seed));

    printSection("Efficiency");
    printField("constant_scale", doubleToString(efficiency_cfg.constant_scale));
    printField("mirror_reflectivity",
               factorDescription(efficiency_cfg.mirror_reflectivity));
    const std::string facet_scale_csv =
        getString(cfg, "efficiency.mirror_reflectivity_scale_csv", "");
    printField("mirror_scale",
               isDisabledText(facet_scale_csv)
                   ? "uniform fallback: " + getString(
                         cfg, "efficiency.mirror_reflectivity_scale", "1")
                   : "per-facet CSV: " + facet_scale_csv);
    printField("filter_transmission",
               factorDescription(efficiency_cfg.filter_transmission));
    printField("atmosphere", factorDescription(efficiency_cfg.atmosphere_transmission));
    printField("atmosphere_model", atmosphereTransmissionDescription(atmosphere_cfg));
    printField("funnel_acceptance",
               efficiency_cfg.use_funnel_acceptance ? "cos(theta)" : "not set -> 1");

    printSection("Errors");
    printField("random_seed", intToString(error_cfg.random_seed));
    printField("structural_deformation",
               isDisabledText(error_cfg.structural_deformation_config) ? "off" : "on");
    if (!isDisabledText(error_cfg.structural_deformation_config)) {
        printField("structural_deformation_config",
                   error_cfg.structural_deformation_config);
        printField("structural_deformation_elevation_deg",
                   doubleToString(telescope_cfg.pointing_el_deg));
    }
    printField("facet_radial_pos_sigma_m",
               doubleToString(error_cfg.facet_radial_position_sigma_m));
    printField("facet_normal_sigma_deg", doubleToString(error_cfg.facet_normal_sigma_deg));
    printField("reflect_dir_sigma_deg",
               doubleToString(error_cfg.reflect_direction_sigma_deg));
    printField("radius_curvature_sigma_m",
               doubleToString(error_cfg.radius_of_curvature_sigma_m));
    printField("reflectivity_scale_sigma",
               doubleToString(error_cfg.reflectivity_scale_sigma));

    printSection("Obstruction");
    printField("enabled", obstruction.enabled ? "true" : "false");
    if (obstruction.enabled) {
        printField("mode", obstruction.mode);
        printField("check_incoming", obstruction.check_incoming ? "true" : "false");
        printField("check_reflected", obstruction.check_reflected ? "true" : "false");
        if (obstruction.mode == "primitives") {
            printField("primitives_csv", obstruction.primitives_csv);
            printField("primitive_count",
                       intToString(static_cast<std::uint64_t>(obstruction.primitives.size())));
        } else {
            printField("mask_csv", obstruction.mask_csv);
            printField("plane_z_m", doubleToString(obstruction.plane_z_m));
            printField("grid",
                       intToString(static_cast<std::uint64_t>(obstruction.nx)) +
                       " x " +
                       intToString(static_cast<std::uint64_t>(obstruction.ny)));
            printField("cell_size_m", doubleToString(obstruction.cell_size_m));
        }
    }

    printSection("Model");
    printField("optics", "facet reflection with configured optical errors");
    printField("speed_of_light_m/ns",
               doubleToString(propagation_cfg.speed_of_light_m_per_ns, 9));
    std::string missing = "crosstalk, afterpulse, dark count";
    if (!detector_cfg.enabled) {
        missing = "explicit microcell saturation, single-p.e. waveform, " +
                  missing;
    }
    if (!light_collector) {
        missing = "collector, " + missing;
    }
    if (!waveform_cfg.enabled) {
        missing = "waveform output, " + missing;
    }
    printField("not included", missing);
}

void printProfileStats(const ProfileStats& stats, double trace_time_s)
{
    printSection("Profile");
    printField("note",
               "eventio_stream_s is wall time for streaming plus callback processing; "
               "sub-stage times are measured inside that callback and are not exclusive");
    if (trace_time_s <= 0.0) {
        trace_time_s = 1.0;
    }
    auto print_item = [trace_time_s](const std::string& name, double seconds) {
        std::ostringstream value;
        value << doubleToString(seconds, 6)
              << " (" << doubleToString(100.0 * seconds / trace_time_s, 2) << "%)";
        printField(name, value.str());
    };
    print_item("eventio_stream_s", stats.eventio_stream_s);
    print_item("transform_s", stats.transform_s);
    print_item("trace_to_plane_s", stats.trace_to_plane_s);
    print_item("obstruction_s", stats.obstruction_s);
    print_item("camera_response_s", stats.camera_response_s);
    print_item("whiteboard_accumulate_s", stats.whiteboard_accumulate_s);
    print_item("camera_accumulate_s", stats.camera_accumulate_s);
    print_item("hdf5_write_s", stats.hdf5_write_s);
}

double poissonUpperTail(int threshold, double lambda)
{
    if (threshold <= 0) {
        return 1.0;
    }
    if (lambda <= 0.0) {
        return 0.0;
    }
    double term = std::exp(-lambda);
    double cdf = term;
    for (int k = 1; k < threshold; ++k) {
        term *= lambda / static_cast<double>(k);
        cdf += term;
    }
    return std::clamp(1.0 - cdf, 0.0, 1.0);
}

double binomialUpperTail(int threshold, std::size_t trials, double p)
{
    if (threshold <= 0) {
        return 1.0;
    }
    if (trials == 0 || p <= 0.0 || static_cast<std::size_t>(threshold) > trials) {
        return 0.0;
    }
    if (p >= 1.0) {
        return 1.0;
    }
    long double tail = 0.0L;
    for (std::size_t k = static_cast<std::size_t>(threshold); k <= trials; ++k) {
        const long double log_term =
            std::lgamma(static_cast<long double>(trials) + 1.0L) -
            std::lgamma(static_cast<long double>(k) + 1.0L) -
            std::lgamma(static_cast<long double>(trials - k) + 1.0L) +
            static_cast<long double>(k) * std::log(static_cast<long double>(p)) +
            static_cast<long double>(trials - k) *
                std::log1p(-static_cast<long double>(p));
        tail += std::exp(log_term);
        if (tail >= 1.0L) {
            return 1.0;
        }
    }
    return static_cast<double>(std::min<long double>(1.0L, tail));
}

std::string scientificDoubleToString(double value, int precision = 6)
{
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(precision) << value;
    return oss.str();
}

void printPureNsbTriggerEstimate(const NsbConfig& nsb_cfg,
                                 const WaveformOutputConfig& waveform_cfg,
                                 const TriggerConfig& trigger_cfg,
                                 std::size_t n_pixels,
                                 std::size_t n_telescopes)
{
    printSection("Pure NSB trigger estimate");
    if (!nsb_cfg.enabled || nsb_cfg.rate_pe_per_ns_per_pixel <= 0.0 ||
        !trigger_cfg.enabled || n_pixels == 0) {
        printField("enabled", "false");
        printField("reason", "NSB, trigger, or camera pixels disabled");
        return;
    }

    const double configured_window_ns =
        trigger_cfg.camera_coincidence_window_ns > 0.0
            ? trigger_cfg.camera_coincidence_window_ns
            : std::max(0.0,
                       waveform_cfg.time_window_end_ns -
                           waveform_cfg.time_window_start_ns);
    const std::size_t window_bins =
        waveform_cfg.enabled && waveform_cfg.time_bin_width_ns > 0.0
            ? std::max<std::size_t>(
                  1,
                  static_cast<std::size_t>(
                      std::ceil(configured_window_ns /
                                waveform_cfg.time_bin_width_ns)))
            : 1;
    const double effective_window_ns =
        waveform_cfg.enabled && waveform_cfg.time_bin_width_ns > 0.0
            ? static_cast<double>(window_bins) * waveform_cfg.time_bin_width_ns
            : configured_window_ns;
    const double full_window_ns =
        waveform_cfg.enabled
            ? std::max(0.0,
                       waveform_cfg.time_window_end_ns -
                           waveform_cfg.time_window_start_ns)
            : effective_window_ns;
    const std::size_t n_windows =
        waveform_cfg.enabled && waveform_cfg.time_bin_width_ns > 0.0
            ? static_cast<std::size_t>(std::ceil(
                  full_window_ns / waveform_cfg.time_bin_width_ns))
            : 1;
    const double lambda =
        nsb_cfg.rate_pe_per_ns_per_pixel * effective_window_ns;
    const double pixel_prob =
        poissonUpperTail(static_cast<int>(std::ceil(trigger_cfg.pixel_threshold_pe)),
                         lambda);
    const double camera_single_window_prob =
        binomialUpperTail(trigger_cfg.camera_multiplicity, n_pixels, pixel_prob);
    const double camera_sliding_upper =
        std::min(1.0, static_cast<double>(n_windows) * camera_single_window_prob);
    const double array_upper =
        n_telescopes > 0
            ? binomialUpperTail(trigger_cfg.array_multiplicity,
                                n_telescopes,
                                camera_sliding_upper)
            : 0.0;

    printField("enabled", "true");
    printField("rate_pe_per_ns_per_pixel",
               doubleToString(nsb_cfg.rate_pe_per_ns_per_pixel));
    printField("effective_window_ns", doubleToString(effective_window_ns));
    printField("pixel_threshold_pe", doubleToString(trigger_cfg.pixel_threshold_pe));
    printField("camera_multiplicity", intToString(trigger_cfg.camera_multiplicity));
    printField("array_multiplicity", intToString(trigger_cfg.array_multiplicity));
    printField("n_pixels", intToString(static_cast<std::uint64_t>(n_pixels)));
    printField("n_telescopes", intToString(static_cast<std::uint64_t>(n_telescopes)));
    printField("pixel_prob_ge_threshold", scientificDoubleToString(pixel_prob, 6));
    printField("camera_single_window_prob",
               scientificDoubleToString(camera_single_window_prob, 6));
    printField("camera_sliding_upper_prob",
               scientificDoubleToString(camera_sliding_upper, 6));
    printField("array_sliding_upper_prob", scientificDoubleToString(array_upper, 6));
}

bool shouldHideInputCardLine(const std::string& line)
{
    const std::string text = lowerCopy(trim(line));
    return startsWith(text, "telfil") || startsWith(text, "direct");
}

void printEventIOMetadataSummary(const EventIOMetadata& metadata)
{
    printSection("EventIO metadata");
    printField("shower_events", intToString(metadata.events.size()));
    printField("telescopes", intToString(metadata.telescopes.size()));
    if (!metadata.events.empty()) {
        int min_event = metadata.events.front().shower_event_id;
        int max_event = metadata.events.front().shower_event_id;
        for (const auto& event : metadata.events) {
            min_event = std::min(min_event, event.shower_event_id);
            max_event = std::max(max_event, event.shower_event_id);
        }
        printField("shower_event_range",
                   intToString(min_event) + ".." + intToString(max_event));

        const auto& first = metadata.events.front();
        printField("first_event_theta_deg", doubleToString(first.theta_deg));
        printField("first_event_phi_deg", doubleToString(first.phi_deg));
        printField("first_event_arrang_deg", doubleToString(first.array_rotation_deg));
        printField("first_event_az_N_to_E_deg",
                   doubleToString(first.azimuth_north_to_east_deg));
    }

    std::size_t hidden_input_lines = 0;
    printField("input_card_lines", intToString(metadata.input_lines.size()));
    for (std::size_t i = 0; i < metadata.input_lines.size(); ++i) {
        if (shouldHideInputCardLine(metadata.input_lines[i])) {
            ++hidden_input_lines;
            continue;
        }
        printField("input[" + intToString(i) + "]", metadata.input_lines[i]);
    }
    if (hidden_input_lines > 0) {
        printField("input_hidden_path_lines", intToString(hidden_input_lines));
    }
}

} // namespace lact
