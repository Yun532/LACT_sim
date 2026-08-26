#include "app/OpticalSimCommon.hpp"
#include "app/PhotonTransport.hpp"
#include "io/CorsikaTraceEventMetadata.hpp"
#include "io/CorsikaTraceHdf5Writer.hpp"
#include "io/CorsikaTraceStats.hpp"
#include "app/CorsikaTraceReport.hpp"
#include "app/CorsikaTraceConfig.hpp"
#include "io/CorsikaTraceAccumulators.hpp"
#include "io/CorsikaTraceCsvOutput.hpp"
#include "io/CorsikaTraceEventIOInput.hpp"
#include "app/PhotonResponseSampler.hpp"
#include "io/EventIOArrayTiming.hpp"
#include "app/TelescopeOpticsCache.hpp"
#include "app/TriggerResponse.hpp"
#include "core/Sha256.hpp"
#include "io/CameraElectronicsEvent.hpp"
#include "io/CorsikaTraceOutputTypes.hpp"

#ifdef LACT_HAS_HDF5
#include <hdf5.h>
#include "io/Hdf5WaveformWriter.hpp"
#endif

#ifdef LACT_HAS_ROOT
#include "io/LactEventRootWriter.hpp"
#endif

#include <algorithm>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <tuple>

using namespace lact;

namespace {


} // namespace

int main(int argc, char** argv) {
    try {
        const auto command = parseConfigCommandLine(argc, argv);
        if (command.help || command.positional.empty() ||
            command.positional.size() > 2) {
            std::cerr
                << "usage: run_corsika_trace <config.txt> "
                   "[corsika_eventio_file] [-C key=value ...]\n";
            return command.help ? 0 : 2;
        }
        const std::string config_path = command.positional[0];
#ifndef LACT_HAS_HESSIO
        throw std::runtime_error(
            "run_corsika_trace requires libhessio. Build external/hessioxxx/source "
            "and reconfigure LACT_sim.");
#else
        std::cout.setf(std::ios::unitbuf);
        std::cout << std::fixed << std::setprecision(6);
        const auto t_start = std::chrono::steady_clock::now();
        auto main_cfg = readKeyValueConfig(config_path);
        applyConfigOverrides(main_cfg, command.overrides);
        ComponentConfigPaths component_paths;
        auto cfg = expandConfig(main_cfg, config_path, component_paths);
        applyConfigOverrides(cfg, command.overrides);

        cfg["source.mode"] = getString(cfg, "source.mode", "EventIO");
        const bool photon_csv_mode =
            isPhotonCsvMode(getString(cfg, "source.mode", "EventIO"));
        if (!isEventIOMode(getString(cfg, "source.mode", "EventIO")) &&
            !photon_csv_mode) {
            throw std::runtime_error(
                "run_corsika_trace requires source.mode=EventIO/CORSIKA or PhotonCsv");
        }

        SyntheticPhotonConfig source_cfg = buildSourceConfig(cfg);
        SourceRuntimeConfig source_runtime_cfg = buildSourceRuntimeConfig(cfg);
        if (cfg.find("source.coordinate_frame") == cfg.end() &&
            cfg.find("source.eventio_coordinate_frame") != cfg.end()) {
            std::cerr << "warning: source.eventio_coordinate_frame is deprecated; "
                         "use source.coordinate_frame instead\n";
        }
        if (command.positional.size() == 2 && !photon_csv_mode) {
            source_runtime_cfg.eventio_path = command.positional[1];
            cfg["source.eventio_path"] = command.positional[1];
        }
        if (photon_csv_mode && source_runtime_cfg.csv_path.empty()) {
            throw std::runtime_error(
                "source.csv_path is required when source.mode=PhotonCsv");
        }
        if (source_runtime_cfg.eventio_path.empty()) {
            source_runtime_cfg.eventio_path =
                getString(cfg, "source.metadata_eventio_path",
                          getString(cfg, "corsika.input",
                                    getString(cfg, "eventio.input", "")));
        }
        if (!photon_csv_mode && source_runtime_cfg.eventio_path.empty()) {
            throw std::runtime_error(
                "source.eventio_path or corsika.input is required");
        }
        if (!source_runtime_cfg.eventio_path.empty()) {
            cfg["source.eventio_path"] = source_runtime_cfg.eventio_path;
        }
        TelescopeConfig telescope_cfg = buildTelescopeConfig(cfg);
        std::vector<MirrorFacet> nominal_facets = buildFacetsFromConfig(cfg);
        ErrorConfig error_cfg = buildErrorConfig(cfg);
        ObstructionMask obstruction = buildObstructionMask(cfg);
        applyStructuralDeformation(nominal_facets, error_cfg, telescope_cfg);
        applyFacetEfficiencyScales(nominal_facets, cfg);
        OutputPlane plane = buildOutputPlane(cfg);
        TelescopeFrame telescope_frame;
        const std::string trace_frame_name = normalizeSourceCoordinateFrame(
            source_runtime_cfg.coordinate_frame);
        if (trace_frame_name == "corsika_nwu_relative" ||
            trace_frame_name == "corsika_nwu_global") {
            telescope_frame = buildCorsikaNwuTelescopeFrame(telescope_cfg);
        } else if (trace_frame_name == "enu_east_relative" ||
                   trace_frame_name == "enu_east_global") {
            telescope_frame = buildEnuEastTelescopeFrame(telescope_cfg);
        } else if (trace_frame_name == "lact_generic_global") {
            telescope_frame = buildTelescopeFrame(telescope_cfg);
        }
        CameraConfig camera_cfg = buildCameraConfig(cfg);
        SipmConfig sipm_cfg = buildSipmConfig(cfg);
        ElectronicsConfig electronics_cfg = buildElectronicsConfig(cfg);
        ElectronicsResponse electronics(electronics_cfg);
        const auto detector_pipeline_cfg = buildDetectorPipelineConfig(cfg);
        validateCameraDetectorCompatibility(camera_cfg, detector_pipeline_cfg);
        NsbConfig nsb_cfg = buildNsbConfig(cfg);
        TriggerConfig trigger_cfg = buildTriggerConfig(cfg);
        CameraGeometry camera = buildCameraGeometry(camera_cfg);
        std::vector<int> camera_pixel_id_axis;
        camera_pixel_id_axis.reserve(camera.size());
        for (const auto& pixel : camera.pixels()) {
            camera_pixel_id_axis.push_back(pixel.id);
        }
        std::sort(camera_pixel_id_axis.begin(), camera_pixel_id_axis.end());
        auto light_collector = buildLightCollector(camera_cfg, camera);
        TelescopeOpticsCache telescope_optics(nominal_facets, error_cfg);
        const MirrorLayout& mirrors = telescope_optics.layoutFor(telescope_cfg.id);
        OpticalEfficiencyConfig efficiency_cfg = buildEfficiencyConfig(cfg);
        resolveNsbSpectralRate(nsb_cfg, efficiency_cfg, camera, telescope_cfg,
                               &detector_pipeline_cfg);
        AtmosphereTransmissionConfig atmosphere_cfg = buildAtmosphereTransmissionConfig(cfg);
        PropagationConfig propagation_cfg = buildPropagationConfig(cfg);
        OpticalEfficiency eff(efficiency_cfg);
        OpticalTracer tracer(propagation_cfg.speed_of_light_m_per_ns,
                             effectiveReflectDirectionSigmaRad(nominal_facets, error_cfg),
                             error_cfg.random_seed);
        const double eventio_mirror_front_z_m = mirrorFrontReferenceZ(mirrors);
        const bool eventio_2d_backproject =
            shouldBackprojectEventIO2d(source_runtime_cfg);
        CorsikaTraceOutputConfig output_cfg = buildCorsikaTraceOutputConfig(cfg);
        WaveformOutputConfig waveform_cfg = buildWaveformOutputConfig(cfg);
        if (waveform_cfg.enabled &&
            waveform_cfg.source == "electronics" &&
            !detector_pipeline_cfg.enabled) {
            throw std::runtime_error(
                "waveform.source=electronics requires "
                "electronics.enabled=true");
        }
        CollectorDebugConfig collector_debug_cfg = buildCollectorDebugConfig(cfg);
        ProfileConfig profile_cfg = buildProfileConfig(cfg);
        AtmosphereHistogramConfig atmosphere_histogram_cfg =
            buildAtmosphereHistogramConfig(cfg);
        std::vector<AtmosphereHistogramBin> atmosphere_histogram =
            makeAtmosphereHistogramBins(atmosphere_histogram_cfg);
        auto eventio_cfg = buildEventIOPhotonConfig(cfg, source_cfg, source_runtime_cfg);
        const PhotonResponseConfig response_cfg = buildPhotonResponseConfig(cfg);
        ProfileStats profile_stats;
        const bool save_csv = outputWantsCsv(output_cfg);
        const bool save_hdf5 = outputWantsHdf5(output_cfg);
        const bool save_lact_root = outputWantsLactRoot(output_cfg);
        if (save_hdf5) {
#ifndef LACT_HAS_HDF5
            throw std::runtime_error(
                "output.format requests HDF5, but this build was configured without HDF5. "
                "Install HDF5 and re-run CMake, or set output.format=csv.");
#endif
        }
        if (save_lact_root) {
#ifndef LACT_HAS_ROOT
            throw std::runtime_error(
                "output.lact_root_enabled=true, but this build was configured without ROOT. "
                "Load/install ROOT >= 6.24 and re-run CMake, or disable output.lact_root_enabled.");
#endif
        }
        cfg["provenance.producer_version"] = LACT_PRODUCER_VERSION;
        if (save_hdf5 || save_lact_root) {
            const std::string provenance_source_path =
                source_runtime_cfg.use_photon_csv
                    ? source_runtime_cfg.csv_path
                    : source_runtime_cfg.eventio_path;
            cfg["provenance.source_path"] = provenance_source_path;
            std::cerr << "run_corsika_trace: hashing input source "
                      << provenance_source_path << "\n";
            cfg["provenance.source_sha256"] =
                sha256File(provenance_source_path);
        }

        std::cout << "========================================\n";
        std::cout << "LACT CORSIKA/EventIO trace\n";
        std::cout << "========================================\n";
        printSection("Configuration files");
        printField("producer_version",
                   getString(cfg, "provenance.producer_version", "source-tree"));
        printField("main", config_path);
        for (const auto& [key, value] : command.overrides) {
            printField("override", key + "=" + value);
        }
        if (!component_paths.telescope.empty()) printField("telescope", component_paths.telescope);
        if (!component_paths.mirror.empty()) printField("mirror", component_paths.mirror);
        if (!component_paths.source.empty()) printField("source", component_paths.source);
        if (!component_paths.output.empty()) printField("output", component_paths.output);
        if (!component_paths.camera.empty()) printField("camera", component_paths.camera);
        if (!component_paths.sipm.empty()) printField("sipm", component_paths.sipm);
        if (!component_paths.electronics.empty()) printField("electronics", component_paths.electronics);
        if (!component_paths.efficiency.empty()) printField("efficiency", component_paths.efficiency);
        if (!component_paths.atmosphere.empty()) printField("atmosphere", component_paths.atmosphere);
        if (!component_paths.nsb.empty()) printField("nsb", component_paths.nsb);
        if (!component_paths.trigger.empty()) printField("trigger", component_paths.trigger);
        if (!component_paths.error.empty()) printField("error", component_paths.error);
        if (!component_paths.obstruction.empty()) {
            printField("obstruction", component_paths.obstruction);
        }
        if (component_paths.source.empty()) {
            printField("source", photon_csv_mode
                                     ? "inline PhotonCsv settings"
                                     : "inline EventIO settings");
        }
        if (cfg.find("provenance.source_sha256") != cfg.end()) {
            printField("source_sha256", cfg.at("provenance.source_sha256"));
        }

        printSection("Run");
        EventIOMetadata metadata;
        if (!eventio_cfg.path.empty()) {
            printField("status", "reading EventIO metadata");
            std::cerr << "run_corsika_trace: reading EventIO metadata from "
                      << eventio_cfg.path << "\n";
            metadata = readEventIOMetadata(eventio_cfg);
        } else {
            printField("status", "using PhotonCsv and configuration metadata");
            std::cerr << "run_corsika_trace: PhotonCsv has no EventIO metadata; "
                         "using telescope/source configuration defaults\n";
        }
        const bool explicit_wavelength_range = hasExplicitMissingWavelengthRange(cfg);
        const bool has_cwavlg = wavelengthRangeFromInputCard(metadata).has_value();
        applyEventIOWavelengthMetadata(eventio_cfg, metadata, cfg);
        applyEventIOAtmosphereMetadata(eventio_cfg, metadata);
        if (cfg.find("atmosphere.detector_altitude_km") == cfg.end() &&
            cfg.find("atmosphere.ground_altitude_km") == cfg.end() &&
            std::isfinite(metadata.observation_altitude_m)) {
            atmosphere_cfg.detector_altitude_km =
                metadata.observation_altitude_m * 1.0e-3;
        }
        AtmosphereTransmission atmosphere(atmosphere_cfg);
        atmosphere_cfg = atmosphere.config();
        PhotonResponseSampler response_sampler(response_cfg, eventio_cfg);
        const std::string missing_wavelength_range_source =
            explicit_wavelength_range ? "cfg"
                                      : (has_cwavlg ? "EventIO input card CWAVLG"
                                                   : "built-in default");

        printCorsikaOpticalConfiguration(cfg,
                                         telescope_cfg,
                                         telescope_frame,
                                         mirrors,
                                         source_cfg,
                                         source_runtime_cfg,
                                         plane,
                                         output_cfg,
                                         camera_cfg,
                                         camera,
                                         light_collector,
                                         sipm_cfg,
                                         electronics_cfg,
                                         detector_pipeline_cfg,
                                         waveform_cfg,
                                         collector_debug_cfg,
                                         nsb_cfg,
                                         trigger_cfg,
                                         response_cfg,
                                         efficiency_cfg,
                                         atmosphere_cfg,
                                         error_cfg,
                                         obstruction,
                                         propagation_cfg,
                                         eventio_cfg,
                                         missing_wavelength_range_source,
                                         eventio_mirror_front_z_m,
                                         eventio_2d_backproject);

        if (!eventio_cfg.path.empty()) {
            printEventIOMetadataSummary(metadata);
        }
        printPureNsbTriggerEstimate(nsb_cfg,
                                    waveform_cfg,
                                    trigger_cfg,
                                    camera.size(),
                                    metadata.telescopes.empty()
                                        ? 1
                                        : metadata.telescopes.size());
	    printField("status", photon_csv_mode
                                  ? "loading PhotonCsv bunches and tracing"
                                  : "streaming EventIO photon bunches and tracing");
        std::cerr << "run_corsika_trace: "
                  << (photon_csv_mode
                          ? "loading PhotonCsv bunches\n"
                          : "streaming photon bunches; no full-file photon preload is used\n");

        std::ofstream whiteboard_out;
        if (!camera_cfg.enabled && save_csv) {
            const std::filesystem::path out_path(output_cfg.hits_csv);
            if (out_path.has_parent_path()) {
                std::filesystem::create_directories(out_path.parent_path());
            }
            whiteboard_out.open(output_cfg.hits_csv);
            if (!whiteboard_out) {
                throw std::runtime_error("failed to write whiteboard CSV: " +
                                         output_cfg.hits_csv);
            }
            writeCorsikaWhiteboardHeader(whiteboard_out,
                                         output_cfg.whiteboard_emitter_info);
        }

        std::ofstream mirror_diagnostic_out;
        if (!output_cfg.mirror_diagnostic_csv.empty()) {
            const std::filesystem::path out_path(
                output_cfg.mirror_diagnostic_csv);
            if (out_path.has_parent_path()) {
                std::filesystem::create_directories(out_path.parent_path());
            }
            mirror_diagnostic_out.open(output_cfg.mirror_diagnostic_csv);
            if (!mirror_diagnostic_out) {
                throw std::runtime_error(
                    "failed to write mirror diagnostic CSV: " +
                    output_cfg.mirror_diagnostic_csv);
            }
            writeCorsikaMirrorDiagnosticHeader(mirror_diagnostic_out);
        }

        std::map<SummaryKey, TraceSummary> summaries;
        std::map<PixelKey, PixelAccumulator> pixels;
        std::map<WaveformKey, WaveformPixelAccumulator> waveforms;
        std::vector<RawWaveformHit> raw_waveform_hits;
        CameraElectronicsEventMap electronics_events;
        std::vector<CollectorDebugPhotonRow> collector_debug_rows;
        std::vector<WhiteboardHdf5Row> whiteboard_hits;

#ifdef LACT_HAS_ROOT
        std::unique_ptr<LactEventRootStreamWriter> lact_root_stream_writer;
        if (save_lact_root) {
            printSection("lact_event ROOT output");
            printField("status", "opening streaming lact_event ROOT writer");
            printField("path", output_cfg.lact_root_path);
            printField("profile", output_cfg.lact_profile);
            lact_root_stream_writer = std::make_unique<LactEventRootStreamWriter>(
                output_cfg,
                waveform_cfg,
                config_path,
                cfg,
                source_runtime_cfg,
                telescope_cfg,
                metadata,
                camera,
                nominal_facets,
                nsb_cfg,
                trigger_cfg);
        }
#endif

        std::uint64_t photon_index = 0;
        int active_shower_event = -1;
        int active_output_event = -1;
        const auto t_trace_start = std::chrono::steady_clock::now();

        auto prepareOutputEvent = [&](int event_id) {
            if (event_id < 0) {
                return;
            }
            std::map<SummaryKey, TraceSummary> event_summaries;
            for (const auto& kv : summaries) {
                if (kv.first.first == event_id) {
                    event_summaries.insert(kv);
                }
            }
            std::map<PixelKey, PixelAccumulator> event_pixels;
            for (auto it = pixels.lower_bound(PixelKey{event_id, std::numeric_limits<int>::min(),
                                                       std::numeric_limits<int>::min()});
                 it != pixels.end() && std::get<0>(it->first) == event_id;
                 ++it) {
                event_pixels.insert(*it);
            }
            std::map<WaveformKey, WaveformPixelAccumulator> event_waveforms;
            for (auto it = waveforms.lower_bound(WaveformKey{
                     event_id, std::numeric_limits<int>::min(),
                     std::numeric_limits<int>::min(), std::numeric_limits<int>::min()});
                 it != waveforms.end() && std::get<0>(it->first) == event_id;
                 ++it) {
                event_waveforms.insert(*it);
            }
            std::vector<RawWaveformHit> event_raw_waveform_hits;
            for (const auto& hit : raw_waveform_hits) {
                if (hit.event_id == event_id) {
                    event_raw_waveform_hits.push_back(hit);
                }
            }

            auto event_electronics = buildCameraElectronicsEvents(
                detector_pipeline_cfg,
                waveform_cfg,
                nsb_cfg,
                source_runtime_cfg,
                camera_pixel_id_axis,
                event_summaries,
                event_raw_waveform_hits);
            if (detector_pipeline_cfg.enabled) {
                std::vector<TelescopeTriggerTime> telescope_times;
                std::vector<int> all_telescope_ids;
                all_telescope_ids.reserve(event_electronics.size());
                for (const auto& item : event_electronics) {
                    all_telescope_ids.push_back(item.second.telescope_id);
                    if (item.second.detector.camera_trigger.triggered) {
                        telescope_times.push_back(TelescopeTriggerTime{
                            item.second.telescope_id,
                            item.second.trigger_time_ns,
                        });
                    }
                }
                const auto geometric_delays = eventIOArrayGeometricDelaysNs(
                    all_telescope_ids,
                    event_id,
                    source_runtime_cfg.event_id_mode,
                    trigger_cfg,
                    telescope_cfg,
                    metadata);
                applyEventIOArrayTimingCorrection(
                    telescope_times,
                    event_id,
                    source_runtime_cfg.event_id_mode,
                    trigger_cfg,
                    telescope_cfg,
                    metadata);
                const auto array_decision =
                    evaluateArrayTrigger(telescope_times, trigger_cfg);
                for (auto& item : event_electronics) {
                    auto& event = item.second;
                    event.geometric_delay_ns =
                        geometric_delays.at(event.telescope_id);
                    const auto corrected = std::find_if(
                        telescope_times.begin(), telescope_times.end(),
                        [&event](const TelescopeTriggerTime& value) {
                            return value.telescope_id == event.telescope_id;
                        });
                    if (corrected != telescope_times.end()) {
                        event.geometric_delay_ns = corrected->geometric_delay_ns;
                        event.coincidence_time_ns =
                            std::isfinite(corrected->coincidence_time_ns)
                                ? corrected->coincidence_time_ns
                                : corrected->trigger_time_ns;
                    }
                    event.array_triggered = array_decision.triggered;
                    const bool telescope_is_coincident = std::binary_search(
                        array_decision.coincident_telescope_ids.begin(),
                        array_decision.coincident_telescope_ids.end(),
                        event.telescope_id);
                    if (!output_cfg.save_only_triggered) {
                        event.selected_for_output =
                            hasFinalIntegratedImageSignal(event.detector);
                    } else if (!trigger_cfg.enabled) {
                        event.selected_for_output = true;
                    } else {
                        event.selected_for_output =
                            event.detector.camera_trigger.triggered &&
                            array_decision.triggered && telescope_is_coincident;
                    }
                }
            }
#ifdef LACT_HAS_ROOT
            if (lact_root_stream_writer) {
            lact_root_stream_writer->writeEvent(event_summaries,
                                                event_pixels,
                                                event_waveforms,
                                                event_raw_waveform_hits,
                                                event_electronics);
            }
#endif
            if (save_hdf5 || save_csv) {
                for (auto& item : event_electronics) {
                    electronics_events.emplace(item.first,
                                               std::move(item.second));
                }
            }
        };

        const bool can_release_streamed_event_data =
            save_lact_root && !save_hdf5 && !save_csv;
        auto releaseStreamedEventData = [&](int event_id) {
            if (!can_release_streamed_event_data || event_id < 0) {
                return;
            }
            for (auto it = pixels.begin(); it != pixels.end();) {
                if (std::get<0>(it->first) == event_id) {
                    it = pixels.erase(it);
                } else {
                    ++it;
                }
            }
            for (auto it = waveforms.begin(); it != waveforms.end();) {
                if (std::get<0>(it->first) == event_id) {
                    it = waveforms.erase(it);
                } else {
                    ++it;
                }
            }
            raw_waveform_hits.erase(
                std::remove_if(raw_waveform_hits.begin(), raw_waveform_hits.end(),
                               [event_id](const RawWaveformHit& hit) {
                                   return hit.event_id == event_id;
                               }),
                raw_waveform_hits.end());
        };

        auto flushOutputEvent = [&](int event_id) {
            prepareOutputEvent(event_id);
            releaseStreamedEventData(event_id);
        };

        auto printRunningEventSummary = [&](int shower_event) {
            if (shower_event < 0) {
                return;
            }
            printSection("Event shower_event=" + intToString(shower_event));
            printField("event_status", "completed");
            printField("telescope_columns",
                       camera_cfg.enabled
                           ? "telescope output_events input_bunches input_photons hit_mirror hit_output hit_camera accepted lost unique_pixels pe time_mean_ns time_rms_ns time_first_ns"
                           : "telescope output_events input_bunches input_photons hit_mirror hit_output signal time_mean_ns time_rms_ns time_first_ns");
            std::uint64_t input_bunches = 0;
            double input_photons = 0.0;
            std::uint64_t blocked = 0;
            std::uint64_t blocked_incoming = 0;
            std::uint64_t blocked_reflected = 0;
            std::uint64_t hit_output_before_obstruction = 0;
            std::uint64_t hit_output = 0;
            std::uint64_t hit_camera = 0;
            std::uint64_t accepted = 0;
            double signal = 0.0;
            std::map<int, TelescopeEventAccumulator> telescopes;
            std::set<int> output_events;
            for (const auto& kv : summaries) {
                const auto& s = kv.second;
                if (showerEventFromOutputEvent(s.event_id,
                                               source_runtime_cfg.event_id_mode) != shower_event) {
                    continue;
                }
                input_bunches += s.input_bunches;
                input_photons += s.input_photons;
                blocked += s.blocked_by_obstruction;
                blocked_incoming += s.blocked_incoming;
                blocked_reflected += s.blocked_reflected;
                hit_output_before_obstruction += s.hit_output_before_obstruction;
                hit_output += s.hit_output_plane;
                hit_camera += s.hit_camera;
                accepted += s.accepted_camera;
                signal += s.weighted_signal;
                output_events.insert(s.event_id);
                auto& tel = telescopes[s.telescope_id];
                tel.telescope_id = s.telescope_id;
                tel.output_events.insert(s.event_id);
                tel.input_bunches += s.input_bunches;
                tel.input_photons += s.input_photons;
                tel.blocked_by_obstruction += s.blocked_by_obstruction;
                tel.blocked_incoming += s.blocked_incoming;
                tel.blocked_reflected += s.blocked_reflected;
                tel.hit_mirror_before_obstruction += s.hit_mirror_before_obstruction;
                tel.hit_output_before_obstruction += s.hit_output_before_obstruction;
                tel.hit_mirror += s.hit_mirror;
                tel.hit_output_plane += s.hit_output_plane;
                tel.hit_camera += s.hit_camera;
                tel.accepted_camera += s.accepted_camera;
                tel.lost_between_pixels += s.lost_between_pixels;
                tel.unique_pixels.insert(s.unique_pixels.begin(), s.unique_pixels.end());
                tel.weighted_signal += s.weighted_signal;
                tel.weighted_time_sum += s.weighted_time_sum;
                tel.weighted_time2_sum += s.weighted_time2_sum;
                tel.first_cherenkov_time_ns =
                    std::min(tel.first_cherenkov_time_ns, s.first_cherenkov_time_ns);
            }
            for (const auto& kv : telescopes) {
                printField("telescope",
                           formatTelescopeEventLine(kv.second, camera_cfg.enabled));
            }
            std::ostringstream value;
            value << "output_events=" << output_events.size()
                  << " telescopes=" << telescopes.size()
                  << " input_bunches=" << input_bunches
                  << " input_photons=" << doubleToString(input_photons, 3)
                  << " blocked=" << blocked
                  << " blocked_incoming=" << blocked_incoming
                  << " blocked_reflected=" << blocked_reflected
                  << " hit_output_before_obstruction=" << hit_output_before_obstruction
                  << " hit_output=" << hit_output;
            if (hit_output_before_obstruction > 0) {
                value << " output_transmission_after_obstruction="
                      << doubleToString(static_cast<double>(hit_output) /
                                        static_cast<double>(hit_output_before_obstruction), 6);
            }
            if (camera_cfg.enabled) {
                value << " hit_camera=" << hit_camera
                      << " accepted=" << accepted
                      << " pe=" << doubleToString(signal, 3);
            } else {
                value << " signal=" << doubleToString(signal, 3);
            }
            printField("event_total", value.str());
        };

        const PhotonTraceContext photon_trace_context{
            &tracer,
            &plane,
            &eff,
            &atmosphere,
            &obstruction,
            &camera,
            light_collector.get(),
            &sipm_cfg,
            &electronics,
            &detector_pipeline_cfg,
            &response_sampler,
            propagation_cfg.speed_of_light_m_per_ns,
            camera_cfg.enabled,
            true,
            nullptr,
            obstruction.mark_only,
            true,
            eventio_2d_backproject,
        };


        auto processBunch = [&](const PhotonBunch& raw_bunch) {
            std::chrono::steady_clock::time_point t_step;
            if (profile_cfg.enabled) {
                t_step = std::chrono::steady_clock::now();
            }
            const PhotonBunch bunch = transformEventIOBunchToTraceFrame(
                raw_bunch, telescope_cfg, metadata, source_runtime_cfg);
            if (profile_cfg.enabled) {
                addElapsed(profile_stats, &ProfileStats::transform_s, t_step);
            }
            const int shower_event =
                showerEventFromOutputEvent(bunch.event_id, source_runtime_cfg.event_id_mode);
            if (active_shower_event < 0) {
                active_shower_event = shower_event;
                printSection("Event shower_event=" + intToString(shower_event));
                printField("event_status", "started");
            } else if (shower_event != active_shower_event) {
                printRunningEventSummary(active_shower_event);
                active_shower_event = shower_event;
                printSection("Event shower_event=" + intToString(shower_event));
                printField("event_status", "started");
            }
            if (active_output_event < 0) {
                active_output_event = bunch.event_id;
            } else if (bunch.event_id != active_output_event) {
                flushOutputEvent(active_output_event);
                active_output_event = bunch.event_id;
            }
            auto& summary = summaries[{bunch.event_id, bunch.telescope_id}];
            summary.event_id = bunch.event_id;
            summary.telescope_id = bunch.telescope_id;
            summary.input_bunches += 1;
            summary.input_photons += bunch.multiplicity;
            const MirrorLayout& telescope_mirrors =
                telescope_optics.layoutFor(bunch.telescope_id);

            const Vec3 global_dir = sourceDirectionInWorld(
                raw_bunch, telescope_cfg, source_runtime_cfg.coordinate_frame);
            auto tracePhotonCandidate = [&](PhotonCandidate candidate) {
                ++photon_index;

                PhotonTraceProfile stage_profile;
                const PhotonTraceBunch trace_bunch{&bunch, &telescope_mirrors,
                                                   global_dir};
                const auto trace = tracePhoton(
                    photon_trace_context, trace_bunch, std::move(candidate),
                    profile_cfg.enabled ? &stage_profile : nullptr);
                if (profile_cfg.enabled) {
                    profile_stats.trace_to_plane_s +=
                        stage_profile.trace_to_plane_s;
                    profile_stats.obstruction_s += stage_profile.obstruction_s;
                    profile_stats.camera_response_s +=
                        stage_profile.camera_response_s;
                }
                const OpticalSurfaceHit& hit = trace.hit;
                const auto diagnose = [&](const char* label) {
                    if (mirror_diagnostic_out) {
                        writeCorsikaMirrorDiagnosticHit(
                            mirror_diagnostic_out, bunch, photon_index, hit,
                            label);
                    }
                };

                accumulateAtmosphereHistogram(
                    atmosphere_histogram,
                    atmosphere_histogram_cfg,
                    bunch.emission_altitude_km,
                    trace.pre_geometry.expected_weight_before_atmosphere,
                    trace.pre_geometry.expected_weight_after_atmosphere,
                    trace.pre_geometry.expected_weight_after_atmosphere);

                if (hit.hit_mirror) {
                    summary.hit_mirror_before_obstruction += 1;
                    if (hit.hit_surface) {
                        summary.hit_output_before_obstruction += 1;
                    }
                    if (trace.stage != PhotonTraceStage::BlockedIncoming) {
                        summary.hit_mirror += 1;
                    }
                }

                switch (trace.stage) {
                    case PhotonTraceStage::AbsorbedBeforeOptics:
                    case PhotonTraceStage::MissedMirror:
                        return;
                    case PhotonTraceStage::BlockedIncoming:
                        summary.blocked_by_obstruction += 1;
                        summary.blocked_incoming += 1;
                        diagnose("blocked_incoming");
                        return;
                    case PhotonTraceStage::ReflectedMissedOutput:
                        diagnose("reflected_missed_output");
                        return;
                    case PhotonTraceStage::BlockedReflected:
                        summary.blocked_by_obstruction += 1;
                        summary.blocked_reflected += 1;
                        diagnose("blocked_reflected");
                        return;
                    default:
                        break;
                }

                summary.hit_output_plane += 1;
                diagnose("reflected_to_output");

                std::chrono::steady_clock::time_point t_accumulate;
                if (profile_cfg.enabled) {
                    t_accumulate = std::chrono::steady_clock::now();
                }

                if (!camera_cfg.enabled) {
                    if (trace.stage != PhotonTraceStage::ReachedWhiteboard) {
                        return;
                    }
                    const double base_signal =
                        hit.weight * hit.relative_efficiency;
                    if (save_hdf5) {
                        whiteboard_hits.push_back(
                            makeWhiteboardHdf5Row(bunch, photon_index, hit));
                    }
                    if (save_csv) {
                        writeCorsikaWhiteboardHit(
                            whiteboard_out, bunch, photon_index, hit,
                            output_cfg.whiteboard_emitter_info);
                    }
                    summary.weighted_signal += base_signal;
                    summary.weighted_time_sum += base_signal * hit.time_ns;
                    summary.weighted_time2_sum +=
                        base_signal * hit.time_ns * hit.time_ns;
                    if (profile_cfg.enabled) {
                        addElapsed(profile_stats,
                                   &ProfileStats::whiteboard_accumulate_s,
                                   t_accumulate);
                    }
                    return;
                }

                appendCollectorDebugPhoton(collector_debug_rows,
                                           collector_debug_cfg,
                                           bunch.event_id,
                                           bunch.telescope_id,
                                           hit);
                if (trace.stage == PhotonTraceStage::LostBetweenPixels) {
                    summary.lost_between_pixels += 1;
                    return;
                }
                summary.hit_camera += 1;
                summary.unique_pixels.insert(hit.pixel_id);
                if (trace.stage == PhotonTraceStage::RejectedPostGeometry) {
                    return;
                }
                if (hit.accepted) {
                    summary.accepted_camera += 1;
                }

                const double signal = hit.weight * hit.relative_efficiency;
                summary.weighted_signal += signal;
                summary.weighted_time_sum += signal * hit.time_ns;
                summary.weighted_time2_sum += signal * hit.time_ns * hit.time_ns;
                summary.first_cherenkov_time_ns =
                    std::min(summary.first_cherenkov_time_ns, hit.time_ns);

                if (profile_cfg.enabled) {
                    t_accumulate = std::chrono::steady_clock::now();
                }
                accumulatePixelHit(pixels, bunch.event_id, bunch.telescope_id,
                                   hit, bunch.origin);
                accumulateWaveformHit(waveforms,
                                      raw_waveform_hits,
                                      waveform_cfg,
                                      detector_pipeline_cfg.enabled,
                                      bunch.event_id,
                                      bunch.telescope_id,
                                      hit,
                                      bunch.origin);
                if (profile_cfg.enabled) {
                    addElapsed(profile_stats,
                               &ProfileStats::camera_accumulate_s,
                               t_accumulate);
                }
            };

            const std::uint64_t candidate_count =
                response_sampler.candidateCount(bunch);
            for (std::uint64_t represented_index = 0;
                 represented_index < candidate_count;
                 ++represented_index) {
                tracePhotonCandidate(
                    response_sampler.candidate(bunch, represented_index));
            }
        };

        const auto t_stream_start = std::chrono::steady_clock::now();
        EventIOStreamStats stream_stats;
        if (photon_csv_mode) {
            PhotonCsvSource csv_source(
                buildPhotonCsvConfig(cfg, source_cfg, source_runtime_cfg));
            PhotonBunch csv_bunch;
            while (csv_source.next(csv_bunch)) {
                processBunch(csv_bunch);
                ++stream_stats.photon_bunches;
                if (csv_bunch.eventio_2d) {
                    ++stream_stats.photon_bunches_2d;
                } else {
                    ++stream_stats.photon_bunches_3d;
                }
            }
            std::ostringstream value;
            value << "photon_bunches=" << stream_stats.photon_bunches
                  << " photon_bunches_2d=" << stream_stats.photon_bunches_2d
                  << " photon_bunches_3d=" << stream_stats.photon_bunches_3d;
            printField("csv_load_done", value.str());
        } else {
            stream_stats = streamEventIOPhotonBunches(
                eventio_cfg,
                processBunch,
                [](const EventIOStreamProgress& progress) {
                    std::ostringstream value;
                    value << "photon_bunches=" << progress.photon_bunches
                          << " photon_bunches_2d=" << progress.photon_bunches_2d
                          << " photon_bunches_3d=" << progress.photon_bunches_3d
                          << " current_shower_event=" << progress.current_shower_event
                          << " elapsed_s=" << doubleToString(progress.elapsed_s, 3);
                    printField(progress.final ? "stream_done" : "stream_progress",
                               value.str());
                });
        }
        if (profile_cfg.enabled) {
            addElapsed(profile_stats, &ProfileStats::eventio_stream_s, t_stream_start);
        }
        flushOutputEvent(active_output_event);
        printRunningEventSummary(active_shower_event);
        const auto t_trace_done = std::chrono::steady_clock::now();

        if (camera_cfg.enabled && save_csv) {
            if (detector_pipeline_cfg.enabled) {
                writeElectronicsPixelCsv(output_cfg.pixel_csv,
                                         camera_pixel_id_axis,
                                         electronics_events,
                                         pixels);
                writeElectronicsTriggerCsv(output_cfg.trigger_csv,
                                           electronics_events);
                if (waveform_cfg.enabled) {
                    writeElectronicsWaveformCsv(output_cfg.waveform_csv,
                                                camera_pixel_id_axis,
                                                electronics_events);
                }
                writeElectronicsHitCsv(
                    output_cfg.primary_pe_csv,
                    output_cfg.fired_pe_csv,
                    camera_pixel_id_axis,
                    electronics_events,
                    detector_pipeline_cfg.save_primary_sequence,
                    detector_pipeline_cfg.save_fired_sequence);
            } else {
                writePixelCsv(output_cfg.pixel_csv, pixels);
            }
        }
        if (collector_debug_cfg.photon_output) {
            printSection("Collector debug output");
            printField("status", "writing collector debug photons");
            printField("path", collector_debug_cfg.photon_csv);
            writeCollectorDebugCsv(collector_debug_cfg, collector_debug_rows);
            printField("rows", intToString(static_cast<std::uint64_t>(
                                   collector_debug_rows.size())));
        }
        if (atmosphere_histogram_cfg.enabled) {
            printSection("Atmosphere histogram");
            printField("status", "writing height histogram");
            printField("path", atmosphere_histogram_cfg.csv_path);
            writeAtmosphereHistogramCsv(atmosphere_histogram_cfg, atmosphere_histogram);
        }
#ifdef LACT_HAS_HDF5
        if (save_hdf5) {
            printSection("HDF5 output");
            printField("status", "writing HDF5 trace file");
            printField("path", output_cfg.hdf5_path);
            const auto t_hdf5_start = std::chrono::steady_clock::now();
            writeNativeTraceHdf5(output_cfg,
                                 waveform_cfg,
                                 config_path,
                                 cfg,
                                 component_paths,
                                 source_runtime_cfg,
                                 telescope_cfg,
                                 metadata,
                                 camera,
                                 nominal_facets,
                                 sipm_cfg,
                                 electronics_cfg,
                                 detector_pipeline_cfg,
                                 efficiency_cfg,
                                 nsb_cfg,
                                 trigger_cfg,
                                 summaries,
                                 pixels,
                                 waveforms,
                                 raw_waveform_hits,
                                 electronics_events,
                                 whiteboard_hits);
            if (profile_cfg.enabled) {
                addElapsed(profile_stats, &ProfileStats::hdf5_write_s, t_hdf5_start);
            }
            printField("status", "HDF5 trace file written");
        }
#endif
#ifdef LACT_HAS_ROOT
        if (save_lact_root) {
            printField("status", "finalizing streaming lact_event ROOT file");
            lact_root_stream_writer->finish();
            printField("status", "lact_event ROOT file written");
        }
#endif
        if (save_csv) {
            writeSummaryCsv(output_cfg.summary_csv, summaries,
                            electronics_events);
        }
        const auto t_done = std::chrono::steady_clock::now();

        printSection("Input");
        printField("config", config_path);
        if (photon_csv_mode) {
            printField("csv_path", source_runtime_cfg.csv_path);
            if (!source_runtime_cfg.eventio_path.empty()) {
                printField("metadata_eventio_path",
                           source_runtime_cfg.eventio_path);
            }
        } else {
            printField("eventio_path", source_runtime_cfg.eventio_path);
        }
        printField("event_id_mode", source_runtime_cfg.event_id_mode);
        printField("output_format", output_cfg.format);
        printField("input_bunches", intToString(stream_stats.photon_bunches));
        printField("input_bunches_2d", intToString(stream_stats.photon_bunches_2d));
        printField("input_bunches_3d", intToString(stream_stats.photon_bunches_3d));
        printField("input_photon_format",
                   stream_stats.photon_bunches_2d > 0 && stream_stats.photon_bunches_3d > 0
                       ? "mixed_2d_3d"
                       : (stream_stats.photon_bunches_3d > 0 ? "3d" : "2d"));
        if (!source_runtime_cfg.eventio_path.empty()) {
            printField("telescopes_in_eventio",
                       intToString(metadata.telescopes.size()));
            printField("shower_events", intToString(metadata.events.size()));
        }
        printField("streams", intToString(summaries.size()));
        for (const auto& tel : metadata.telescopes) {
            std::ostringstream label;
            label << "tel[" << tel.telescope_id << "]";
            std::ostringstream value;
            value << "pos_m=(" << doubleToString(tel.x_m) << ", "
                  << doubleToString(tel.y_m) << ", "
                  << doubleToString(tel.z_m) << "), radius_m="
                  << doubleToString(tel.radius_m);
            printField(label.str(), value.str());
        }
        printEventSummary(summaries,
                          camera_cfg.enabled,
                          source_runtime_cfg.event_id_mode,
                          metadata);
        printSection("Detailed stream summary");
        if (save_hdf5) {
            printField("hdf5_path", output_cfg.hdf5_path);
        }
        if (save_lact_root) {
            printField("lact_root_path", output_cfg.lact_root_path);
        }
        if (save_csv) {
            printField("summary_csv", output_cfg.summary_csv);
            if (camera_cfg.enabled) {
                printField("pixel_csv", output_cfg.pixel_csv);
            } else {
                printField("hits_csv", output_cfg.hits_csv);
            }
        }
        printField("note", save_csv
                               ? "per-array/per-telescope details are written to CSV and/or HDF5, not expanded in the log"
                               : "per-array/per-telescope details are written to HDF5, not expanded in the log");
        printSection("Timing");
        printField("trace_time_s", doubleToString(elapsedSeconds(t_trace_start, t_trace_done)));
        printField("total_time_s", doubleToString(elapsedSeconds(t_start, t_done)));
        if (profile_cfg.enabled) {
            printProfileStats(profile_stats, elapsedSeconds(t_trace_start, t_trace_done));
        }
        printSection("Machine-readable summary");
        std::cout << "input_bunches=" << stream_stats.photon_bunches << "\n";
        std::cout << "shower_events=" << metadata.events.size() << "\n";
        std::cout << "streams=" << summaries.size() << "\n";
        std::cout << "camera_enabled=" << (camera_cfg.enabled ? 1 : 0) << "\n";
        std::cout << "camera_mode="
                  << (camera_cfg.enabled
                          ? camera_cfg.mode
                          : (camera_cfg.whiteboard
                                 ? "whiteboard"
                                 : "implicit_whiteboard_legacy"))
                  << "\n";
        if (save_hdf5) {
            std::cout << "hdf5_path=" << output_cfg.hdf5_path << "\n";
        }
        if (save_lact_root) {
            std::cout << "lact_root_path=" << output_cfg.lact_root_path << "\n";
        }
        if (save_csv) {
            std::cout << "summary_csv=" << output_cfg.summary_csv << "\n";
            if (camera_cfg.enabled) {
                std::cout << "pixel_csv=" << output_cfg.pixel_csv << "\n";
            } else {
                std::cout << "hits_csv=" << output_cfg.hits_csv << "\n";
            }
        }
        return 0;
#endif
    } catch (const std::exception& ex) {
        std::cerr << "run_corsika_trace error: " << ex.what() << "\n";
        return 1;
    }
}
