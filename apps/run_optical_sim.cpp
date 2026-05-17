#include "app/OpticalSimCommon.hpp"

using namespace lact;

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: run_optical_sim <config.txt>\n";
        return 2;
    }

    try {
        const auto t_start = std::chrono::steady_clock::now();
        auto main_cfg = readKeyValueConfig(argv[1]);
        ComponentConfigPaths component_paths;
        auto cfg = expandConfig(main_cfg, argv[1], component_paths);
        const auto t_config_read = std::chrono::steady_clock::now();

        SyntheticPhotonConfig source_cfg = buildSourceConfig(cfg);
        SourceRuntimeConfig source_runtime_cfg = buildSourceRuntimeConfig(cfg);
        if (source_runtime_cfg.use_eventio) {
            throw std::runtime_error(
                "run_optical_sim no longer accepts source.mode=EventIO. "
                "Use run_corsika_trace for CORSIKA/EventIO input so the "
                "CORSIKA NWU coordinate transform and event metadata are handled consistently.");
        }
        TelescopeConfig telescope_cfg = buildTelescopeConfig(cfg);
        TelescopeFrame telescope_frame = buildTelescopeFrame(telescope_cfg);
        std::vector<MirrorFacet> facets = buildFacetsFromConfig(cfg);
        ErrorConfig error_cfg = buildErrorConfig(cfg);
        ObstructionMask obstruction = buildObstructionMask(cfg);
        applyStructuralDeformation(facets, error_cfg, telescope_cfg);
        OutputPlane plane = buildOutputPlane(cfg);
        CameraConfig camera_cfg = buildCameraConfig(cfg);
        SipmConfig sipm_cfg = buildSipmConfig(cfg);
        ElectronicsConfig electronics_cfg = buildElectronicsConfig(cfg);
        ElectronicsResponse electronics(electronics_cfg);
        CameraGeometry camera = buildCameraGeometry(camera_cfg);
        auto light_collector = buildLightCollector(camera_cfg, camera);
        applyFacetErrors(facets, error_cfg);
        applyTelescopeFrame(facets, plane, telescope_frame);
        MirrorLayout mirrors = makeMirrorLayoutFromFacets(facets);
        OpticalEfficiencyConfig efficiency_cfg = buildEfficiencyConfig(cfg);
        PropagationConfig propagation_cfg = buildPropagationConfig(cfg);
        OpticalEfficiency eff(efficiency_cfg);
        std::string output_csv = getString(cfg, "output.csv", "surface_hits.csv");
        std::string output_pixel_csv = getString(cfg, "output.pixel_csv", "camera_pixel_image.csv");
        const std::string output_mode = lowerCopy(trim(getString(cfg, "output.mode", "hits")));
        const bool save_pixel_csv = camera_cfg.enabled &&
            (output_mode == "pixel" || output_mode == "pixels" ||
             output_mode == "both" || cfg.find("output.pixel_csv") != cfg.end());
        const bool save_hits_csv = !(output_mode == "pixel" || output_mode == "pixels");
        const auto t_setup_done = std::chrono::steady_clock::now();

        const std::string mirror_mode = lowerCopy(getString(cfg, "mirror.mode", "generated"));
        const std::string mirror_csv = getString(cfg, "mirror.csv_path", "");

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "========================================\n";
        std::cout << "LACT optical simulation\n";
        std::cout << "========================================\n";

        printSection("Configuration files");
        printField("main", argv[1]);
        if (!component_paths.telescope.empty()) {
            printField("telescope", component_paths.telescope);
        }
        if (!component_paths.mirror.empty()) {
            printField("mirror", component_paths.mirror);
        }
        if (!component_paths.source.empty()) {
            printField("source", component_paths.source);
        }
        if (!component_paths.output.empty()) {
            printField("output", component_paths.output);
        }
        if (!component_paths.camera.empty()) {
            printField("camera", component_paths.camera);
        }
        if (!component_paths.sipm.empty()) {
            printField("sipm", component_paths.sipm);
        }
        if (!component_paths.electronics.empty()) {
            printField("electronics", component_paths.electronics);
        }
        if (!component_paths.efficiency.empty()) {
            printField("efficiency", component_paths.efficiency);
        }
        if (!component_paths.atmosphere.empty()) {
            printField("atmosphere", component_paths.atmosphere);
        }
        if (!component_paths.error.empty()) {
            printField("error", component_paths.error);
        }
        if (!component_paths.obstruction.empty()) {
            printField("obstruction", component_paths.obstruction);
        }

        printSection("Telescope");
        printField("id", intToString(static_cast<std::uint64_t>(telescope_cfg.id)));
        printField("name", telescope_cfg.name);
        printField("position_m", vec3ToString(telescope_cfg.position_m));
        printField("pointing_az_deg", doubleToString(telescope_cfg.pointing_az_deg));
        printField("pointing_el_deg", doubleToString(telescope_cfg.pointing_el_deg));
        printField("focal_length_m", doubleToString(telescope_cfg.focal_length_m));
        printField("coordinate_system", telescope_cfg.coordinate_system);
        printField("coordinate_transform", "local telescope frame -> global frame");
        printField("frame_x_axis", vec3ToString(telescope_frame.x_axis));
        printField("frame_y_axis", vec3ToString(telescope_frame.y_axis));
        printField("frame_z_axis", vec3ToString(telescope_frame.z_axis));

        printSection("Mirror");
        printField("mode", mirror_mode);
        if (!mirror_csv.empty()) {
            printField("csv", mirror_csv);
        }
        if (mirror_mode == "elevation_series" || mirror_mode == "series") {
            printField("series_elevation_deg",
                       doubleToString(getDouble(cfg, "mirror.series_elevation_deg",
                                                telescope_cfg.pointing_el_deg)));
            printField("series_angles_deg",
                       getString(cfg, "mirror.series_angles_deg", ""));
            printField("series_csv_pattern",
                       getString(cfg, "mirror.series_csv_pattern", ""));
        }
        printField("facets", intToString(mirrors.size()));

        printSection("Source");
        printField("mode", source_runtime_cfg.use_photon_csv
                               ? "PhotonCsv"
                               : (source_runtime_cfg.use_eventio
                                      ? "EventIO"
                                      : sourceModeName(source_cfg.mode)));
        if (!source_runtime_cfg.use_photon_csv && !source_runtime_cfg.use_eventio) {
            printField("n_bunches",
                       intToString(static_cast<std::uint64_t>(source_cfg.n_bunches)));
        }
        printField("multiplicity", doubleToString(source_cfg.multiplicity));
        if (source_runtime_cfg.use_photon_csv || source_runtime_cfg.use_eventio) {
            if (source_runtime_cfg.use_photon_csv) {
                printField("csv", source_runtime_cfg.csv_path);
            } else {
                printField("eventio_path", source_runtime_cfg.eventio_path);
                printField("event_id_mode", source_runtime_cfg.event_id_mode);
            }
            printField("local_telescope_frame",
                       source_runtime_cfg.csv_local_telescope_frame ? "true" : "false");
            printField("filter_telescope_id",
                       source_runtime_cfg.filter_telescope_id
                           ? intToString(source_runtime_cfg.selected_telescope_id)
                           : "off");
            printField("filter_event_id",
                       source_runtime_cfg.filter_event_id
                           ? intToString(source_runtime_cfg.selected_event_id)
                           : "off");
            printField("filter_shower_event_id",
                       source_runtime_cfg.filter_shower_event_id
                           ? intToString(source_runtime_cfg.selected_shower_event_id)
                           : "off");
            printField("max_shower_events",
                       source_runtime_cfg.max_shower_events > 0
                           ? intToString(source_runtime_cfg.max_shower_events)
                           : "off");
            printField("default_wavelength_nm", doubleToString(source_cfg.wavelength_nm));
            printField("default_weight", doubleToString(source_cfg.photon_weight));
        } else if (source_cfg.mode == SyntheticMode::ParallelBeam) {
            printField("source_plane_z", doubleToString(source_cfg.source_plane_z));
            printField("beam_radius_m", doubleToString(source_cfg.beam_radius_m));
            if (cfg.find("source.beam_theta_deg") != cfg.end() ||
                cfg.find("source.beam_phi_deg") != cfg.end()) {
                printField("beam_theta_deg",
                           getString(cfg, "source.beam_theta_deg", "0"));
                printField("beam_phi_deg",
                           getString(cfg, "source.beam_phi_deg", "0"));
            }
            printField("beam_direction", vec3ToString(source_cfg.beam_direction));
        } else if (source_cfg.mode == SyntheticMode::PointSource) {
            printField("source_position", vec3ToString(source_cfg.source_position));
            printField("aperture_z", doubleToString(source_cfg.aperture_z));
            printField("aperture_radius_m", doubleToString(source_cfg.aperture_radius_m));
            const double ideal_entrance_plane_z = 0.0;
            const double ideal_distance_m =
                std::abs(source_cfg.source_position.z - ideal_entrance_plane_z);
            const double actual_distance_m =
                std::abs(source_cfg.source_position.z - source_cfg.aperture_z);
            printField("ideal_entrance_plane_z", doubleToString(ideal_entrance_plane_z));
            printField("ideal_entrance_plane_distance_m", doubleToString(ideal_distance_m));
            printField("configured_target_plane_distance_m", doubleToString(actual_distance_m));
            if (std::abs(actual_distance_m - ideal_distance_m) > 1e-9) {
                printField("warning",
                           "PointSource target plane differs from ideal local entrance plane (z=0)");
            }
        }
        printField("random_seed", intToString(source_cfg.random_seed));

        printSection("Output plane");
        printField("point", vec3ToString(plane.point));
        printField("normal", vec3ToString(plane.normal));
        printField("mode", output_mode);
        if (save_hits_csv) {
            printField("hits_csv", output_csv);
        }
        if (save_pixel_csv) {
            printField("pixel_csv", output_pixel_csv);
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
                printField("collector_exit_m",
                           doubleToString(camera_cfg.collector_exit_size_m));
                printField("collector_height_m",
                           doubleToString(camera_cfg.collector_height_m));
            }
        } else {
            printField("mode", "whiteboard only");
        }

        printSection("SiPM");
        printField("size_m", doubleToString(sipm_cfg.size_m));
        printField("pde", factorDescription(efficiency_cfg.sipm_pde));

        printSection("Electronics");
        printField("response", "reserved; SiPM PDE is handled by sipm.pde");

        printSection("Efficiency");
        printField("constant_scale", doubleToString(efficiency_cfg.constant_scale));
        printField("mirror_reflectivity",
                   factorDescription(efficiency_cfg.mirror_reflectivity));
        printField("filter_transmission",
                   factorDescription(efficiency_cfg.filter_transmission));
        printField("atmosphere",
                   factorDescription(efficiency_cfg.atmosphere_transmission));
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
        printField("facet_normal_sigma_deg",
                   doubleToString(error_cfg.facet_normal_sigma_deg));
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
            printField("mark_only", obstruction.mark_only ? "true" : "false");
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
        printField("optics", "ideal reflection only");
        printField("speed_of_light_m/ns",
                   doubleToString(propagation_cfg.speed_of_light_m_per_ns, 9));
        printField("not included",
                   light_collector
                       ? "mirror roughness, misalignment, camera electronics"
                       : "mirror roughness, misalignment, collector, camera electronics");

        printSection("Run");
        printField("setup_time_s", doubleToString(elapsedSeconds(t_start, t_setup_done)));
        printField("config_read_s", doubleToString(elapsedSeconds(t_start, t_config_read)));
        printField("status", "tracing started");

        std::unique_ptr<PhotonSource> source;
        std::size_t reserve_hits = 0;
        if (source_runtime_cfg.use_photon_csv) {
            auto csv_cfg = buildPhotonCsvConfig(cfg, source_cfg, source_runtime_cfg);
            auto csv_source = std::make_unique<PhotonCsvSource>(csv_cfg);
            reserve_hits = csv_source->size();
            source = std::move(csv_source);
        } else if (source_runtime_cfg.use_eventio) {
            throw std::runtime_error(
                "run_optical_sim no longer accepts source.mode=EventIO; use run_corsika_trace.");
        } else {
            reserve_hits = static_cast<std::size_t>(std::max(0, source_cfg.n_bunches));
            source = std::make_unique<SyntheticPhotonSource>(source_cfg);
        }
        OpticalTracer tracer(propagation_cfg.speed_of_light_m_per_ns,
                             error_cfg.reflect_direction_sigma_deg * DEG_TO_RAD,
                             error_cfg.random_seed);

        std::vector<OpticalSurfaceHit> hits;
        if (save_hits_csv) {
            hits.reserve(reserve_hits);
        }
        std::map<PixelKey, PixelAccumulator> pixels;

        PhotonBunch bunch;
        int n_total = 0;
        int n_hit_mirror_before_obstruction = 0;
        int n_hit_surface_before_obstruction = 0;
        int n_hit_mirror = 0;
        int n_hit_surface = 0;
        int n_hit_camera = 0;
        int n_accepted = 0;
        int n_blocked = 0;
        int n_blocked_incoming = 0;
        int n_blocked_reflected = 0;
        std::set<int> unique_hit_pixels;
        double sum_w = 0.0;
        double sum_r2 = 0.0;

        const auto t_trace_start = std::chrono::steady_clock::now();
        while (source->next(bunch)) {
            ++n_total;

            Photon photon = bunch.photon;
            photon.normalizeDirection();
            photon.weight *= bunch.multiplicity;
            if ((!source_runtime_cfg.use_photon_csv && !source_runtime_cfg.use_eventio) ||
                source_runtime_cfg.csv_local_telescope_frame) {
                applyTelescopeFrame(photon, telescope_frame);
            }

            OpticalSurfaceHit hit = tracer.traceToPlane(photon, mirrors, plane, eff);
            if (hit.hit_mirror) {
                ++n_hit_mirror_before_obstruction;
                if (hit.hit_surface) {
                    ++n_hit_surface_before_obstruction;
                }
                hit.obstruction_blocked_incoming =
                    incomingSegmentBlockedByObstruction(photon.pos, hit.mirror_point,
                                                        obstruction, &telescope_frame);
                if (hit.obstruction_blocked_incoming) {
                    ++n_blocked;
                    ++n_blocked_incoming;
                    if (!obstruction.mark_only) {
                        continue;
                    }
                } else {
                    ++n_hit_mirror;
                }
            }
            if (hit.hit_surface) {
                hit.obstruction_blocked_reflected =
                    segmentBlockedByObstruction(hit.mirror_point, hit.surface_point,
                                                obstruction, &telescope_frame);
                if (hit.obstruction_blocked_reflected) {
                    ++n_blocked_reflected;
                    if (!hit.obstruction_blocked_incoming) {
                        ++n_blocked;
                    }
                    if (!obstruction.mark_only) {
                        continue;
                    }
                }
                hit.obstruction_blocked = hit.obstruction_blocked_incoming ||
                                          hit.obstruction_blocked_reflected;
                const bool physically_reaches_output = !hit.obstruction_blocked;
                if (camera_cfg.enabled && physically_reaches_output) {
                    applyCameraResponse(camera, light_collector.get(), plane, sipm_cfg,
                                        electronics, hit);
                    if (hit.hit_camera) {
                        ++n_hit_camera;
                        unique_hit_pixels.insert(hit.pixel_id);
                    }
                    if (hit.accepted) {
                        ++n_accepted;
                    }
                    if (save_pixel_csv) {
                        accumulatePixelHit(pixels, bunch.event_id, bunch.telescope_id, hit);
                    }
                }
                if (physically_reaches_output) {
                    ++n_hit_surface;
                }
                if (save_hits_csv) {
                    hits.push_back(hit);
                }

                if (physically_reaches_output) {
                    double w = hit.weight * hit.relative_efficiency;
                    double r2 = hit.u_m * hit.u_m + hit.v_m * hit.v_m;
                    sum_w += w;
                    sum_r2 += w * r2;
                }
            }
        }
        const auto t_trace_done = std::chrono::steady_clock::now();

        if (save_hits_csv && !writeSurfaceHitsCSV(output_csv, hits)) {
            throw std::runtime_error("failed to write output CSV: " + output_csv);
        }
        if (save_pixel_csv) {
            writePixelCsv(output_pixel_csv, pixels);
        }
        const auto t_write_done = std::chrono::steady_clock::now();

        const double mirror_fraction = n_total > 0
            ? static_cast<double>(n_hit_mirror) / static_cast<double>(n_total)
            : 0.0;
        const double mirror_before_fraction = n_total > 0
            ? static_cast<double>(n_hit_mirror_before_obstruction) / static_cast<double>(n_total)
            : 0.0;
        const double surface_fraction = n_total > 0
            ? static_cast<double>(n_hit_surface) / static_cast<double>(n_total)
            : 0.0;
        const double surface_before_fraction = n_total > 0
            ? static_cast<double>(n_hit_surface_before_obstruction) / static_cast<double>(n_total)
            : 0.0;
        const double mirror_obstruction_transmission = n_hit_mirror_before_obstruction > 0
            ? static_cast<double>(n_hit_mirror) /
              static_cast<double>(n_hit_mirror_before_obstruction)
            : 0.0;
        const double output_obstruction_transmission = n_hit_surface_before_obstruction > 0
            ? static_cast<double>(n_hit_surface) /
              static_cast<double>(n_hit_surface_before_obstruction)
            : 0.0;
        const bool has_sampling_area =
            !source_runtime_cfg.use_photon_csv && !source_runtime_cfg.use_eventio;
        const double sampling_radius_m = source_cfg.mode == SyntheticMode::ParallelBeam
            ? source_cfg.beam_radius_m
            : source_cfg.aperture_radius_m;
        const double source_sampling_area_m2 = has_sampling_area
            ? std::acos(-1.0) * sampling_radius_m * sampling_radius_m
            : std::nan("");
        const double mirror_area_before_obstruction_m2 =
            source_sampling_area_m2 * mirror_before_fraction;
        const double mirror_area_after_incoming_obstruction_m2 =
            source_sampling_area_m2 * mirror_fraction;
        const double output_area_before_obstruction_m2 =
            source_sampling_area_m2 * surface_before_fraction;
        const double output_area_after_obstruction_m2 =
            source_sampling_area_m2 * surface_fraction;
        const double output_area_loss_from_obstruction_m2 =
            output_area_before_obstruction_m2 - output_area_after_obstruction_m2;
        const double camera_fraction = n_total > 0
            ? static_cast<double>(n_hit_camera) / static_cast<double>(n_total)
            : 0.0;
        const double accepted_fraction = n_total > 0
            ? static_cast<double>(n_accepted) / static_cast<double>(n_total)
            : 0.0;
        const double camera_fill_fraction = n_hit_surface > 0
            ? static_cast<double>(n_hit_camera) / static_cast<double>(n_hit_surface)
            : 0.0;
        const double rms = sum_w > 0.0 ? std::sqrt(sum_r2 / sum_w) : std::nan("");

        printSection("Results");
        printField("total_photons", intToString(static_cast<std::uint64_t>(n_total)));
        printField("blocked_by_obstruction",
                   intToString(static_cast<std::uint64_t>(n_blocked)));
        printField("blocked_incoming",
                   intToString(static_cast<std::uint64_t>(n_blocked_incoming)));
        printField("blocked_reflected",
                   intToString(static_cast<std::uint64_t>(n_blocked_reflected)));
        printField("hit_mirror_before_obstruction",
                   intToString(static_cast<std::uint64_t>(n_hit_mirror_before_obstruction)));
        printField("hit_output_before_obstruction",
                   intToString(static_cast<std::uint64_t>(n_hit_surface_before_obstruction)));
        printField("hit_mirror", intToString(static_cast<std::uint64_t>(n_hit_mirror)));
        printField("hit_output_plane", intToString(static_cast<std::uint64_t>(n_hit_surface)));
        if (camera_cfg.enabled) {
            printField("hit_camera", intToString(static_cast<std::uint64_t>(n_hit_camera)));
            printField("accepted_camera", intToString(static_cast<std::uint64_t>(n_accepted)));
            printField("lost_between_pixels",
                       intToString(static_cast<std::uint64_t>(n_hit_surface - n_hit_camera)));
            printField("unique_hit_pixels",
                       intToString(static_cast<std::uint64_t>(unique_hit_pixels.size())));
        }
        printField("hit_mirror_before_obstruction_fraction", doubleToString(mirror_before_fraction));
        printField("hit_mirror_fraction", doubleToString(mirror_fraction));
        printField("hit_output_before_obstruction_fraction", doubleToString(surface_before_fraction));
        printField("hit_output_fraction", doubleToString(surface_fraction));
        printField("mirror_transmission_after_incoming_obstruction",
                   doubleToString(mirror_obstruction_transmission));
        printField("mirror_loss_fraction_from_incoming_obstruction",
                   doubleToString(1.0 - mirror_obstruction_transmission));
        printField("output_transmission_after_obstruction",
                   doubleToString(output_obstruction_transmission));
        printField("output_loss_fraction_from_obstruction",
                   doubleToString(1.0 - output_obstruction_transmission));
        if (has_sampling_area) {
            printField("source_sampling_area_m2", doubleToString(source_sampling_area_m2));
            printField("mirror_collecting_area_before_obstruction_m2",
                       doubleToString(mirror_area_before_obstruction_m2));
            printField("mirror_collecting_area_after_incoming_obstruction_m2",
                       doubleToString(mirror_area_after_incoming_obstruction_m2));
            printField("output_collecting_area_before_obstruction_m2",
                       doubleToString(output_area_before_obstruction_m2));
            printField("output_collecting_area_after_obstruction_m2",
                       doubleToString(output_area_after_obstruction_m2));
            printField("output_collecting_area_loss_from_obstruction_m2",
                       doubleToString(output_area_loss_from_obstruction_m2));
        }
        if (camera_cfg.enabled) {
            printField("hit_camera_fraction", doubleToString(camera_fraction));
            printField("accepted_camera_fraction", doubleToString(accepted_fraction));
            printField("camera_fill_fraction", doubleToString(camera_fill_fraction));
        }
        printField("weighted_spot_rms_m", doubleToString(rms));
        printField("weighted_spot_rms_mm", doubleToString(rms * 1000.0));
        if (save_hits_csv) {
            printField("output_csv", output_csv);
        }
        if (save_pixel_csv) {
            printField("pixel_csv", output_pixel_csv);
        }

        printSection("Timing");
        printField("trace_time_s", doubleToString(elapsedSeconds(t_trace_start, t_trace_done)));
        printField("csv_write_time_s", doubleToString(elapsedSeconds(t_trace_done, t_write_done)));
        printField("total_time_s", doubleToString(elapsedSeconds(t_start, t_write_done)));

        printSection("Machine-readable summary");
        std::cout << "mirror_facets=" << mirrors.size() << "\n";
        std::cout << "total_photons=" << n_total << "\n";
        std::cout << "blocked_by_obstruction=" << n_blocked << "\n";
        std::cout << "blocked_incoming=" << n_blocked_incoming << "\n";
        std::cout << "blocked_reflected=" << n_blocked_reflected << "\n";
        std::cout << "hit_mirror_before_obstruction="
                  << n_hit_mirror_before_obstruction << "\n";
        std::cout << "hit_output_before_obstruction="
                  << n_hit_surface_before_obstruction << "\n";
        std::cout << "hit_mirror=" << n_hit_mirror << "\n";
        std::cout << "hit_output_plane=" << n_hit_surface << "\n";
        std::cout << "mirror_transmission_after_incoming_obstruction="
                  << mirror_obstruction_transmission << "\n";
        std::cout << "output_transmission_after_obstruction="
                  << output_obstruction_transmission << "\n";
        if (has_sampling_area) {
            std::cout << "source_sampling_area_m2=" << source_sampling_area_m2 << "\n";
            std::cout << "mirror_collecting_area_before_obstruction_m2="
                      << mirror_area_before_obstruction_m2 << "\n";
            std::cout << "mirror_collecting_area_after_incoming_obstruction_m2="
                      << mirror_area_after_incoming_obstruction_m2 << "\n";
            std::cout << "output_collecting_area_before_obstruction_m2="
                      << output_area_before_obstruction_m2 << "\n";
            std::cout << "output_collecting_area_after_obstruction_m2="
                      << output_area_after_obstruction_m2 << "\n";
            std::cout << "output_collecting_area_loss_from_obstruction_m2="
                      << output_area_loss_from_obstruction_m2 << "\n";
        }
        if (camera_cfg.enabled) {
            std::cout << "hit_camera=" << n_hit_camera << "\n";
            std::cout << "accepted_camera=" << n_accepted << "\n";
            std::cout << "lost_between_pixels=" << (n_hit_surface - n_hit_camera) << "\n";
            std::cout << "unique_hit_pixels=" << unique_hit_pixels.size() << "\n";
        }
        std::cout << "weighted_spot_rms_m=" << rms << "\n";
        if (save_hits_csv) {
            std::cout << "output_csv=" << output_csv << "\n";
        }
        if (save_pixel_csv) {
            std::cout << "pixel_csv=" << output_pixel_csv << "\n";
        }
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "run_optical_sim error: " << ex.what() << "\n";
        return 1;
    }
}
