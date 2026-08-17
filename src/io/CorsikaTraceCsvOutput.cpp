#include "io/CorsikaTraceCsvOutput.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace lact {

void writeAtmosphereHistogramCsv(const AtmosphereHistogramConfig& cfg,
                                 const std::vector<AtmosphereHistogramBin>& bins)
{
    if (!cfg.enabled) {
        return;
    }
    const std::filesystem::path out_path(cfg.csv_path);
    if (out_path.has_parent_path()) {
        std::filesystem::create_directories(out_path.parent_path());
    }
    std::ofstream out(cfg.csv_path);
    if (!out) {
        throw std::runtime_error("failed to write atmosphere histogram CSV: " +
                                 cfg.csv_path);
    }
    out << "altitude_low_km,altitude_high_km,altitude_center_km,bunches,"
        << "before_weight,after_weight,theory_weight,transmission,"
        << "theory_transmission,relative_error\n";
    for (const auto& bin : bins) {
        const double center = 0.5 * (bin.low_km + bin.high_km);
        const double transmission =
            bin.before_weight > 0.0 ? bin.after_weight / bin.before_weight : 0.0;
        const double theory_transmission =
            bin.before_weight > 0.0 ? bin.theory_weight / bin.before_weight : 0.0;
        const double relative_error =
            bin.theory_weight > 0.0
                ? (bin.after_weight - bin.theory_weight) / bin.theory_weight
                : 0.0;
        out << bin.low_km << ','
            << bin.high_km << ','
            << center << ','
            << bin.bunches << ','
            << bin.before_weight << ','
            << bin.after_weight << ','
            << bin.theory_weight << ','
            << transmission << ','
            << theory_transmission << ','
            << relative_error << '\n';
    }
}

void writeCollectorDebugCsv(const CollectorDebugConfig& cfg,
                            const std::vector<CollectorDebugPhotonRow>& rows)
{
    if (!cfg.photon_output) {
        return;
    }
    const std::filesystem::path out_path(cfg.photon_csv);
    if (out_path.has_parent_path()) {
        std::filesystem::create_directories(out_path.parent_path());
    }
    std::ofstream out(cfg.photon_csv);
    if (!out) {
        throw std::runtime_error("failed to write collector debug CSV: " +
                                 cfg.photon_csv);
    }
    out << "event_id,telescope_id,pixel_id,accepted,collector_reflections,"
        << "collector_reflection_limit_reached,time_ns,"
        << "collector_path_length_m,collector_time_delay_ns,"
        << "wavelength_nm,photon_weight,"
        << "relative_efficiency,signal_weight,exit_x_m,exit_y_m,exit_z_m,"
        << "dir_u,dir_v,dir_w,sipm_gap_rejected,"
        << "sipm_channel_gap_rejected,"
        << "sipm_grid_column,sipm_grid_row,sipm_channel_id,"
        << "sipm_microcell_id\n";
    out << std::setprecision(10);
    for (const auto& row : rows) {
        out << row.event_id << ','
            << row.telescope_id << ','
            << row.pixel_id << ','
            << row.accepted << ','
            << row.collector_reflections << ','
            << row.collector_reflection_limit_reached << ','
            << row.time_ns << ','
            << row.collector_path_length_m << ','
            << row.collector_time_delay_ns << ','
            << row.wavelength_nm << ','
            << row.photon_weight << ','
            << row.relative_efficiency << ','
            << row.signal_weight << ','
            << row.exit_x_m << ','
            << row.exit_y_m << ','
            << row.exit_z_m << ','
            << row.dir_u << ','
            << row.dir_v << ','
            << row.dir_w << ','
            << row.sipm_gap_rejected << ','
            << row.sipm_channel_gap_rejected << ','
            << row.sipm_grid_column << ','
            << row.sipm_grid_row << ','
            << row.sipm_channel_id << ','
            << row.sipm_microcell_id << '\n';
    }
}

void writeCorsikaWhiteboardHeader(std::ofstream& ofs, bool include_emitter_info)
{
    ofs << std::setprecision(10);
    ofs << "event_id,telescope_id,photon_index,mirror_id,"
        << "surface_x_m,surface_y_m,surface_z_m,"
        << "mirror_x_m,mirror_y_m,mirror_z_m,"
        << "input_x_m,input_y_m,input_z_m,"
        << "input_dir_x,input_dir_y,input_dir_z,"
        << "u_m,v_m,dir_x,dir_y,dir_z,"
        << "time_ns,wavelength_nm,weight,relative_efficiency,"
        << "signal_weight";
    if (include_emitter_info) {
        ofs << ",has_emitter,emitter_mass_gev,emitter_charge,"
            << "emitter_energy_gev,emitter_time_ns";
    }
    ofs << "\n";
}

void writeCorsikaWhiteboardHit(std::ofstream& ofs,
                               const PhotonBunch& bunch,
                               std::uint64_t photon_index,
                               const OpticalSurfaceHit& hit,
                               bool include_emitter_info)
{
    const double signal = hit.weight * hit.relative_efficiency;
    ofs << bunch.event_id << ","
        << bunch.telescope_id << ","
        << photon_index << ","
        << hit.mirror_id << ","
        << hit.surface_point.x << ","
        << hit.surface_point.y << ","
        << hit.surface_point.z << ","
        << hit.mirror_point.x << ","
        << hit.mirror_point.y << ","
        << hit.mirror_point.z << ","
        << bunch.photon.pos.x << ","
        << bunch.photon.pos.y << ","
        << bunch.photon.pos.z << ","
        << bunch.photon.dir.x << ","
        << bunch.photon.dir.y << ","
        << bunch.photon.dir.z << ","
        << hit.u_m << ","
        << hit.v_m << ","
        << hit.out_dir.x << ","
        << hit.out_dir.y << ","
        << hit.out_dir.z << ","
        << hit.time_ns << ","
        << hit.wavelength_nm << ","
        << hit.weight << ","
        << hit.relative_efficiency << ","
        << signal;
    if (include_emitter_info) {
        ofs << ","
            << (bunch.has_emitter ? 1 : 0) << ","
            << bunch.emitter_mass_gev << ","
            << bunch.emitter_charge << ","
            << bunch.emitter_energy_gev << ","
            << bunch.emitter_time_ns;
    }
    ofs << "\n";
}

void writeCorsikaMirrorDiagnosticHeader(std::ofstream& ofs)
{
    ofs << std::setprecision(10);
    ofs << "event_id,telescope_id,photon_index,mirror_id,"
        << "mirror_x_m,mirror_y_m,mirror_z_m,"
        << "weight,relative_efficiency,status\n";
}

void writeCorsikaMirrorDiagnosticHit(std::ofstream& ofs,
                                     const PhotonBunch& bunch,
                                     std::uint64_t photon_index,
                                     const OpticalSurfaceHit& hit,
                                     const char* status)
{
    ofs << bunch.event_id << ","
        << bunch.telescope_id << ","
        << photon_index << ","
        << hit.mirror_id << ","
        << hit.mirror_point.x << ","
        << hit.mirror_point.y << ","
        << hit.mirror_point.z << ","
        << hit.weight << ","
        << hit.relative_efficiency << ","
        << status << "\n";
}

WhiteboardHdf5Row makeWhiteboardHdf5Row(const PhotonBunch& bunch,
                                        std::uint64_t photon_index,
                                        const OpticalSurfaceHit& hit)
{
    const double signal = hit.weight * hit.relative_efficiency;
    return WhiteboardHdf5Row{
        static_cast<std::int64_t>(bunch.event_id),
        static_cast<std::int32_t>(bunch.telescope_id),
        static_cast<std::int64_t>(photon_index),
        static_cast<std::int32_t>(hit.mirror_id),
        static_cast<float>(hit.surface_point.x),
        static_cast<float>(hit.surface_point.y),
        static_cast<float>(hit.surface_point.z),
        static_cast<float>(hit.u_m),
        static_cast<float>(hit.v_m),
        static_cast<float>(hit.out_dir.x),
        static_cast<float>(hit.out_dir.y),
        static_cast<float>(hit.out_dir.z),
        static_cast<float>(hit.time_ns),
        static_cast<float>(hit.wavelength_nm),
        static_cast<float>(hit.weight),
        static_cast<float>(hit.relative_efficiency),
        static_cast<float>(signal),
        static_cast<std::uint8_t>(bunch.has_emitter ? 1 : 0),
        static_cast<float>(bunch.emitter_mass_gev),
        static_cast<float>(bunch.emitter_charge),
        static_cast<float>(bunch.emitter_energy_gev),
        static_cast<float>(bunch.emitter_time_ns),
    };
}

void writeSummaryCsv(const std::string& path,
                     const std::map<SummaryKey, TraceSummary>& summaries,
                     const CameraElectronicsEventMap& electronics_events)
{
    const std::filesystem::path out_path(path);
    if (out_path.has_parent_path()) {
        std::filesystem::create_directories(out_path.parent_path());
    }
    std::ofstream ofs(path);
    if (!ofs) {
        throw std::runtime_error("failed to write summary CSV: " + path);
    }
    ofs << std::setprecision(17);
    ofs << "event_id,telescope_id,input_bunches,input_photons,"
        << "blocked_by_obstruction,blocked_incoming,blocked_reflected,"
        << "hit_mirror_before_obstruction,hit_output_before_obstruction,"
        << "hit_mirror,hit_output_plane,hit_camera,accepted_camera,"
        << "mirror_transmission_after_incoming_obstruction,"
        << "output_transmission_after_obstruction,"
        << "lost_between_pixels,unique_hit_pixels,pe,signal,"
        << "time_mean_ns,time_rms_ns,time_first_ns,reference_time_ns,"
        << "primary_cherenkov_pe,primary_nsb_pe,fired_cherenkov_pe,"
        << "fired_nsb_pe,fired_dark_pe,gap_lost_pe,saturation_lost_pe,"
        << "triggered,n_pixels_above_threshold,trigger_time_ns,"
        << "geometric_delay_ns,coincidence_time_ns,array_triggered\n";
    for (const auto& kv : summaries) {
        const auto& s = kv.second;
        const auto canonical = electronics_events.find(kv.first);
        if (canonical != electronics_events.end() &&
            !canonical->second.selected_for_output) {
            continue;
        }
        const double mean = s.weighted_signal > 0.0
            ? s.weighted_time_sum / s.weighted_signal
            : 0.0;
        const double var = s.weighted_signal > 0.0
            ? std::max(0.0, s.weighted_time2_sum / s.weighted_signal - mean * mean)
            : 0.0;
        const double mirror_transmission = s.hit_mirror_before_obstruction > 0
            ? static_cast<double>(s.hit_mirror) /
              static_cast<double>(s.hit_mirror_before_obstruction)
            : 0.0;
        const double output_transmission = s.hit_output_before_obstruction > 0
            ? static_cast<double>(s.hit_output_plane) /
              static_cast<double>(s.hit_output_before_obstruction)
            : 0.0;
        double primary_cherenkov_pe = 0.0;
        double primary_nsb_pe = 0.0;
        double fired_cherenkov_pe = 0.0;
        double fired_nsb_pe = 0.0;
        double fired_dark_pe = 0.0;
        double gap_lost_pe = 0.0;
        double saturation_lost_pe = 0.0;
        if (canonical != electronics_events.end()) {
            for (const auto& pixel : canonical->second.detector.pixels) {
                primary_cherenkov_pe += pixel.primary_cherenkov_pe;
                primary_nsb_pe += pixel.primary_nsb_pe;
                fired_cherenkov_pe += pixel.fired_cherenkov_pe;
                fired_nsb_pe += pixel.fired_nsb_pe;
                fired_dark_pe += pixel.fired_dark_pe;
                gap_lost_pe += pixel.gap_lost_pe;
                saturation_lost_pe += pixel.saturation_lost_pe;
            }
        }
        const bool triggered = canonical != electronics_events.end() &&
            canonical->second.detector.camera_trigger.triggered;
        const double nan = std::numeric_limits<double>::quiet_NaN();
        ofs << s.event_id << ","
            << s.telescope_id << ","
            << s.input_bunches << ","
            << s.input_photons << ","
            << s.blocked_by_obstruction << ","
            << s.blocked_incoming << ","
            << s.blocked_reflected << ","
            << s.hit_mirror_before_obstruction << ","
            << s.hit_output_before_obstruction << ","
            << s.hit_mirror << ","
            << s.hit_output_plane << ","
            << s.hit_camera << ","
            << s.accepted_camera << ","
            << mirror_transmission << ","
            << output_transmission << ","
            << s.lost_between_pixels << ","
            << s.unique_pixels.size() << ","
            << s.weighted_signal << ","
            << s.weighted_signal << ","
            << mean << ","
            << std::sqrt(var) << ","
            << (std::isfinite(s.first_cherenkov_time_ns)
                    ? s.first_cherenkov_time_ns : nan) << ","
            << (canonical != electronics_events.end()
                    ? canonical->second.reference_time_ns : 0.0) << ","
            << primary_cherenkov_pe << ","
            << primary_nsb_pe << ","
            << fired_cherenkov_pe << ","
            << fired_nsb_pe << ","
            << fired_dark_pe << ","
            << gap_lost_pe << ","
            << saturation_lost_pe << ","
            << (triggered ? 1 : 0) << ","
            << (canonical != electronics_events.end()
                    ? canonical->second.detector.camera_trigger
                          .max_pixels_above_threshold : 0) << ","
            << (canonical != electronics_events.end()
                    ? canonical->second.trigger_time_ns : nan) << ","
            << (canonical != electronics_events.end()
                    ? canonical->second.geometric_delay_ns : nan) << ","
            << (canonical != electronics_events.end()
                    ? canonical->second.coincidence_time_ns : nan) << ","
            << (canonical != electronics_events.end() &&
                    canonical->second.array_triggered ? 1 : 0) << "\n";
    }
}

void ensureCsvParent(const std::string& path)
{
    const std::filesystem::path output_path(path);
    if (output_path.has_parent_path()) {
        std::filesystem::create_directories(output_path.parent_path());
    }
}

void writeElectronicsPixelCsv(
    const std::string& path,
    const std::vector<int>& pixel_id_axis,
    const CameraElectronicsEventMap& events,
    const std::map<PixelKey, PixelAccumulator>& optical_pixels)
{
    ensureCsvParent(path);
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("failed to write electronics pixel CSV: " + path);
    }
    output << std::setprecision(17);
    output << "event_id,telescope_id,pixel_id,reference_time_ns,"
           << "primary_cherenkov_pe,primary_nsb_pe,primary_dark_pe,"
           << "fired_cherenkov_pe,fired_nsb_pe,fired_dark_pe,fired_pe,"
           << "gap_lost_pe,saturation_lost_pe,time_mean_ns,time_rms_ns\n";
    for (const auto& item : events) {
        const auto& event = item.second;
        if (!event.selected_for_output) continue;
        if (event.detector.pixels.size() != pixel_id_axis.size()) {
            throw std::runtime_error("electronics CSV pixel axis mismatch");
        }
        for (std::size_t column = 0; column < pixel_id_axis.size(); ++column) {
            const int pixel_id = pixel_id_axis[column];
            const auto& pixel = event.detector.pixels[column];
            const double fired_pe = pixel.fired_cherenkov_pe +
                pixel.fired_nsb_pe + pixel.fired_dark_pe;
            const double primary_pe = pixel.primary_cherenkov_pe +
                pixel.primary_nsb_pe + pixel.primary_dark_pe;
            if (fired_pe <= 0.0 && primary_pe <= 0.0) continue;
            const auto optical = optical_pixels.find(
                PixelKey{event.event_id, event.telescope_id, pixel_id});
            double mean = std::numeric_limits<double>::quiet_NaN();
            double rms = std::numeric_limits<double>::quiet_NaN();
            if (optical != optical_pixels.end() && optical->second.signal > 0.0) {
                const auto& source = optical->second;
                mean = source.time_sum / source.signal;
                const double variance = std::max(
                    0.0, source.time2_sum / source.signal - mean * mean);
                rms = std::sqrt(variance);
            }
            output << event.event_id << ',' << event.telescope_id << ','
                   << pixel_id << ',' << event.reference_time_ns << ','
                   << pixel.primary_cherenkov_pe << ','
                   << pixel.primary_nsb_pe << ',' << pixel.primary_dark_pe << ','
                   << pixel.fired_cherenkov_pe << ',' << pixel.fired_nsb_pe << ','
                   << pixel.fired_dark_pe << ',' << fired_pe << ','
                   << pixel.gap_lost_pe << ',' << pixel.saturation_lost_pe << ','
                   << mean << ',' << rms << '\n';
        }
    }
}

void writeElectronicsWaveformCsv(
    const std::string& path,
    const std::vector<int>& pixel_id_axis,
    const CameraElectronicsEventMap& events)
{
    ensureCsvParent(path);
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("failed to write electronics waveform CSV: " + path);
    }
    output << std::setprecision(17);
    output << "event_id,telescope_id,reference_time_ns,pixel_id,time_bin,"
           << "time_center_ns,global_time_ns,sample_value,sample_unit\n";
    for (const auto& item : events) {
        const auto& event = item.second;
        if (!event.selected_for_output) continue;
        const auto& result = event.detector;
        const std::size_t expected_samples = result.n_pixels * result.n_samples;
        if (result.waveform.size() != expected_samples ||
            result.time_centers_ns.size() != result.n_samples ||
            pixel_id_axis.size() != result.n_pixels) {
            throw std::runtime_error(
                "cannot write electronics waveform CSV: detector waveform "
                "dimensions are inconsistent (is single-p.e. shaping enabled?)");
        }
        for (std::size_t column = 0; column < result.n_pixels; ++column) {
            for (std::size_t bin = 0; bin < result.n_samples; ++bin) {
                const float value = static_cast<float>(
                    result.waveform[column * result.n_samples + bin]);
                if (value == 0.0) continue;
                const double relative_time = result.time_centers_ns[bin];
                output << event.event_id << ',' << event.telescope_id << ','
                       << event.reference_time_ns << ','
                       << pixel_id_axis.at(column) << ',' << bin << ','
                       << relative_time << ','
                       << event.reference_time_ns + relative_time << ','
                       << value << ',' << result.sample_unit << '\n';
            }
        }
    }
}

void writeElectronicsTriggerCsv(const std::string& path,
                                const CameraElectronicsEventMap& events)
{
    ensureCsvParent(path);
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("failed to write electronics trigger CSV: " + path);
    }
    output << std::setprecision(17);
    output << "event_id,telescope_id,reference_time_ns,triggered,"
           << "n_pixels_above_threshold,trigger_time_ns,geometric_delay_ns,"
           << "coincidence_time_ns,array_triggered\n";
    for (const auto& item : events) {
        const auto& event = item.second;
        if (!event.selected_for_output) continue;
        output << event.event_id << ',' << event.telescope_id << ','
               << event.reference_time_ns << ','
               << (event.detector.camera_trigger.triggered ? 1 : 0) << ','
               << event.detector.camera_trigger.max_pixels_above_threshold << ','
               << event.trigger_time_ns << ',' << event.geometric_delay_ns << ','
               << event.coincidence_time_ns << ','
               << (event.array_triggered ? 1 : 0) << '\n';
    }
}

void writeElectronicsHitCsv(const std::string& primary_path,
                            const std::string& fired_path,
                            const std::vector<int>& pixel_id_axis,
                            const CameraElectronicsEventMap& events,
                            bool write_primary,
                            bool write_fired)
{
    std::ofstream primary;
    std::ofstream fired;
    if (write_primary) {
        ensureCsvParent(primary_path);
        primary.open(primary_path);
        if (!primary) throw std::runtime_error(
            "failed to write primary p.e. CSV: " + primary_path);
        primary << std::setprecision(17)
                << "event_id,telescope_id,reference_time_ns,pixel_id,time_ns,"
                << "global_time_ns,sensor_x_m,sensor_y_m,wavelength_nm,"
                << "primary_pe,origin\n";
    }
    if (write_fired) {
        ensureCsvParent(fired_path);
        fired.open(fired_path);
        if (!fired) throw std::runtime_error(
            "failed to write fired p.e. CSV: " + fired_path);
        fired << std::setprecision(17)
              << "event_id,telescope_id,reference_time_ns,pixel_id,time_ns,"
              << "global_time_ns,channel_id,microcell_id,fired_pe,"
              << "charge_factor,time_jitter_ns,waveform_time_ns,"
              << "global_waveform_time_ns,origin\n";
    }
    for (const auto& item : events) {
        const auto& event = item.second;
        if (!event.selected_for_output) continue;
        if (write_primary) {
            for (const auto& hit : event.detector.primary_hits) {
                primary << hit.event_id << ',' << hit.telescope_id << ','
                        << event.reference_time_ns << ','
                        << pixel_id_axis.at(static_cast<std::size_t>(hit.pixel_id))
                        << ',' << hit.time_ns << ','
                        << event.reference_time_ns + hit.time_ns << ','
                        << hit.sensor_x_m << ',' << hit.sensor_y_m << ','
                        << hit.wavelength_nm << ',' << hit.primary_pe << ','
                        << static_cast<int>(hit.origin) << '\n';
            }
        }
        if (write_fired) {
            for (const auto& hit : event.detector.fired_hits) {
                fired << hit.event_id << ',' << hit.telescope_id << ','
                      << event.reference_time_ns << ','
                      << pixel_id_axis.at(static_cast<std::size_t>(hit.pixel_id))
                      << ',' << hit.time_ns << ','
                      << event.reference_time_ns + hit.time_ns << ','
                      << hit.channel_id << ',' << hit.microcell_id << ','
                      << hit.fired_pe << ',' << hit.charge_factor << ','
                      << hit.time_jitter_ns << ','
                      << hit.time_ns + hit.time_jitter_ns << ','
                      << event.reference_time_ns + hit.time_ns +
                             hit.time_jitter_ns
                      << ',' << static_cast<int>(hit.origin)
                      << '\n';
            }
        }
    }
}


} // namespace lact
