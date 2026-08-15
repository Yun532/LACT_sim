#include "app/CorsikaTraceConfig.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lact {

CorsikaTraceOutputConfig buildCorsikaTraceOutputConfig(
    const std::map<std::string, std::string>& cfg)
{
    CorsikaTraceOutputConfig out;
    out.hits_csv = getString(cfg, "output.hits_csv",
                             getString(cfg, "output.whiteboard_csv", out.hits_csv));
    out.pixel_csv = getString(cfg, "output.pixel_csv", out.pixel_csv);
    out.summary_csv = getString(cfg, "output.summary_csv", out.summary_csv);
    out.waveform_csv = getString(cfg, "output.waveform_csv", out.waveform_csv);
    out.trigger_csv = getString(cfg, "output.trigger_csv", out.trigger_csv);
    out.primary_pe_csv = getString(
        cfg, "output.primary_pe_csv", out.primary_pe_csv);
    out.fired_pe_csv = getString(
        cfg, "output.fired_pe_csv", out.fired_pe_csv);
    out.mirror_diagnostic_csv =
        getString(cfg, "output.mirror_diagnostic_csv",
                  out.mirror_diagnostic_csv);
    out.hdf5_path = getString(cfg, "output.hdf5_path",
                              getString(cfg, "output.h5_path", out.hdf5_path));
    out.lact_root_enabled =
        getBool(cfg, "output.lact_root_enabled", out.lact_root_enabled);
    out.lact_root_path =
        getString(cfg, "output.lact_root_path", out.lact_root_path);
    out.lact_profile =
        lowerCopy(trim(getString(cfg, "output.lact_profile", out.lact_profile)));
    out.format = lowerCopy(trim(getString(cfg, "output.format", out.format)));
    out.hdf5_storage =
        lowerCopy(trim(getString(cfg, "output.hdf5_storage", out.hdf5_storage)));
    out.hdf5_waveform_storage = lowerCopy(trim(getString(
        cfg, "output.hdf5_waveform_storage", out.hdf5_waveform_storage)));
    out.hdf5_write_components =
        getBool(cfg, "output.hdf5_write_components", out.hdf5_write_components);
    out.hdf5_write_waveforms =
        getBool(cfg, "output.hdf5_write_waveforms", out.hdf5_write_waveforms);
    out.lact_root_write_components =
        getBool(cfg, "output.lact_root_write_components",
                out.lact_root_write_components);
    out.lact_root_auto_flush_mb =
        getDouble(cfg, "output.lact_root_auto_flush_mb",
                  out.lact_root_auto_flush_mb);
    out.lact_root_flush_events =
        getInt(cfg, "output.lact_root_flush_events",
               out.lact_root_flush_events);
    out.save_only_triggered =
        getBool(cfg, "output.save_only_triggered", out.save_only_triggered);
    out.whiteboard_emitter_info =
        getBool(cfg, "output.whiteboard_emitter_info",
                getBool(cfg, "output.include_emitter_info",
                        out.whiteboard_emitter_info));
    out.write_pixel_time_stats =
        getBool(cfg,
                "output.write_pixel_time_stats",
                getBool(cfg,
                        "output.hdf5_write_pixel_time_stats",
                        out.write_pixel_time_stats));
    if (out.format.empty()) {
        out.format = "hdf5";
    }
    if (!(out.format == "hdf5" || out.format == "h5" ||
          out.format == "csv" || out.format == "both" ||
          out.format == "root" || out.format == "none")) {
        throw std::runtime_error("output.format must be hdf5, csv, both, root, or none");
    }
    if (!std::isfinite(out.lact_root_auto_flush_mb) ||
        out.lact_root_auto_flush_mb < 0.0) {
        throw std::runtime_error("output.lact_root_auto_flush_mb must be finite and >= 0");
    }
    if (out.lact_root_flush_events < 0) {
        throw std::runtime_error("output.lact_root_flush_events must be >= 0");
    }
    if (out.lact_profile.empty()) {
        out.lact_profile = "image_pe";
    }
    if (!(out.lact_profile == "image_pe" ||
          out.lact_profile == "timeseries_pe" ||
          out.lact_profile == "debug_full")) {
        throw std::runtime_error(
            "output.lact_profile must be image_pe, timeseries_pe, or debug_full");
    }
    if (out.hdf5_storage.empty()) {
        out.hdf5_storage = "dense";
    }
    if (!(out.hdf5_storage == "sparse" ||
          out.hdf5_storage == "dense" ||
          out.hdf5_storage == "both")) {
        throw std::runtime_error("output.hdf5_storage must be sparse, dense, or both");
    }
    if (!(out.hdf5_waveform_storage == "sparse" ||
          out.hdf5_waveform_storage == "dense")) {
        throw std::runtime_error(
            "output.hdf5_waveform_storage must be sparse or dense");
    }
    return out;
}

std::string normalizeWaveformTimeReference(std::string value)
{
    value = lowerCopy(trim(value));
    if (value.empty() || value == "absolute") {
        return "absolute";
    }
    if (value == "image_mean" || value == "image-mean") {
        return "image_mean";
    }
    if (value == "image_first" ||
        value == "image-first" ||
        value == "first_cherenkov" ||
        value == "first-cherenkov") {
        return "image_first";
    }
    throw std::runtime_error(
        "waveform.time_reference must be absolute, image_mean, or image_first");
}

WaveformOutputConfig buildWaveformOutputConfig(
    const std::map<std::string, std::string>& cfg)
{
    WaveformOutputConfig out;
    out.enabled = getBool(cfg, "waveform.enabled", out.enabled);
    out.source = lowerCopy(trim(getString(cfg, "waveform.source", out.source)));
    out.time_reference = normalizeWaveformTimeReference(
        getString(cfg, "waveform.time_reference", out.time_reference));
    out.time_bin_width_ns =
        getDouble(cfg, "waveform.time_bin_width_ns", out.time_bin_width_ns);
    out.time_window_start_ns =
        getDouble(cfg, "waveform.time_window_start_ns", out.time_window_start_ns);
    out.time_window_end_ns =
        getDouble(cfg, "waveform.time_window_end_ns", out.time_window_end_ns);
    if (out.source.empty()) {
        out.source = "none";
    }
    if (out.enabled) {
        if (!(out.source == "photon_count" || out.source == "pe" ||
              out.source == "electronics")) {
            throw std::runtime_error(
                "waveform.source must be photon_count, pe, or electronics when waveform.enabled=true");
        }
        if (out.time_bin_width_ns <= 0.0) {
            throw std::runtime_error("waveform.time_bin_width_ns must be > 0");
        }
        if (out.time_window_end_ns <= out.time_window_start_ns) {
            throw std::runtime_error(
                "waveform.time_window_end_ns must be greater than waveform.time_window_start_ns");
        }
        const double bin_count =
            (out.time_window_end_ns - out.time_window_start_ns) /
            out.time_bin_width_ns;
        const double rounded_bin_count = std::round(bin_count);
        const double bin_tolerance =
            1.0e-10 * std::max(1.0, std::abs(bin_count));
        if (std::abs(bin_count - rounded_bin_count) > bin_tolerance) {
            throw std::runtime_error(
                "waveform time interval must be an integer multiple of time_bin_width_ns");
        }
    }
    return out;
}

CollectorDebugConfig buildCollectorDebugConfig(
    const std::map<std::string, std::string>& cfg)
{
    CollectorDebugConfig out;
    out.photon_output = getBool(cfg, "collector.debug_photon_output", out.photon_output);
    out.photon_csv = getString(cfg, "collector.debug_photon_csv", out.photon_csv);
    out.max_photons = getUInt64(cfg, "collector.debug_max_photons", out.max_photons);
    return out;
}

ProfileConfig buildProfileConfig(const std::map<std::string, std::string>& cfg)
{
    ProfileConfig out;
    out.enabled = getBool(cfg, "profile.enabled", out.enabled);
    return out;
}

AtmosphereHistogramConfig buildAtmosphereHistogramConfig(
    const std::map<std::string, std::string>& cfg)
{
    AtmosphereHistogramConfig out;
    out.csv_path = getString(
        cfg,
        "atmosphere.height_histogram_csv",
        getString(cfg, "atmosphere.histogram_csv", ""));
    out.enabled = !isDisabledText(out.csv_path);
    out.min_altitude_km =
        getDouble(cfg, "atmosphere.histogram_min_altitude_km", out.min_altitude_km);
    out.max_altitude_km =
        getDouble(cfg, "atmosphere.histogram_max_altitude_km", out.max_altitude_km);
    out.bin_width_km =
        getDouble(cfg, "atmosphere.histogram_bin_width_km", out.bin_width_km);
    if (out.enabled) {
        if (out.bin_width_km <= 0.0) {
            throw std::runtime_error("atmosphere histogram bin width must be positive");
        }
        if (out.max_altitude_km <= out.min_altitude_km) {
            throw std::runtime_error("atmosphere histogram max altitude must exceed min altitude");
        }
    }
    return out;
}

} // namespace lact
