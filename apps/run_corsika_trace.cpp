#include "app/OpticalSimCommon.hpp"
#include "app/PhotonResponseSampler.hpp"
#include "io/EventIOArrayTiming.hpp"
#include "app/TelescopeOpticsCache.hpp"
#include "app/TriggerResponse.hpp"
#include "core/Sha256.hpp"
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

struct CollectorDebugConfig {
    bool photon_output = false;
    std::string photon_csv = "collector_debug_photons.csv";
    std::uint64_t max_photons = 100000;
};

struct WavelengthRange {
    double min_nm = 0.0;
    double max_nm = 0.0;
};

struct ProfileConfig {
    bool enabled = false;
};

struct AtmosphereHistogramConfig {
    bool enabled = false;
    std::string csv_path = "atmosphere_height_histogram.csv";
    double min_altitude_km = 4.4;
    double max_altitude_km = 100.0;
    double bin_width_km = 1.0;
};

struct AtmosphereHistogramBin {
    double low_km = 0.0;
    double high_km = 0.0;
    std::uint64_t bunches = 0;
    double before_weight = 0.0;
    double after_weight = 0.0;
    double theory_weight = 0.0;
};

struct ProfileStats {
    double eventio_stream_s = 0.0;
    double transform_s = 0.0;
    double trace_to_plane_s = 0.0;
    double obstruction_s = 0.0;
    double camera_response_s = 0.0;
    double whiteboard_accumulate_s = 0.0;
    double camera_accumulate_s = 0.0;
    double hdf5_write_s = 0.0;
};

struct TelescopeEventAccumulator {
    int telescope_id = 0;
    std::set<int> output_events;
    std::uint64_t input_bunches = 0;
    double input_photons = 0.0;
    std::uint64_t blocked_by_obstruction = 0;
    std::uint64_t blocked_incoming = 0;
    std::uint64_t blocked_reflected = 0;
    std::uint64_t hit_mirror_before_obstruction = 0;
    std::uint64_t hit_output_before_obstruction = 0;
    std::uint64_t hit_mirror = 0;
    std::uint64_t hit_output_plane = 0;
    std::uint64_t hit_camera = 0;
    std::uint64_t accepted_camera = 0;
    std::uint64_t lost_between_pixels = 0;
    std::set<int> unique_pixels;
    double weighted_signal = 0.0;
    double weighted_time_sum = 0.0;
    double weighted_time2_sum = 0.0;
    double first_cherenkov_time_ns = std::numeric_limits<double>::infinity();
};

struct CollectorDebugPhotonRow {
    int event_id = 0;
    int telescope_id = 0;
    int pixel_id = -1;
    int accepted = 0;
    int collector_reflections = 0;
    int collector_reflection_limit_reached = 0;
    double time_ns = 0.0;
    double collector_path_length_m = 0.0;
    double collector_time_delay_ns = 0.0;
    double wavelength_nm = 0.0;
    double photon_weight = 0.0;
    double relative_efficiency = 0.0;
    double signal_weight = 0.0;
    double exit_x_m = 0.0;
    double exit_y_m = 0.0;
    double exit_z_m = 0.0;
    double dir_u = 0.0;
    double dir_v = 0.0;
    double dir_w = 0.0;
};

struct WhiteboardHdf5Row {
    std::int64_t event_id;
    std::int32_t telescope_id;
    std::int64_t photon_index;
    std::int32_t mirror_id;
    float surface_x_m;
    float surface_y_m;
    float surface_z_m;
    float u_m;
    float v_m;
    float dir_x;
    float dir_y;
    float dir_z;
    float time_ns;
    float wavelength_nm;
    float weight;
    float relative_efficiency;
    float signal_weight;
    std::uint8_t has_emitter;
    float emitter_mass_gev;
    float emitter_charge;
    float emitter_energy_gev;
    float emitter_time_ns;
};

PhotonBunch transformEventIOBunchToTraceFrame(
    const PhotonBunch& input,
    const TelescopeConfig& telescope_cfg,
    const EventIOMetadata& metadata,
    const SourceRuntimeConfig& source_runtime_cfg)
{
    TelescopeConfig telescope = telescope_cfg;
    const std::string frame_name = normalizeSourceCoordinateFrame(
        source_runtime_cfg.coordinate_frame);
    if ((frame_name == "corsika_nwu_global" ||
         frame_name == "enu_east_global" ||
         frame_name == "lact_generic_global") &&
        source_runtime_cfg.use_eventio_telescope_position) {
        if (auto tel = metadata.telescopeById(input.telescope_id)) {
            if (frame_name == "enu_east_global") {
                // EventIO telescope metadata is CORSIKA NWU. Convert it to
                // the ENU frame expected by this PhotonCsv input.
                telescope.position_m = {-tel->y_m, tel->x_m, tel->z_m};
            } else {
                telescope.position_m = {tel->x_m, tel->y_m, tel->z_m};
            }
        }
    }
    PhotonBunch out = transformBunchToTelescopeLocal(input, telescope, frame_name);
    if (source_runtime_cfg.use_eventio || out.eventio_2d) {
        applyEventIOReferenceZOffset(
            out, source_runtime_cfg.eventio_reference_z_m);
    }
    return out;
}

int showerEventFromOutputEvent(int event_id, const std::string& event_id_mode)
{
    const std::string mode = lowerCopy(trim(event_id_mode));
    if (mode == "event_array100" || mode == "runid") {
        return event_id / 100;
    }
    return event_id;
}

int arrayIdFromOutputEvent(int event_id, const std::string& event_id_mode)
{
    const std::string mode = lowerCopy(trim(event_id_mode));
    if (mode == "event_array100" || mode == "runid") {
        return event_id % 100;
    }
    return 0;
}

double mirrorFrontReferenceZ(const MirrorLayout& mirrors)
{
    double z = -std::numeric_limits<double>::infinity();
    for (const auto& tile : mirrors.tiles()) {
        z = std::max(z, tile.center.z);
    }
    return std::isfinite(z) ? z : -16.0;
}

bool shouldBackprojectEventIO2d(const SourceRuntimeConfig& source_runtime_cfg)
{
    const std::string mode = lowerCopy(trim(source_runtime_cfg.eventio_2d_plane_mode));
    if (mode == "forward") {
        return false;
    }
    // A 2D EventIO bunch supplies an anchor on the unperturbed photon line,
    // not a creation point. In both auto and explicit backproject modes the
    // mirror may therefore lie at either sign of t relative to that anchor.
    return true;
}

struct OutputEventMetadata {
    int event_id = 0;
    int shower_event = 0;
    int array_id = 0;
    double energy_gev = 0.0;
    double core_x_north_m = 0.0;
    double core_y_west_m = 0.0;
    double azimuth_north_to_east_deg = 0.0;
    double array_time_offset_ns = 0.0;
    double area_weight_m2 = 0.0;
    bool has_explicit_area_weight = false;
    bool found = false;
    bool used_array_offset = false;
};

OutputEventMetadata outputEventMetadata(int event_id,
                                        const std::string& event_id_mode,
                                        const EventIOMetadata& metadata)
{
    OutputEventMetadata out;
    out.event_id = event_id;
    const auto identity = metadata.output_event_identity.find(event_id);
    if (identity != metadata.output_event_identity.end()) {
        out.shower_event = identity->second.first;
        out.array_id = identity->second.second;
    } else {
        out.shower_event = showerEventFromOutputEvent(event_id, event_id_mode);
        out.array_id = arrayIdFromOutputEvent(event_id, event_id_mode);
    }

    auto event_it = std::find_if(
        metadata.events.begin(), metadata.events.end(),
        [&out](const EventIOEventHeader& event) {
            return event.shower_event_id == out.shower_event;
        });
    if (event_it == metadata.events.end()) {
        return out;
    }

    out.found = true;
    out.energy_gev = event_it->energy_gev;
    out.core_x_north_m = event_it->core_x_m;
    out.core_y_west_m = event_it->core_y_m;
    out.azimuth_north_to_east_deg = event_it->azimuth_north_to_east_deg;

    if (auto offsets = metadata.arrayOffsetsForShower(out.shower_event)) {
        out.array_time_offset_ns = offsets->time_offset_ns;
        const std::size_t offset_index = static_cast<std::size_t>(out.array_id);
        if (out.array_id >= 0 && offset_index < offsets->x_m.size() &&
            offset_index < offsets->y_m.size()) {
            // MC_TELOFF stores the offset of the detector array with respect
            // to the shower core.  The core position in the telescope/input
            // array frame is therefore the opposite vector.
            out.core_x_north_m = -offsets->x_m[offset_index];
            out.core_y_west_m = -offsets->y_m[offset_index];
            out.used_array_offset = true;
        }
        if (out.array_id >= 0 && offset_index < offsets->weight.size()) {
            out.area_weight_m2 = offsets->weight[offset_index];
            out.has_explicit_area_weight = offsets->has_explicit_weights;
        }
    }
    return out;
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

std::optional<WavelengthRange> wavelengthRangeFromInputCard(
    const EventIOMetadata& metadata)
{
    for (const auto& line : metadata.input_lines) {
        std::istringstream iss(line);
        std::string key;
        if (!(iss >> key)) {
            continue;
        }
        if (lowerCopy(key) != "cwavlg") {
            continue;
        }
        WavelengthRange range;
        if (iss >> range.min_nm >> range.max_nm &&
            std::isfinite(range.min_nm) && std::isfinite(range.max_nm) &&
            range.min_nm > 0.0 && range.max_nm > range.min_nm) {
            return range;
        }
    }
    return std::nullopt;
}

bool hasExplicitMissingWavelengthRange(const std::map<std::string, std::string>& cfg)
{
    return cfg.find("source.missing_wavelength_min_nm") != cfg.end() ||
           cfg.find("source.missing_wavelength_max_nm") != cfg.end() ||
           cfg.find("source.wavelength_min_nm") != cfg.end() ||
           cfg.find("source.wavelength_max_nm") != cfg.end();
}

void applyEventIOWavelengthMetadata(EventIOPhotonConfig& eventio_cfg,
                                    const EventIOMetadata& metadata,
                                    const std::map<std::string, std::string>& cfg)
{
    if (hasExplicitMissingWavelengthRange(cfg)) {
        return;
    }
    if (auto range = wavelengthRangeFromInputCard(metadata)) {
        eventio_cfg.missing_wavelength_min_nm = range->min_nm;
        eventio_cfg.missing_wavelength_max_nm = range->max_nm;
    }
}

void applyEventIOAtmosphereMetadata(EventIOPhotonConfig& eventio_cfg,
                                    const EventIOMetadata& metadata)
{
    if (std::isfinite(metadata.observation_altitude_m)) {
        eventio_cfg.observation_altitude_km =
            metadata.observation_altitude_m * 1.0e-3;
    }
}

CorsikaTraceOutputConfig buildCorsikaTraceOutputConfig(
    const std::map<std::string, std::string>& cfg)
{
    CorsikaTraceOutputConfig out;
    out.hits_csv = getString(cfg, "output.hits_csv",
                             getString(cfg, "output.whiteboard_csv", out.hits_csv));
    out.pixel_csv = getString(cfg, "output.pixel_csv", out.pixel_csv);
    out.summary_csv = getString(cfg, "output.summary_csv", out.summary_csv);
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

std::vector<AtmosphereHistogramBin> makeAtmosphereHistogramBins(
    const AtmosphereHistogramConfig& cfg)
{
    std::vector<AtmosphereHistogramBin> bins;
    if (!cfg.enabled) {
        return bins;
    }
    for (double low = cfg.min_altitude_km; low < cfg.max_altitude_km;
         low += cfg.bin_width_km) {
        bins.push_back({low, std::min(cfg.max_altitude_km, low + cfg.bin_width_km)});
    }
    return bins;
}

void accumulateAtmosphereHistogram(std::vector<AtmosphereHistogramBin>& bins,
                                   const AtmosphereHistogramConfig& cfg,
                                   double altitude_km,
                                   double before_weight,
                                   double after_weight,
                                   double theory_weight)
{
    if (!cfg.enabled || bins.empty() || !std::isfinite(altitude_km)) {
        return;
    }
    if (altitude_km < cfg.min_altitude_km || altitude_km >= cfg.max_altitude_km) {
        return;
    }
    std::size_t index = static_cast<std::size_t>(
        (altitude_km - cfg.min_altitude_km) / cfg.bin_width_km);
    if (index >= bins.size()) {
        index = bins.size() - 1;
    }
    auto& bin = bins[index];
    bin.bunches += 1;
    bin.before_weight += before_weight;
    bin.after_weight += after_weight;
    bin.theory_weight += theory_weight;
}

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

void addElapsed(ProfileStats& stats,
                double ProfileStats::*field,
                const std::chrono::steady_clock::time_point& start)
{
    stats.*field +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

bool outputWantsCsv(const CorsikaTraceOutputConfig& cfg)
{
    return cfg.format == "csv" || cfg.format == "both";
}

bool outputWantsHdf5(const CorsikaTraceOutputConfig& cfg)
{
    return cfg.format == "hdf5" || cfg.format == "h5" || cfg.format == "both";
}

bool outputWantsLactRoot(const CorsikaTraceOutputConfig& cfg)
{
    return cfg.lact_root_enabled;
}

std::size_t waveformBinCount(const WaveformOutputConfig& cfg)
{
    if (!cfg.enabled) {
        return 0;
    }
    const double span = cfg.time_window_end_ns - cfg.time_window_start_ns;
    return static_cast<std::size_t>(std::ceil(span / cfg.time_bin_width_ns));
}

int waveformBinForTime(const WaveformOutputConfig& cfg, double time_ns)
{
    if (!cfg.enabled ||
        time_ns < cfg.time_window_start_ns ||
        time_ns >= cfg.time_window_end_ns) {
        return -1;
    }
    const auto bin = static_cast<int>(
        std::floor((time_ns - cfg.time_window_start_ns) / cfg.time_bin_width_ns));
    const auto n_bins = static_cast<int>(waveformBinCount(cfg));
    return bin >= 0 && bin < n_bins ? bin : -1;
}

bool waveformUsesImageMeanReference(const WaveformOutputConfig& cfg)
{
    return cfg.enabled && cfg.time_reference == "image_mean";
}

bool waveformUsesImageReference(const WaveformOutputConfig& cfg)
{
    return cfg.enabled &&
        (cfg.time_reference == "image_mean" || cfg.time_reference == "image_first");
}

void accumulateWaveformHit(std::map<WaveformKey, WaveformPixelAccumulator>& waveform,
                           std::vector<RawWaveformHit>& raw_waveform_hits,
                           const WaveformOutputConfig& cfg,
                           bool capture_detector_hit,
                           int event_id,
                           int telescope_id,
                           const OpticalSurfaceHit& hit)
{
    if (!hit.hit_camera || hit.pixel_id < 0) {
        return;
    }
    const double pe = hit.weight * hit.relative_efficiency;
    if (capture_detector_hit) {
        raw_waveform_hits.push_back(RawWaveformHit{
            event_id,
            telescope_id,
            hit.pixel_id,
            hit.time_ns,
            1,
            pe,
            hit.collector_exit_x_m,
            hit.collector_exit_y_m,
            hit.wavelength_nm,
            electronics::HitOrigin::Cherenkov,
        });
    }
    if (!cfg.enabled) {
        return;
    }
    if (waveformUsesImageReference(cfg)) {
        if (!capture_detector_hit) {
            raw_waveform_hits.push_back(RawWaveformHit{
                event_id,
                telescope_id,
                hit.pixel_id,
                hit.time_ns,
                1,
                pe,
                hit.collector_exit_x_m,
                hit.collector_exit_y_m,
                hit.wavelength_nm,
                electronics::HitOrigin::Cherenkov,
            });
        }
        return;
    }
    const int bin = waveformBinForTime(cfg, hit.time_ns);
    if (bin < 0) {
        return;
    }
    auto& acc = waveform[{event_id, telescope_id, hit.pixel_id, bin}];
    acc.event_id = event_id;
    acc.telescope_id = telescope_id;
    acc.pixel_id = hit.pixel_id;
    acc.time_bin = bin;
    acc.photon_count += 1;
    acc.pe += pe;
}

void appendCollectorDebugPhoton(std::vector<CollectorDebugPhotonRow>& rows,
                                const CollectorDebugConfig& cfg,
                                int event_id,
                                int telescope_id,
                                const OpticalSurfaceHit& hit)
{
    if (!cfg.photon_output || !hit.collector_enabled ||
        rows.size() >= cfg.max_photons) {
        return;
    }
    rows.push_back(CollectorDebugPhotonRow{
        event_id,
        telescope_id,
        hit.pixel_id,
        hit.accepted ? 1 : 0,
        hit.collector_reflections,
        hit.collector_reflection_limit_reached ? 1 : 0,
        hit.time_ns,
        hit.collector_path_length_m,
        hit.collector_time_delay_ns,
        hit.wavelength_nm,
        hit.weight,
        hit.relative_efficiency,
        hit.weight * hit.relative_efficiency,
        hit.collector_exit_x_m,
        hit.collector_exit_y_m,
        hit.collector_exit_z_m,
        hit.collector_dir_u,
        hit.collector_dir_v,
        hit.collector_dir_w,
    });
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
        << "dir_u,dir_v,dir_w\n";
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
            << row.dir_w << '\n';
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
                     const std::map<SummaryKey, TraceSummary>& summaries)
{
    const std::filesystem::path out_path(path);
    if (out_path.has_parent_path()) {
        std::filesystem::create_directories(out_path.parent_path());
    }
    std::ofstream ofs(path);
    if (!ofs) {
        throw std::runtime_error("failed to write summary CSV: " + path);
    }
    ofs << std::setprecision(10);
    ofs << "event_id,telescope_id,input_bunches,input_photons,"
        << "blocked_by_obstruction,blocked_incoming,blocked_reflected,"
        << "hit_mirror_before_obstruction,hit_output_before_obstruction,"
        << "hit_mirror,hit_output_plane,hit_camera,accepted_camera,"
        << "mirror_transmission_after_incoming_obstruction,"
        << "output_transmission_after_obstruction,"
        << "lost_between_pixels,unique_hit_pixels,pe,signal,"
        << "time_mean_ns,time_rms_ns\n";
    for (const auto& kv : summaries) {
        const auto& s = kv.second;
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
            << std::sqrt(var) << "\n";
    }
}

#ifdef LACT_HAS_HDF5
int hdf5PixelShapeCode(PixelShape shape)
{
    if (shape == PixelShape::Square) return 1;
    if (shape == PixelShape::Hexagonal) return 2;
    if (shape == PixelShape::Circular) return 3;
    return 0;
}

int hdf5FacetShapeCode(ApertureShape shape)
{
    if (shape == ApertureShape::Square) return 1;
    if (shape == ApertureShape::Hexagon) return 2;
    if (shape == ApertureShape::Circular) return 3;
    return 0;
}

void h5Check(herr_t status, const std::string& message)
{
    if (status < 0) {
        throw std::runtime_error("HDF5 write failed: " + message);
    }
}

void writeStringAttribute(hid_t object, const std::string& name, const std::string& value)
{
    hid_t space = H5Screate(H5S_SCALAR);
    if (space < 0) {
        throw std::runtime_error("HDF5 write failed: create scalar dataspace");
    }
    hid_t type = H5Tcopy(H5T_C_S1);
    if (type < 0) {
        H5Sclose(space);
        throw std::runtime_error("HDF5 write failed: create string type");
    }
    H5Tset_size(type, std::max<std::size_t>(1, value.size() + 1));
    H5Tset_strpad(type, H5T_STR_NULLTERM);
    hid_t attr = H5Acreate2(object, name.c_str(), type, space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) {
        H5Tclose(type);
        H5Sclose(space);
        throw std::runtime_error("HDF5 write failed: create attribute " + name);
    }
    h5Check(H5Awrite(attr, type, value.c_str()), "write attribute " + name);
    H5Aclose(attr);
    H5Tclose(type);
    H5Sclose(space);
}

void writeStringDataset(hid_t group, const std::string& name, const std::string& value)
{
    hid_t space = H5Screate(H5S_SCALAR);
    if (space < 0) {
        throw std::runtime_error("HDF5 write failed: create scalar dataspace");
    }
    hid_t type = H5Tcopy(H5T_C_S1);
    if (type < 0) {
        H5Sclose(space);
        throw std::runtime_error("HDF5 write failed: create string type");
    }
    H5Tset_size(type, std::max<std::size_t>(1, value.size() + 1));
    H5Tset_strpad(type, H5T_STR_NULLTERM);
    hid_t ds = H5Dcreate2(group, name.c_str(), type, space,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (ds < 0) {
        H5Tclose(type);
        H5Sclose(space);
        throw std::runtime_error("HDF5 write failed: create dataset " + name);
    }
    h5Check(H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, value.c_str()),
            "write dataset " + name);
    H5Dclose(ds);
    H5Tclose(type);
    H5Sclose(space);
}

std::string readTextIfExists(const std::string& path)
{
    if (path.empty() || !std::filesystem::exists(path)) {
        return "";
    }
    std::ifstream ifs(path);
    std::ostringstream ss;
    ss << ifs.rdbuf();
    return ss.str();
}

template <typename Row>
void writeCompound1D(hid_t group,
                     const std::string& name,
                     hid_t type,
                     const std::vector<Row>& rows)
{
    hsize_t dims[1] = {rows.size()};
    hid_t space = H5Screate_simple(1, dims, nullptr);
    if (space < 0) {
        throw std::runtime_error("HDF5 write failed: create dataspace " + name);
    }
    hid_t ds = H5Dcreate2(group, name.c_str(), type, space,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (ds < 0) {
        H5Sclose(space);
        throw std::runtime_error("HDF5 write failed: create dataset " + name);
    }
    if (!rows.empty()) {
        h5Check(H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, rows.data()),
                "write dataset " + name);
    }
    H5Dclose(ds);
    H5Sclose(space);
}

template <typename T>
void writePlain1D(hid_t group,
                  const std::string& name,
                  hid_t type,
                  const std::vector<T>& values)
{
    hsize_t dims[1] = {values.size()};
    hid_t space = H5Screate_simple(1, dims, nullptr);
    if (space < 0) {
        throw std::runtime_error("HDF5 write failed: create dataspace " + name);
    }
    hid_t ds = H5Dcreate2(group, name.c_str(), type, space,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (ds < 0) {
        H5Sclose(space);
        throw std::runtime_error("HDF5 write failed: create dataset " + name);
    }
    if (!values.empty()) {
        h5Check(H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()),
                "write dataset " + name);
    }
    H5Dclose(ds);
    H5Sclose(space);
}

template <typename T>
void writePlain2D(hid_t group,
                  const std::string& name,
                  hid_t type,
                  const std::vector<T>& values,
                  hsize_t n_rows,
                  hsize_t n_cols)
{
    hsize_t dims[2] = {n_rows, n_cols};
    hid_t space = H5Screate_simple(2, dims, nullptr);
    if (space < 0) {
        throw std::runtime_error("HDF5 write failed: create dataspace " + name);
    }
    hid_t ds = H5Dcreate2(group, name.c_str(), type, space,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (ds < 0) {
        H5Sclose(space);
        throw std::runtime_error("HDF5 write failed: create dataset " + name);
    }
    if (!values.empty()) {
        h5Check(H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()),
                "write dataset " + name);
    }
    H5Dclose(ds);
    H5Sclose(space);
}

template <typename T>
void writePlain3D(hid_t group,
                  const std::string& name,
                  hid_t type,
                  const std::vector<T>& values,
                  hsize_t n0,
                  hsize_t n1,
                  hsize_t n2)
{
    hsize_t dims[3] = {n0, n1, n2};
    hid_t space = H5Screate_simple(3, dims, nullptr);
    if (space < 0) {
        throw std::runtime_error("HDF5 write failed: create dataspace " + name);
    }
    hid_t ds = H5Dcreate2(group, name.c_str(), type, space,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (ds < 0) {
        H5Sclose(space);
        throw std::runtime_error("HDF5 write failed: create dataset " + name);
    }
    if (!values.empty()) {
        h5Check(H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()),
                "write dataset " + name);
    }
    H5Dclose(ds);
    H5Sclose(space);
}

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
                          const OpticalEfficiencyConfig& efficiency_cfg,
                          const NsbConfig& nsb_cfg,
                          const TriggerConfig& trigger_cfg,
                          const std::map<SummaryKey, TraceSummary>& summaries,
                          const std::map<PixelKey, PixelAccumulator>& pixels,
                          const std::map<WaveformKey, WaveformPixelAccumulator>& waveforms,
                          const std::vector<RawWaveformHit>& raw_waveform_hits,
                          const std::vector<WhiteboardHdf5Row>& whiteboard_hits)
{
    const bool write_sparse =
        output_cfg.hdf5_storage == "sparse" || output_cfg.hdf5_storage == "both";
    const bool write_dense =
        output_cfg.hdf5_storage == "dense" || output_cfg.hdf5_storage == "both";

    const std::filesystem::path out_path(output_cfg.hdf5_path);
    if (out_path.has_parent_path()) {
        std::filesystem::create_directories(out_path.parent_path());
    }

    hid_t file = H5Fcreate(output_cfg.hdf5_path.c_str(), H5F_ACC_TRUNC,
                           H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
        throw std::runtime_error("failed to create HDF5 file: " + output_cfg.hdf5_path);
    }

    try {
        writeStringAttribute(file, "format", "LACT_sim trace HDF5");
        writeStringAttribute(file, "format_version", "0.1-cpp");
        writeStringAttribute(
            file, "producer_version",
            getString(cfg, "provenance.producer_version", "source-tree"));
        writeStringAttribute(
            file, "source_path",
            getString(cfg, "provenance.source_path", ""));
        writeStringAttribute(
            file, "source_sha256",
            getString(cfg, "provenance.source_sha256", ""));
        writeStringAttribute(file, "image_storage", output_cfg.hdf5_storage);
        writeStringAttribute(file, "hdf5_write_components",
                             output_cfg.hdf5_write_components ? "true" : "false");
        writeStringAttribute(file, "lact_root_write_components",
                             output_cfg.lact_root_write_components ? "true" : "false");
        writeStringAttribute(file, "save_only_triggered",
                             output_cfg.save_only_triggered ? "true" : "false");
        writeStringAttribute(file, "write_pixel_time_stats",
                             output_cfg.write_pixel_time_stats ? "true" : "false");
        writeStringAttribute(file, "waveform_enabled",
                             waveform_cfg.enabled ? "true" : "false");
        writeStringAttribute(file, "hdf5_write_waveforms",
                             output_cfg.hdf5_write_waveforms ? "true" : "false");
        writeStringAttribute(file, "hdf5_waveform_storage",
                             output_cfg.hdf5_waveform_storage);
        writeStringAttribute(file, "event_id_mode", source_runtime_cfg.event_id_mode);
        writeStringAttribute(file, "source_eventio_path", source_runtime_cfg.eventio_path);

        hid_t config_group = H5Gcreate2(file, "config", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        writeStringAttribute(config_group, "main_config_path", main_config_path);
        for (const auto& kv : cfg) {
            writeStringAttribute(config_group, kv.first, kv.second);
        }
        const std::string main_text = readTextIfExists(main_config_path);
        if (!main_text.empty()) {
            writeStringDataset(config_group, "main_config_text", main_text);
        }
        hid_t components = H5Gcreate2(config_group, "components",
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        const std::vector<std::pair<std::string, std::string>> component_items = {
            {"telescope", component_paths.telescope},
            {"mirror", component_paths.mirror},
            {"source", component_paths.source},
            {"output", component_paths.output},
            {"camera", component_paths.camera},
            {"sipm", component_paths.sipm},
            {"electronics", component_paths.electronics},
            {"efficiency", component_paths.efficiency},
            {"atmosphere", component_paths.atmosphere},
            {"error", component_paths.error},
            {"obstruction", component_paths.obstruction},
        };
        for (const auto& item : component_items) {
            if (item.second.empty()) {
                continue;
            }
            writeStringAttribute(components, item.first, item.second);
            const std::string text = readTextIfExists(item.second);
            if (!text.empty()) {
                writeStringDataset(components, item.first + "_text", text);
            }
        }
        H5Gclose(components);
        H5Gclose(config_group);

        hid_t metadata_group = H5Gcreate2(file, "metadata",
                                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t coordinates_group = H5Gcreate2(metadata_group, "coordinates",
                                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        writeStringAttribute(coordinates_group, "array_position_frame",
                             "CORSIKA IACT NWU horizontal frame");
        writeStringAttribute(coordinates_group, "array_x_m",
                             "CORSIKA magnetic-North-positive telescope position coordinate");
        writeStringAttribute(coordinates_group, "array_y_m",
                             "West-positive telescope position coordinate");
        writeStringAttribute(coordinates_group, "array_z_m",
                             "Up-positive telescope position coordinate");
        writeStringAttribute(coordinates_group, "pointing_az_deg",
                             "CORSIKA magnetic-North-to-East azimuth; 0=+array_x, 90=East/-array_y");
        writeStringAttribute(coordinates_group, "pointing_el_deg",
                             "Elevation above local horizon; zenith angle = 90 - elevation");
        writeStringAttribute(coordinates_group, "source_coordinate_frame",
                             source_runtime_cfg.coordinate_frame);
        writeStringAttribute(coordinates_group, "eventio_photon_frame",
                             source_runtime_cfg.coordinate_frame);
        writeStringAttribute(coordinates_group, "eventio_corsika_iact_positions",
                             "Photon bunch x/y/z are telescope-relative CORSIKA IACT coordinates before rotation to telescope-local optics");
        writeStringAttribute(coordinates_group, "eventio_teloff_core_note",
                             "MC_TELOFF is detector-array offset with respect to shower core; /events/corsika stores core = -MC_TELOFF in NWU coordinates");
        writeStringAttribute(coordinates_group, "eventio_telpos_note",
                             "Telescope positions are hessio MC_TELPOS detector coordinates; array_z_up_m may include the detector sphere/radius convention used by the producer");
        writeStringAttribute(coordinates_group, "camera_plane_coordinates",
                             "camera x/y are output-plane u/v coordinates; for the LACT focal plane u=-mirror-local x (East at north pointing) and v=+mirror-local y (sky-up)");
        H5Gclose(coordinates_group);

        hid_t sipm_group = H5Gcreate2(metadata_group, "sipm",
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        writeStringAttribute(sipm_group, "size_m", doubleToString(sipm_cfg.size_m));
        writeStringAttribute(sipm_group, "pde", factorDescription(efficiency_cfg.sipm_pde));
        H5Gclose(sipm_group);

        hid_t efficiency_group = H5Gcreate2(metadata_group, "efficiency",
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        writeStringAttribute(efficiency_group, "constant_scale",
                             doubleToString(efficiency_cfg.constant_scale));
        writeStringAttribute(efficiency_group, "mirror_reflectivity",
                             factorDescription(efficiency_cfg.mirror_reflectivity));
        writeStringAttribute(efficiency_group, "filter_transmission",
                             factorDescription(efficiency_cfg.filter_transmission));
        writeStringAttribute(efficiency_group, "atmosphere",
                             factorDescription(efficiency_cfg.atmosphere_transmission));
        writeStringAttribute(efficiency_group, "funnel_acceptance",
                             efficiency_cfg.use_funnel_acceptance ? "cos(theta)" : "not set -> 1");
        H5Gclose(efficiency_group);

        hid_t electronics_group = H5Gcreate2(metadata_group, "electronics",
                                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        writeStringAttribute(electronics_group, "model", "integrated_pe_placeholder");
        writeStringAttribute(electronics_group, "response",
                             "reserved; SiPM PDE is handled by sipm.pde");
        H5Gclose(electronics_group);

        hid_t waveform_group_meta = H5Gcreate2(metadata_group, "waveform",
                                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        writeStringAttribute(waveform_group_meta, "enabled",
                             waveform_cfg.enabled ? "true" : "false");
        writeStringAttribute(waveform_group_meta, "source", waveform_cfg.source);
        writeStringAttribute(waveform_group_meta, "time_reference",
                             waveform_cfg.time_reference);
        writeStringAttribute(waveform_group_meta, "time_bin_width_ns",
                             doubleToString(waveform_cfg.time_bin_width_ns));
        writeStringAttribute(waveform_group_meta, "time_window_start_ns",
                             doubleToString(waveform_cfg.time_window_start_ns));
        writeStringAttribute(waveform_group_meta, "time_window_end_ns",
                             doubleToString(waveform_cfg.time_window_end_ns));
        writeStringAttribute(waveform_group_meta, "note",
                             "proxy waveform before real electronics; when NSB is enabled it is sampled per time bin");
        H5Gclose(waveform_group_meta);

        hid_t nsb_group = H5Gcreate2(metadata_group, "nsb",
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        writeStringAttribute(nsb_group, "enabled", nsb_cfg.enabled ? "true" : "false");
        writeStringAttribute(nsb_group, "model", nsb_cfg.model);
        writeStringAttribute(nsb_group, "rate_pe_per_ns_per_pixel",
                             doubleToString(nsb_cfg.rate_pe_per_ns_per_pixel));
        writeStringAttribute(nsb_group, "window_ns", doubleToString(nsb_cfg.window_ns));
        writeStringAttribute(nsb_group, "seed", intToString(nsb_cfg.seed));
        writeStringAttribute(nsb_group, "spectrum_csv", nsb_cfg.spectrum_csv);
        writeStringAttribute(nsb_group, "spectrum_unit", nsb_cfg.spectrum_unit);
        writeStringAttribute(nsb_group, "effective_area_m2",
                             doubleToString(nsb_cfg.effective_area_m2));
        writeStringAttribute(nsb_group, "pixel_solid_angle_sr",
                             doubleToString(nsb_cfg.pixel_solid_angle_sr));
        writeStringAttribute(nsb_group, "computed_from_spectrum",
                             nsb_cfg.computed_from_spectrum ? "true" : "false");
        writeStringAttribute(nsb_group, "spectral_integral_pe_s_sr_m2",
                             doubleToString(nsb_cfg.spectral_integral_pe_s_sr_m2));
        H5Gclose(nsb_group);

        hid_t trigger_group_meta = H5Gcreate2(metadata_group, "trigger",
                                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        writeStringAttribute(trigger_group_meta, "enabled",
                             trigger_cfg.enabled ? "true" : "false");
        writeStringAttribute(trigger_group_meta, "model", "simple_multiplicity");
        writeStringAttribute(trigger_group_meta, "pixel_threshold_pe",
                             doubleToString(trigger_cfg.pixel_threshold_pe));
        writeStringAttribute(trigger_group_meta, "camera_multiplicity",
                             intToString(trigger_cfg.camera_multiplicity));
        writeStringAttribute(trigger_group_meta, "array_multiplicity",
                             intToString(trigger_cfg.array_multiplicity));
        writeStringAttribute(trigger_group_meta, "coincidence_window_ns",
                             doubleToString(trigger_cfg.coincidence_window_ns));
        writeStringAttribute(trigger_group_meta, "camera_coincidence_window_ns",
                             doubleToString(trigger_cfg.camera_coincidence_window_ns));
        writeStringAttribute(trigger_group_meta, "array_coincidence_window_ns",
                             doubleToString(trigger_cfg.array_coincidence_window_ns));
        writeStringAttribute(trigger_group_meta, "array_time_correction",
                             trigger_cfg.array_time_correction);
        writeStringAttribute(trigger_group_meta, "array_wavefront_speed_m_per_ns",
                             doubleToString(
                                 trigger_cfg.array_time_correction == "plane_wave"
                                     ? resolveEventIOArrayWavefrontSpeedMPerNs(
                                           trigger_cfg, metadata)
                                     : trigger_cfg.array_wavefront_speed_m_per_ns,
                                 15));
        H5Gclose(trigger_group_meta);
        H5Gclose(metadata_group);

        struct CameraRow {
            std::int32_t pixel_id;
            float x_m;
            float y_m;
            float size_m;
            std::int16_t shape_code;
        };
        std::vector<CameraRow> camera_rows;
        camera_rows.reserve(camera.size());
        for (const auto& pixel : camera.pixels()) {
            camera_rows.push_back(CameraRow{
                static_cast<std::int32_t>(pixel.id),
                static_cast<float>(pixel.center.x),
                static_cast<float>(pixel.center.y),
                static_cast<float>(pixel.size),
                static_cast<std::int16_t>(hdf5PixelShapeCode(pixel.shape)),
            });
        }
        hid_t camera_group = H5Gcreate2(file, "camera", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        writeStringAttribute(camera_group, "shape_code_map",
                             "0=unknown,1=square,2=hexagon,3=circular");
        hid_t camera_type = H5Tcreate(H5T_COMPOUND, sizeof(CameraRow));
        H5Tinsert(camera_type, "pixel_id", HOFFSET(CameraRow, pixel_id), H5T_NATIVE_INT32);
        H5Tinsert(camera_type, "x_m", HOFFSET(CameraRow, x_m), H5T_NATIVE_FLOAT);
        H5Tinsert(camera_type, "y_m", HOFFSET(CameraRow, y_m), H5T_NATIVE_FLOAT);
        H5Tinsert(camera_type, "size_m", HOFFSET(CameraRow, size_m), H5T_NATIVE_FLOAT);
        H5Tinsert(camera_type, "shape_code", HOFFSET(CameraRow, shape_code), H5T_NATIVE_INT16);
        writeCompound1D(camera_group, "pixels", camera_type, camera_rows);
        H5Tclose(camera_type);
        H5Gclose(camera_group);

        struct FacetRow {
            std::int32_t mirror_id;
            float center_x_m;
            float center_y_m;
            float center_z_m;
            float normal_x;
            float normal_y;
            float normal_z;
            float radius_of_curvature_m;
            float size1_m;
            float size2_m;
            float aperture_rotation_rad;
            std::int16_t shape_code;
        };
        std::vector<FacetRow> facet_rows;
        facet_rows.reserve(facets.size());
        for (const auto& facet : facets) {
            facet_rows.push_back(FacetRow{
                static_cast<std::int32_t>(facet.id),
                static_cast<float>(facet.center.x),
                static_cast<float>(facet.center.y),
                static_cast<float>(facet.center.z),
                static_cast<float>(facet.normal.x),
                static_cast<float>(facet.normal.y),
                static_cast<float>(facet.normal.z),
                static_cast<float>(facet.radius_of_curvature),
                static_cast<float>(facet.size1),
                static_cast<float>(facet.size2),
                static_cast<float>(facet.aperture_rotation_rad),
                static_cast<std::int16_t>(hdf5FacetShapeCode(facet.aperture_shape)),
            });
        }
        hid_t mirror_group = H5Gcreate2(file, "mirrors", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        writeStringAttribute(mirror_group, "shape_code_map",
                             "0=unknown,1=square,2=hexagon,3=circular");
        hid_t facet_type = H5Tcreate(H5T_COMPOUND, sizeof(FacetRow));
        H5Tinsert(facet_type, "mirror_id", HOFFSET(FacetRow, mirror_id), H5T_NATIVE_INT32);
        H5Tinsert(facet_type, "center_x_m", HOFFSET(FacetRow, center_x_m), H5T_NATIVE_FLOAT);
        H5Tinsert(facet_type, "center_y_m", HOFFSET(FacetRow, center_y_m), H5T_NATIVE_FLOAT);
        H5Tinsert(facet_type, "center_z_m", HOFFSET(FacetRow, center_z_m), H5T_NATIVE_FLOAT);
        H5Tinsert(facet_type, "normal_x", HOFFSET(FacetRow, normal_x), H5T_NATIVE_FLOAT);
        H5Tinsert(facet_type, "normal_y", HOFFSET(FacetRow, normal_y), H5T_NATIVE_FLOAT);
        H5Tinsert(facet_type, "normal_z", HOFFSET(FacetRow, normal_z), H5T_NATIVE_FLOAT);
        H5Tinsert(facet_type, "radius_of_curvature_m",
                  HOFFSET(FacetRow, radius_of_curvature_m), H5T_NATIVE_FLOAT);
        H5Tinsert(facet_type, "size1_m", HOFFSET(FacetRow, size1_m), H5T_NATIVE_FLOAT);
        H5Tinsert(facet_type, "size2_m", HOFFSET(FacetRow, size2_m), H5T_NATIVE_FLOAT);
        H5Tinsert(facet_type, "aperture_rotation_rad",
                  HOFFSET(FacetRow, aperture_rotation_rad), H5T_NATIVE_FLOAT);
        H5Tinsert(facet_type, "shape_code", HOFFSET(FacetRow, shape_code), H5T_NATIVE_INT16);
        writeCompound1D(mirror_group, "facets", facet_type, facet_rows);
        H5Tclose(facet_type);
        H5Gclose(mirror_group);

        struct TelescopeRow {
            std::int32_t telescope_id;
            double x_m;
            double y_m;
            double z_m;
            double array_x_north_m;
            double array_y_west_m;
            double array_z_up_m;
            double radius_m;
            double pointing_az_deg;
            double pointing_el_deg;
            double focal_length_m;
        };
        std::vector<TelescopeRow> telescope_rows;
        if (!metadata.telescopes.empty()) {
            telescope_rows.reserve(metadata.telescopes.size());
            for (const auto& tel : metadata.telescopes) {
                telescope_rows.push_back(TelescopeRow{
                    static_cast<std::int32_t>(tel.telescope_id),
                    tel.x_m,
                    tel.y_m,
                    tel.z_m,
                    tel.x_m,
                    tel.y_m,
                    tel.z_m,
                    tel.radius_m,
                    telescope_cfg.pointing_az_deg,
                    telescope_cfg.pointing_el_deg,
                    telescope_cfg.focal_length_m,
                });
            }
        } else {
            telescope_rows.push_back(TelescopeRow{
                static_cast<std::int32_t>(telescope_cfg.id),
                telescope_cfg.position_m.x,
                telescope_cfg.position_m.y,
                telescope_cfg.position_m.z,
                telescope_cfg.position_m.x,
                telescope_cfg.position_m.y,
                telescope_cfg.position_m.z,
                0.0,
                telescope_cfg.pointing_az_deg,
                telescope_cfg.pointing_el_deg,
                telescope_cfg.focal_length_m,
            });
        }
        hid_t telescope_group = H5Gcreate2(file, "telescopes",
                                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        writeStringAttribute(telescope_group, "coordinate_frame",
                             "CORSIKA IACT horizontal frame");
        writeStringAttribute(telescope_group, "x_m_compat",
                             "same as array_x_north_m; kept for compatibility");
        writeStringAttribute(telescope_group, "y_m_compat",
                             "same as array_y_west_m; kept for compatibility");
        writeStringAttribute(telescope_group, "z_m_compat",
                             "same as array_z_up_m; kept for compatibility");
        writeStringAttribute(telescope_group, "pointing_convention",
                             "azimuth North-to-East, elevation above horizon");
        hid_t telescope_type = H5Tcreate(H5T_COMPOUND, sizeof(TelescopeRow));
        H5Tinsert(telescope_type, "telescope_id",
                  HOFFSET(TelescopeRow, telescope_id), H5T_NATIVE_INT32);
        H5Tinsert(telescope_type, "x_m", HOFFSET(TelescopeRow, x_m), H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_type, "y_m", HOFFSET(TelescopeRow, y_m), H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_type, "z_m", HOFFSET(TelescopeRow, z_m), H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_type, "array_x_north_m",
                  HOFFSET(TelescopeRow, array_x_north_m), H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_type, "array_y_west_m",
                  HOFFSET(TelescopeRow, array_y_west_m), H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_type, "array_z_up_m",
                  HOFFSET(TelescopeRow, array_z_up_m), H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_type, "radius_m",
                  HOFFSET(TelescopeRow, radius_m), H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_type, "pointing_az_deg",
                  HOFFSET(TelescopeRow, pointing_az_deg), H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_type, "pointing_el_deg",
                  HOFFSET(TelescopeRow, pointing_el_deg), H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_type, "focal_length_m",
                  HOFFSET(TelescopeRow, focal_length_m), H5T_NATIVE_DOUBLE);
        writeCompound1D(telescope_group, "table", telescope_type, telescope_rows);
        H5Tclose(telescope_type);
        H5Gclose(telescope_group);

        std::set<SummaryKey> image_keys;
        for (const auto& kv : summaries) {
            image_keys.insert(kv.first);
        }
        for (const auto& kv : pixels) {
            image_keys.insert({std::get<0>(kv.first), std::get<1>(kv.first)});
        }

        struct SparsePixelRow {
            std::int32_t pixel_id;
            std::int32_t photon_count;
            float pe;
            float signal;
            float time_mean_ns;
            float time_rms_ns;
        };
        struct ImageIndexRow {
            std::int32_t image_index;
            std::int64_t event_id;
            std::int32_t telescope_id;
            std::int64_t start;
            std::int32_t count;
            double total_photons;
            double total_pe;
            double total_signal;
            double time_mean_ns;
            double time_rms_ns;
            double time_first_ns;
        };
        std::vector<SparsePixelRow> sparse_rows;
        std::vector<ImageIndexRow> image_rows;
        image_rows.reserve(image_keys.size());

        std::int32_t image_index = 0;
        for (const auto& key : image_keys) {
            const auto start = static_cast<std::int64_t>(sparse_rows.size());
            const int event_id = key.first;
            const int telescope_id = key.second;
            double total_signal = 0.0;
            double total_pe = 0.0;
            double total_photons = 0.0;
            const PixelKey pixel_begin{
                event_id, telescope_id, std::numeric_limits<int>::min()};
            const PixelKey pixel_end{
                event_id, telescope_id, std::numeric_limits<int>::max()};
            for (auto it = pixels.lower_bound(pixel_begin);
                 it != pixels.end() && it->first <= pixel_end;
                 ++it) {
                const auto& kv = *it;
                const auto& p = kv.second;
                const double mean = p.signal > 0.0 ? p.time_sum / p.signal : 0.0;
                const double var = p.signal > 0.0
                    ? std::max(0.0, p.time2_sum / p.signal - mean * mean)
                    : 0.0;
                sparse_rows.push_back(SparsePixelRow{
                    static_cast<std::int32_t>(p.pixel_id),
                    static_cast<std::int32_t>(p.photon_count),
                    static_cast<float>(p.pe),
                    static_cast<float>(p.signal),
                    static_cast<float>(mean),
                    static_cast<float>(std::sqrt(var)),
                });
                total_signal += p.signal;
                total_pe += p.pe;
                total_photons += static_cast<double>(p.photon_count);
            }

            double mean = 0.0;
            double rms = 0.0;
            double first = 0.0;
            auto summary_it = summaries.find(key);
            if (summary_it != summaries.end()) {
                const auto& s = summary_it->second;
                total_signal = s.weighted_signal;
                total_pe = s.weighted_signal;
                total_photons = static_cast<double>(s.hit_camera);
                if (std::isfinite(s.first_cherenkov_time_ns)) {
                    first = s.first_cherenkov_time_ns;
                }
                if (s.weighted_signal > 0.0) {
                    const double m = s.weighted_time_sum / s.weighted_signal;
                    const double v = std::max(0.0, s.weighted_time2_sum / s.weighted_signal - m * m);
                    mean = m;
                    rms = std::sqrt(v);
                }
            }
            const auto count = static_cast<std::int32_t>(
                static_cast<std::int64_t>(sparse_rows.size()) - start);
            image_rows.push_back(ImageIndexRow{
                image_index++,
                static_cast<std::int64_t>(event_id),
                static_cast<std::int32_t>(telescope_id),
                start,
                count,
                total_photons,
                total_pe,
                total_signal,
                mean,
                rms,
                first,
            });
        }

        std::vector<std::int32_t> pixel_id_axis;
        std::map<std::int32_t, std::size_t> pixel_to_col;
        std::vector<float> dense_signal;
        std::vector<float> dense_pe;
        std::vector<float> dense_cherenkov_pe;
        std::vector<float> dense_nsb_pe;
        std::vector<float> dense_time_mean_ns;
        std::vector<float> dense_time_rms_ns;
        std::vector<std::int32_t> dense_photon_count;
        const bool have_dense_images = write_dense && !camera_rows.empty();
        const bool have_camera_axis = !camera_rows.empty();
        if (have_camera_axis) {
            pixel_id_axis.reserve(camera_rows.size());
            for (std::size_t i = 0; i < camera_rows.size(); ++i) {
                pixel_id_axis.push_back(camera_rows[i].pixel_id);
                pixel_to_col[camera_rows[i].pixel_id] = i;
            }
        }

        if (have_dense_images) {
            const std::size_t n_images = image_rows.size();
            const std::size_t n_pixels = camera_rows.size();
            dense_signal.assign(n_images * n_pixels, 0.0f);
            dense_pe.assign(n_images * n_pixels, 0.0f);
            dense_photon_count.assign(n_images * n_pixels, 0);
            if (output_cfg.write_pixel_time_stats) {
                dense_time_mean_ns.assign(n_images * n_pixels, 0.0f);
                dense_time_rms_ns.assign(n_images * n_pixels, 0.0f);
            }
            for (const auto& image : image_rows) {
                const std::size_t row = static_cast<std::size_t>(image.image_index);
                const std::int64_t begin = image.start;
                const std::int64_t end = image.start + image.count;
                for (std::int64_t i = begin; i < end; ++i) {
                    const auto& pixel = sparse_rows[static_cast<std::size_t>(i)];
                    const auto col_it = pixel_to_col.find(pixel.pixel_id);
                    if (col_it == pixel_to_col.end()) {
                        continue;
                    }
                    const std::size_t index = row * n_pixels + col_it->second;
                    dense_signal[index] = pixel.signal;
                    dense_pe[index] = pixel.pe;
                    dense_photon_count[index] = pixel.photon_count;
                    if (output_cfg.write_pixel_time_stats) {
                        dense_time_mean_ns[index] = pixel.time_mean_ns;
                        dense_time_rms_ns[index] = pixel.time_rms_ns;
                    }
                }
            }

            dense_cherenkov_pe = dense_pe;
            dense_nsb_pe.assign(n_images * n_pixels, 0.0f);
            if (nsb_cfg.enabled && nsb_cfg.rate_pe_per_ns_per_pixel > 0.0 &&
                nsb_cfg.window_ns > 0.0) {
                if (waveform_cfg.enabled && waveform_cfg.source == "pe") {
                    const std::size_t n_bins = waveformBinCount(waveform_cfg);
                    for (const auto& image : image_rows) {
                        const std::size_t row =
                            static_cast<std::size_t>(image.image_index);
                        generateTimeBinnedNsbPe(
                            nsb_cfg,
                            waveform_cfg,
                            static_cast<int>(image.event_id),
                            static_cast<int>(image.telescope_id),
                            n_pixels,
                            n_bins,
                            [&](std::size_t col, std::size_t, float nsb_pe) {
                                const std::size_t index = row * n_pixels + col;
                                dense_nsb_pe[index] += nsb_pe;
                                dense_pe[index] += nsb_pe;
                                dense_signal[index] += nsb_pe;
                            });
                    }
                } else {
                    for (const auto& image : image_rows) {
                        const std::size_t row =
                            static_cast<std::size_t>(image.image_index);
                        generateIntegratedNsbPe(
                            nsb_cfg,
                            static_cast<int>(image.event_id),
                            static_cast<int>(image.telescope_id),
                            n_pixels,
                            nsb_cfg.window_ns,
                            [&](std::size_t col, float nsb_pe) {
                                const std::size_t i = row * n_pixels + col;
                                dense_nsb_pe[i] = nsb_pe;
                                dense_pe[i] += nsb_pe;
                                dense_signal[i] += nsb_pe;
                            });
                    }
                }
            }

            for (auto& image : image_rows) {
                const std::size_t row = static_cast<std::size_t>(image.image_index);
                double total_pe = 0.0;
                double total_signal = 0.0;
                for (std::size_t col = 0; col < n_pixels; ++col) {
                    const std::size_t index = row * n_pixels + col;
                    total_pe += dense_pe[index];
                    total_signal += dense_signal[index];
                }
                image.total_pe = total_pe;
                image.total_signal = total_signal;
            }
        }

        struct TelescopeTriggerRow {
            std::int64_t event_id;
            std::int32_t telescope_id;
            std::int8_t triggered;
            std::int32_t n_pixels_above_threshold;
            double total_pe;
            double trigger_time_ns;
            double trigger_first_time_ns;
            double trigger_max_multiplicity_time_ns;
            double geometric_delay_ns = std::numeric_limits<double>::quiet_NaN();
            double coincidence_time_ns = std::numeric_limits<double>::quiet_NaN();
        };
        struct ArrayTriggerRow {
            std::int64_t event_id;
            std::int8_t array_triggered;
            std::int32_t n_triggered_telescopes;
        };
        std::vector<TelescopeTriggerRow> telescope_trigger_rows;
        std::map<SummaryKey, std::size_t> telescope_trigger_row_index;
        std::map<int, std::vector<TelescopeTriggerTime>>
            telescope_trigger_times_by_event;
        telescope_trigger_rows.reserve(image_rows.size());
        for (const auto& image : image_rows) {
            double total_pe = image.total_pe;
            std::size_t trigger_pixels = 0;
            std::size_t trigger_bins = 1;
            double trigger_bin_width_ns = 1.0;
            double first_trigger_bin_center_ns =
                std::isfinite(image.time_mean_ns) ? image.time_mean_ns : 0.0;
            std::function<double(std::size_t, std::size_t)> pe_at;
            std::unordered_map<std::size_t, double> trigger_waveform_pe;
            std::vector<double> sparse_trigger_pe;
            const bool use_time_trigger =
                waveform_cfg.enabled && waveform_cfg.source == "pe" &&
                !camera_rows.empty();
            if (use_time_trigger) {
                trigger_pixels = camera_rows.size();
                trigger_bins = waveformBinCount(waveform_cfg);
                trigger_bin_width_ns = waveform_cfg.time_bin_width_ns;
                double reference_time_ns = 0.0;
                if (waveform_cfg.time_reference == "image_first" &&
                    std::isfinite(image.time_first_ns)) {
                    reference_time_ns = image.time_first_ns;
                } else if (waveform_cfg.time_reference == "image_mean" &&
                           std::isfinite(image.time_mean_ns)) {
                    reference_time_ns = image.time_mean_ns;
                }
                first_trigger_bin_center_ns =
                    reference_time_ns + waveform_cfg.time_window_start_ns +
                    0.5 * waveform_cfg.time_bin_width_ns;
                auto add_trigger_pe = [&](int pixel_id, int bin, double pe) {
                    const auto col_it = pixel_to_col.find(pixel_id);
                    if (col_it == pixel_to_col.end() || bin < 0 ||
                        static_cast<std::size_t>(bin) >= trigger_bins || pe == 0.0) {
                        return;
                    }
                    trigger_waveform_pe[
                        col_it->second * trigger_bins + static_cast<std::size_t>(bin)] += pe;
                };
                if (waveformUsesImageReference(waveform_cfg)) {
                    for (const auto& hit : raw_waveform_hits) {
                        if (hit.event_id != image.event_id ||
                            hit.telescope_id != image.telescope_id) {
                            continue;
                        }
                        add_trigger_pe(
                            hit.pixel_id,
                            waveformBinForTime(
                                waveform_cfg, hit.time_ns - reference_time_ns),
                            hit.pe);
                    }
                } else {
                    const WaveformKey begin_key{
                        static_cast<int>(image.event_id),
                        static_cast<int>(image.telescope_id),
                        std::numeric_limits<int>::min(),
                        std::numeric_limits<int>::min()};
                    const WaveformKey end_key{
                        static_cast<int>(image.event_id),
                        static_cast<int>(image.telescope_id),
                        std::numeric_limits<int>::max(),
                        std::numeric_limits<int>::max()};
                    for (auto it = waveforms.lower_bound(begin_key);
                         it != waveforms.end() && it->first <= end_key;
                         ++it) {
                        add_trigger_pe(
                            it->second.pixel_id, it->second.time_bin, it->second.pe);
                    }
                }

                // Use the same deterministic NSB realization as the dense
                // image and serialized waveform.  A separately seeded
                // per-cell sampler would have the same distribution but
                // would make the trigger impossible to reproduce from the
                // saved waveform.
                generateTimeBinnedNsbPe(
                    nsb_cfg,
                    waveform_cfg,
                    static_cast<int>(image.event_id),
                    static_cast<int>(image.telescope_id),
                    trigger_pixels,
                    trigger_bins,
                    [&](std::size_t col, std::size_t bin, float pe) {
                        trigger_waveform_pe[col * trigger_bins + bin] += pe;
                    });
                pe_at = [&](std::size_t col, std::size_t bin) {
                    const auto found = trigger_waveform_pe.find(
                        col * trigger_bins + bin);
                    return found == trigger_waveform_pe.end()
                        ? 0.0
                        : found->second;
                };
            } else if (have_dense_images) {
                const std::size_t row = static_cast<std::size_t>(image.image_index);
                const std::size_t n_pixels = camera_rows.size();
                total_pe = 0.0;
                for (std::size_t col = 0; col < n_pixels; ++col) {
                    total_pe += dense_pe[row * n_pixels + col];
                }
                trigger_pixels = n_pixels;
                pe_at = [&, row, n_pixels](std::size_t col, std::size_t) {
                    return static_cast<double>(dense_pe[row * n_pixels + col]);
                };
            } else {
                const std::int64_t begin = image.start;
                const std::int64_t end = image.start + image.count;
                trigger_pixels = camera_rows.size();
                sparse_trigger_pe.assign(trigger_pixels, 0.0);
                for (std::int64_t i = begin; i < end; ++i) {
                    const auto& pixel = sparse_rows[static_cast<std::size_t>(i)];
                    const auto col_it = pixel_to_col.find(pixel.pixel_id);
                    if (col_it != pixel_to_col.end()) {
                        sparse_trigger_pe[col_it->second] += pixel.pe;
                    }
                }
                pe_at = [&](std::size_t col, std::size_t) {
                    return sparse_trigger_pe[col];
                };
            }
            const auto camera_trigger = evaluateBinnedPeTrigger(
                trigger_pixels,
                trigger_bins,
                trigger_bin_width_ns,
                first_trigger_bin_center_ns,
                trigger_cfg, pe_at);
            if (camera_trigger.triggered) {
                telescope_trigger_times_by_event[static_cast<int>(image.event_id)]
                    .push_back(TelescopeTriggerTime{
                        static_cast<int>(image.telescope_id),
                        camera_trigger.trigger_time_ns});
            }
            telescope_trigger_row_index[{
                static_cast<int>(image.event_id),
                static_cast<int>(image.telescope_id)}] =
                telescope_trigger_rows.size();
            telescope_trigger_rows.push_back(TelescopeTriggerRow{
                image.event_id,
                image.telescope_id,
                static_cast<std::int8_t>(camera_trigger.triggered ? 1 : 0),
                static_cast<std::int32_t>(
                    camera_trigger.n_pixels_above_threshold),
                total_pe,
                camera_trigger.trigger_time_ns,
                camera_trigger.first_trigger_time_ns,
                camera_trigger.max_multiplicity_time_ns,
            });
        }
        std::set<int> trigger_event_ids;
        for (const auto& key : image_keys) {
            trigger_event_ids.insert(key.first);
        }
        std::vector<ArrayTriggerRow> array_trigger_rows;
        std::map<int, ArrayTriggerDecision> array_trigger_decisions;
        array_trigger_rows.reserve(trigger_event_ids.size());
        for (const int event_id : trigger_event_ids) {
            applyEventIOArrayTimingCorrection(
                telescope_trigger_times_by_event[event_id], event_id,
                source_runtime_cfg.event_id_mode, trigger_cfg,
                telescope_cfg, metadata);
            for (const auto& trigger_time :
                 telescope_trigger_times_by_event[event_id]) {
                const auto row_index = telescope_trigger_row_index.find({
                    event_id, trigger_time.telescope_id});
                if (row_index == telescope_trigger_row_index.end()) {
                    throw std::runtime_error(
                        "array timing result has no HDF5 telescope trigger row");
                }
                auto& row = telescope_trigger_rows[row_index->second];
                row.geometric_delay_ns = trigger_time.geometric_delay_ns;
                row.coincidence_time_ns =
                    std::isfinite(trigger_time.coincidence_time_ns)
                        ? trigger_time.coincidence_time_ns
                        : trigger_time.trigger_time_ns;
            }
            const auto decision = evaluateArrayTrigger(
                telescope_trigger_times_by_event[event_id], trigger_cfg);
            array_trigger_decisions[event_id] = decision;
            const int n_triggered = static_cast<int>(
                telescope_trigger_times_by_event[event_id].size());
            array_trigger_rows.push_back(ArrayTriggerRow{
                static_cast<std::int64_t>(event_id),
                static_cast<std::int8_t>(decision.triggered ? 1 : 0),
                static_cast<std::int32_t>(n_triggered),
            });
        }

        if (output_cfg.save_only_triggered && trigger_cfg.enabled) {
            std::set<SummaryKey> triggered_image_keys;
            for (const auto& row : telescope_trigger_rows) {
                const auto& array_decision =
                    array_trigger_decisions[static_cast<int>(row.event_id)];
                const bool telescope_is_coincident = std::binary_search(
                    array_decision.coincident_telescope_ids.begin(),
                    array_decision.coincident_telescope_ids.end(),
                    static_cast<int>(row.telescope_id));
                if (row.triggered && array_decision.triggered &&
                    telescope_is_coincident) {
                    triggered_image_keys.insert({
                        static_cast<int>(row.event_id),
                        static_cast<int>(row.telescope_id),
                    });
                }
            }

            std::vector<ImageIndexRow> filtered_image_rows;
            std::vector<SparsePixelRow> filtered_sparse_rows;
            std::vector<float> filtered_dense_signal;
            std::vector<float> filtered_dense_pe;
            std::vector<float> filtered_dense_cherenkov_pe;
            std::vector<float> filtered_dense_nsb_pe;
            std::vector<float> filtered_dense_time_mean_ns;
            std::vector<float> filtered_dense_time_rms_ns;
            std::vector<std::int32_t> filtered_dense_photon_count;

            const std::size_t n_pixels = camera_rows.size();
            if (have_dense_images) {
                filtered_dense_signal.reserve(triggered_image_keys.size() * n_pixels);
                filtered_dense_pe.reserve(triggered_image_keys.size() * n_pixels);
                filtered_dense_photon_count.reserve(triggered_image_keys.size() * n_pixels);
                if (output_cfg.hdf5_write_components) {
                    filtered_dense_cherenkov_pe.reserve(triggered_image_keys.size() * n_pixels);
                    filtered_dense_nsb_pe.reserve(triggered_image_keys.size() * n_pixels);
                }
                if (output_cfg.write_pixel_time_stats) {
                    filtered_dense_time_mean_ns.reserve(triggered_image_keys.size() * n_pixels);
                    filtered_dense_time_rms_ns.reserve(triggered_image_keys.size() * n_pixels);
                }
            }

            for (const auto& image : image_rows) {
                const SummaryKey key{
                    static_cast<int>(image.event_id),
                    static_cast<int>(image.telescope_id),
                };
                if (triggered_image_keys.find(key) == triggered_image_keys.end()) {
                    continue;
                }

                ImageIndexRow filtered = image;
                filtered.image_index = static_cast<std::int32_t>(filtered_image_rows.size());
                filtered.start = static_cast<std::int64_t>(filtered_sparse_rows.size());

                const std::int64_t begin = image.start;
                const std::int64_t end = image.start + image.count;
                for (std::int64_t i = begin; i < end; ++i) {
                    filtered_sparse_rows.push_back(sparse_rows[static_cast<std::size_t>(i)]);
                }
                filtered.count = static_cast<std::int32_t>(
                    static_cast<std::int64_t>(filtered_sparse_rows.size()) - filtered.start);

                if (have_dense_images) {
                    const std::size_t old_row = static_cast<std::size_t>(image.image_index);
                    const std::size_t old_begin = old_row * n_pixels;
                    const std::size_t old_end = old_begin + n_pixels;
                    filtered_dense_signal.insert(filtered_dense_signal.end(),
                                                 dense_signal.begin() + old_begin,
                                                 dense_signal.begin() + old_end);
                    filtered_dense_pe.insert(filtered_dense_pe.end(),
                                             dense_pe.begin() + old_begin,
                                             dense_pe.begin() + old_end);
                    filtered_dense_photon_count.insert(filtered_dense_photon_count.end(),
                                                       dense_photon_count.begin() + old_begin,
                                                       dense_photon_count.begin() + old_end);
                    if (output_cfg.hdf5_write_components) {
                        filtered_dense_cherenkov_pe.insert(filtered_dense_cherenkov_pe.end(),
                                                           dense_cherenkov_pe.begin() + old_begin,
                                                           dense_cherenkov_pe.begin() + old_end);
                        filtered_dense_nsb_pe.insert(filtered_dense_nsb_pe.end(),
                                                     dense_nsb_pe.begin() + old_begin,
                                                     dense_nsb_pe.begin() + old_end);
                    }
                    if (output_cfg.write_pixel_time_stats) {
                        filtered_dense_time_mean_ns.insert(filtered_dense_time_mean_ns.end(),
                                                           dense_time_mean_ns.begin() + old_begin,
                                                           dense_time_mean_ns.begin() + old_end);
                        filtered_dense_time_rms_ns.insert(filtered_dense_time_rms_ns.end(),
                                                          dense_time_rms_ns.begin() + old_begin,
                                                          dense_time_rms_ns.begin() + old_end);
                    }
                }

                filtered_image_rows.push_back(filtered);
            }

            image_rows.swap(filtered_image_rows);
            sparse_rows.swap(filtered_sparse_rows);
            if (have_dense_images) {
                dense_signal.swap(filtered_dense_signal);
                dense_pe.swap(filtered_dense_pe);
                dense_photon_count.swap(filtered_dense_photon_count);
                if (output_cfg.hdf5_write_components) {
                    dense_cherenkov_pe.swap(filtered_dense_cherenkov_pe);
                    dense_nsb_pe.swap(filtered_dense_nsb_pe);
                }
                if (output_cfg.write_pixel_time_stats) {
                    dense_time_mean_ns.swap(filtered_dense_time_mean_ns);
                    dense_time_rms_ns.swap(filtered_dense_time_rms_ns);
                }
            }
        }

        struct EventRow {
            std::int32_t event_index;
            std::int64_t event_id;
        };
        struct CorsikaEventRow {
            std::int64_t event_id;
            std::int32_t shower_event_id;
            std::int32_t array_id;
            std::int32_t has_explicit_area_weight;
            std::int32_t primary_type;
            double energy_gev;
            double theta_deg;
            double phi_deg;
            double azimuth_north_to_east_deg;
            double core_x_north_m;
            double core_y_west_m;
            double array_rotation_deg;
            double array_time_offset_ns;
            double area_weight_m2;
            double h_first_int_m;
            double x_max_g_cm2;
            double h_max_m;
            double starting_grammage_g_cm2;
            double ground_gammas;
            double ground_electrons;
            double ground_hadrons;
            double ground_muons;
        };
        struct CorsikaShowerRow {
            std::int32_t shower_event_id;
            std::int32_t primary_type;
            double energy_gev;
            double theta_deg;
            double phi_deg;
            double azimuth_north_to_east_deg;
            double core_x_north_m;
            double core_y_west_m;
            double array_rotation_deg;
            double h_first_int_m;
            double x_max_g_cm2;
            double h_max_m;
            double starting_grammage_g_cm2;
            double ground_gammas;
            double ground_electrons;
            double ground_hadrons;
            double ground_muons;
        };
        std::vector<EventRow> event_rows;
        std::vector<CorsikaEventRow> corsika_event_rows;
        std::set<int> event_ids;
        for (const auto& image : image_rows) {
            event_ids.insert(static_cast<int>(image.event_id));
        }
        int event_index = 0;
        for (const int event_id : event_ids) {
            event_rows.push_back(EventRow{
                event_index++,
                static_cast<std::int64_t>(event_id),
            });
            const OutputEventMetadata event_meta = outputEventMetadata(
                event_id, source_runtime_cfg.event_id_mode, metadata);
            const int shower_event = event_meta.shower_event;
            const int array_id = event_meta.array_id;
            auto event_it = std::find_if(
                metadata.events.begin(), metadata.events.end(),
                [shower_event](const EventIOEventHeader& event) {
                    return event.shower_event_id == shower_event;
                });
            if (event_it != metadata.events.end()) {
                corsika_event_rows.push_back(CorsikaEventRow{
                    static_cast<std::int64_t>(event_id),
                    static_cast<std::int32_t>(shower_event),
                    static_cast<std::int32_t>(array_id),
                    static_cast<std::int32_t>(event_meta.has_explicit_area_weight),
                    static_cast<std::int32_t>(event_it->primary_type),
                    event_it->energy_gev,
                    event_it->theta_deg,
                    event_it->phi_deg,
                    event_it->azimuth_north_to_east_deg,
                    event_meta.core_x_north_m,
                    event_meta.core_y_west_m,
                    event_it->array_rotation_deg,
                    event_meta.array_time_offset_ns,
                    event_meta.area_weight_m2,
                    event_it->h_first_int_m,
                    event_it->x_max_g_cm2,
                    event_it->h_max_m,
                    event_it->starting_grammage_g_cm2,
                    event_it->ground_gammas,
                    event_it->ground_electrons,
                    event_it->ground_hadrons,
                    event_it->ground_muons,
                });
            }
        }
        hid_t events_group = H5Gcreate2(file, "events", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t event_type = H5Tcreate(H5T_COMPOUND, sizeof(EventRow));
        H5Tinsert(event_type, "event_index", HOFFSET(EventRow, event_index), H5T_NATIVE_INT32);
        H5Tinsert(event_type, "event_id", HOFFSET(EventRow, event_id), H5T_NATIVE_INT64);
        writeCompound1D(events_group, "table", event_type, event_rows);
        H5Tclose(event_type);
        if (!corsika_event_rows.empty()) {
            hid_t corsika_event_type =
                H5Tcreate(H5T_COMPOUND, sizeof(CorsikaEventRow));
            H5Tinsert(corsika_event_type, "event_id",
                      HOFFSET(CorsikaEventRow, event_id), H5T_NATIVE_INT64);
            H5Tinsert(corsika_event_type, "shower_event_id",
                      HOFFSET(CorsikaEventRow, shower_event_id), H5T_NATIVE_INT32);
            H5Tinsert(corsika_event_type, "array_id",
                      HOFFSET(CorsikaEventRow, array_id), H5T_NATIVE_INT32);
            H5Tinsert(corsika_event_type, "has_explicit_area_weight",
                      HOFFSET(CorsikaEventRow, has_explicit_area_weight),
                      H5T_NATIVE_INT32);
            H5Tinsert(corsika_event_type, "primary_type",
                      HOFFSET(CorsikaEventRow, primary_type), H5T_NATIVE_INT32);
            H5Tinsert(corsika_event_type, "energy_gev",
                      HOFFSET(CorsikaEventRow, energy_gev), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "theta_deg",
                      HOFFSET(CorsikaEventRow, theta_deg), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "phi_deg",
                      HOFFSET(CorsikaEventRow, phi_deg), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "azimuth_north_to_east_deg",
                      HOFFSET(CorsikaEventRow, azimuth_north_to_east_deg),
                      H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "core_x_north_m",
                      HOFFSET(CorsikaEventRow, core_x_north_m), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "core_y_west_m",
                      HOFFSET(CorsikaEventRow, core_y_west_m), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "array_rotation_deg",
                      HOFFSET(CorsikaEventRow, array_rotation_deg), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "array_time_offset_ns",
                      HOFFSET(CorsikaEventRow, array_time_offset_ns), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "area_weight_m2",
                      HOFFSET(CorsikaEventRow, area_weight_m2), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "h_first_int_m",
                      HOFFSET(CorsikaEventRow, h_first_int_m), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "x_max_g_cm2",
                      HOFFSET(CorsikaEventRow, x_max_g_cm2), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "h_max_m",
                      HOFFSET(CorsikaEventRow, h_max_m), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "starting_grammage_g_cm2",
                      HOFFSET(CorsikaEventRow, starting_grammage_g_cm2), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "ground_gammas",
                      HOFFSET(CorsikaEventRow, ground_gammas), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "ground_electrons",
                      HOFFSET(CorsikaEventRow, ground_electrons), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "ground_hadrons",
                      HOFFSET(CorsikaEventRow, ground_hadrons), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_event_type, "ground_muons",
                      HOFFSET(CorsikaEventRow, ground_muons), H5T_NATIVE_DOUBLE);
            writeCompound1D(events_group, "corsika", corsika_event_type,
                            corsika_event_rows);
            H5Tclose(corsika_event_type);
        }
        if (!metadata.events.empty()) {
            std::vector<CorsikaShowerRow> corsika_shower_rows;
            corsika_shower_rows.reserve(metadata.events.size());
            for (const auto& event : metadata.events) {
                corsika_shower_rows.push_back(CorsikaShowerRow{
                    static_cast<std::int32_t>(event.shower_event_id),
                    static_cast<std::int32_t>(event.primary_type),
                    event.energy_gev,
                    event.theta_deg,
                    event.phi_deg,
                    event.azimuth_north_to_east_deg,
                    event.core_x_m,
                    event.core_y_m,
                    event.array_rotation_deg,
                    event.h_first_int_m,
                    event.x_max_g_cm2,
                    event.h_max_m,
                    event.starting_grammage_g_cm2,
                    event.ground_gammas,
                    event.ground_electrons,
                    event.ground_hadrons,
                    event.ground_muons,
                });
            }
            hid_t corsika_shower_type =
                H5Tcreate(H5T_COMPOUND, sizeof(CorsikaShowerRow));
            H5Tinsert(corsika_shower_type, "shower_event_id",
                      HOFFSET(CorsikaShowerRow, shower_event_id), H5T_NATIVE_INT32);
            H5Tinsert(corsika_shower_type, "primary_type",
                      HOFFSET(CorsikaShowerRow, primary_type), H5T_NATIVE_INT32);
            H5Tinsert(corsika_shower_type, "energy_gev",
                      HOFFSET(CorsikaShowerRow, energy_gev), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "theta_deg",
                      HOFFSET(CorsikaShowerRow, theta_deg), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "phi_deg",
                      HOFFSET(CorsikaShowerRow, phi_deg), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "azimuth_north_to_east_deg",
                      HOFFSET(CorsikaShowerRow, azimuth_north_to_east_deg),
                      H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "core_x_north_m",
                      HOFFSET(CorsikaShowerRow, core_x_north_m), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "core_y_west_m",
                      HOFFSET(CorsikaShowerRow, core_y_west_m), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "array_rotation_deg",
                      HOFFSET(CorsikaShowerRow, array_rotation_deg), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "h_first_int_m",
                      HOFFSET(CorsikaShowerRow, h_first_int_m), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "x_max_g_cm2",
                      HOFFSET(CorsikaShowerRow, x_max_g_cm2), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "h_max_m",
                      HOFFSET(CorsikaShowerRow, h_max_m), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "starting_grammage_g_cm2",
                      HOFFSET(CorsikaShowerRow, starting_grammage_g_cm2), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "ground_gammas",
                      HOFFSET(CorsikaShowerRow, ground_gammas), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "ground_electrons",
                      HOFFSET(CorsikaShowerRow, ground_electrons), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "ground_hadrons",
                      HOFFSET(CorsikaShowerRow, ground_hadrons), H5T_NATIVE_DOUBLE);
            H5Tinsert(corsika_shower_type, "ground_muons",
                      HOFFSET(CorsikaShowerRow, ground_muons), H5T_NATIVE_DOUBLE);
            writeCompound1D(events_group, "corsika_showers", corsika_shower_type,
                            corsika_shower_rows);
            H5Tclose(corsika_shower_type);
        }
        H5Gclose(events_group);

        hid_t images_group = H5Gcreate2(file, "images", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t image_type = H5Tcreate(H5T_COMPOUND, sizeof(ImageIndexRow));
        H5Tinsert(image_type, "image_index",
                  HOFFSET(ImageIndexRow, image_index), H5T_NATIVE_INT32);
        H5Tinsert(image_type, "event_id", HOFFSET(ImageIndexRow, event_id), H5T_NATIVE_INT64);
        H5Tinsert(image_type, "telescope_id",
                  HOFFSET(ImageIndexRow, telescope_id), H5T_NATIVE_INT32);
        H5Tinsert(image_type, "start", HOFFSET(ImageIndexRow, start), H5T_NATIVE_INT64);
        H5Tinsert(image_type, "count", HOFFSET(ImageIndexRow, count), H5T_NATIVE_INT32);
        H5Tinsert(image_type, "total_photons",
                  HOFFSET(ImageIndexRow, total_photons), H5T_NATIVE_DOUBLE);
        H5Tinsert(image_type, "total_pe", HOFFSET(ImageIndexRow, total_pe), H5T_NATIVE_DOUBLE);
        H5Tinsert(image_type, "total_signal",
                  HOFFSET(ImageIndexRow, total_signal), H5T_NATIVE_DOUBLE);
        H5Tinsert(image_type, "time_mean_ns",
                  HOFFSET(ImageIndexRow, time_mean_ns), H5T_NATIVE_DOUBLE);
        H5Tinsert(image_type, "time_rms_ns",
                  HOFFSET(ImageIndexRow, time_rms_ns), H5T_NATIVE_DOUBLE);
        H5Tinsert(image_type, "time_first_ns",
                  HOFFSET(ImageIndexRow, time_first_ns), H5T_NATIVE_DOUBLE);
        writeCompound1D(images_group, "index", image_type, image_rows);
        H5Tclose(image_type);

        if (write_sparse) {
            hid_t sparse_group = H5Gcreate2(images_group, "sparse",
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            hid_t sparse_type = H5Tcreate(H5T_COMPOUND, sizeof(SparsePixelRow));
            H5Tinsert(sparse_type, "pixel_id",
                      HOFFSET(SparsePixelRow, pixel_id), H5T_NATIVE_INT32);
            H5Tinsert(sparse_type, "photon_count",
                      HOFFSET(SparsePixelRow, photon_count), H5T_NATIVE_INT32);
            H5Tinsert(sparse_type, "pe", HOFFSET(SparsePixelRow, pe), H5T_NATIVE_FLOAT);
            H5Tinsert(sparse_type, "signal", HOFFSET(SparsePixelRow, signal), H5T_NATIVE_FLOAT);
            H5Tinsert(sparse_type, "time_mean_ns",
                      HOFFSET(SparsePixelRow, time_mean_ns), H5T_NATIVE_FLOAT);
            H5Tinsert(sparse_type, "time_rms_ns",
                      HOFFSET(SparsePixelRow, time_rms_ns), H5T_NATIVE_FLOAT);
            writeCompound1D(sparse_group, "pixels", sparse_type, sparse_rows);
            H5Tclose(sparse_type);
            H5Gclose(sparse_group);
        }

        if (have_dense_images) {
            const std::size_t n_images = image_rows.size();
            const std::size_t n_pixels = camera_rows.size();
            hid_t dense_group = H5Gcreate2(images_group, "dense",
                                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            writePlain1D(dense_group, "pixel_id_axis", H5T_NATIVE_INT32, pixel_id_axis);
            writePlain2D(dense_group,
                         "signal",
                         H5T_NATIVE_FLOAT,
                         dense_signal,
                         static_cast<hsize_t>(n_images),
                         static_cast<hsize_t>(n_pixels));
            writePlain2D(dense_group,
                         "pe",
                         H5T_NATIVE_FLOAT,
                         dense_pe,
                         static_cast<hsize_t>(n_images),
                         static_cast<hsize_t>(n_pixels));
            writePlain2D(dense_group,
                         "photon_count",
                         H5T_NATIVE_INT32,
                         dense_photon_count,
                         static_cast<hsize_t>(n_images),
                         static_cast<hsize_t>(n_pixels));
            if (output_cfg.hdf5_write_components) {
                writePlain2D(dense_group,
                             "cherenkov_pe",
                             H5T_NATIVE_FLOAT,
                             dense_cherenkov_pe,
                             static_cast<hsize_t>(n_images),
                             static_cast<hsize_t>(n_pixels));
                writePlain2D(dense_group,
                             "nsb_pe",
                             H5T_NATIVE_FLOAT,
                             dense_nsb_pe,
                             static_cast<hsize_t>(n_images),
                             static_cast<hsize_t>(n_pixels));
            }
            if (output_cfg.write_pixel_time_stats) {
                writePlain2D(dense_group,
                             "time_mean_ns",
                             H5T_NATIVE_FLOAT,
                             dense_time_mean_ns,
                             static_cast<hsize_t>(n_images),
                             static_cast<hsize_t>(n_pixels));
                writePlain2D(dense_group,
                             "time_rms_ns",
                             H5T_NATIVE_FLOAT,
                             dense_time_rms_ns,
                             static_cast<hsize_t>(n_images),
                             static_cast<hsize_t>(n_pixels));
            }
            H5Gclose(dense_group);
        }
        H5Gclose(images_group);

        if (waveform_cfg.enabled && output_cfg.hdf5_write_waveforms &&
            have_camera_axis) {
            std::vector<Hdf5WaveformImage> waveform_images;
            waveform_images.reserve(image_rows.size());
            for (const auto& image : image_rows) {
                waveform_images.push_back(Hdf5WaveformImage{
                    image.image_index,
                    static_cast<int>(image.event_id),
                    static_cast<int>(image.telescope_id),
                    static_cast<double>(image.time_first_ns),
                    static_cast<double>(image.time_mean_ns),
                });
            }
            writeHdf5Waveforms(file,
                               output_cfg,
                               waveform_cfg,
                               nsb_cfg,
                               pixel_id_axis,
                               waveform_images,
                               waveforms,
                               raw_waveform_hits);
        }

        hid_t trigger_group = H5Gcreate2(file, "trigger",
                                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t telescope_trigger_type =
            H5Tcreate(H5T_COMPOUND, sizeof(TelescopeTriggerRow));
        H5Tinsert(telescope_trigger_type, "event_id",
                  HOFFSET(TelescopeTriggerRow, event_id), H5T_NATIVE_INT64);
        H5Tinsert(telescope_trigger_type, "telescope_id",
                  HOFFSET(TelescopeTriggerRow, telescope_id), H5T_NATIVE_INT32);
        H5Tinsert(telescope_trigger_type, "triggered",
                  HOFFSET(TelescopeTriggerRow, triggered), H5T_NATIVE_SCHAR);
        H5Tinsert(telescope_trigger_type, "n_pixels_above_threshold",
                  HOFFSET(TelescopeTriggerRow, n_pixels_above_threshold),
                  H5T_NATIVE_INT32);
        H5Tinsert(telescope_trigger_type, "total_pe",
                  HOFFSET(TelescopeTriggerRow, total_pe), H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_trigger_type, "trigger_time_ns",
                  HOFFSET(TelescopeTriggerRow, trigger_time_ns), H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_trigger_type, "trigger_first_time_ns",
                  HOFFSET(TelescopeTriggerRow, trigger_first_time_ns),
                  H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_trigger_type, "trigger_max_multiplicity_time_ns",
                  HOFFSET(TelescopeTriggerRow,
                          trigger_max_multiplicity_time_ns),
                  H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_trigger_type, "geometric_delay_ns",
                  HOFFSET(TelescopeTriggerRow, geometric_delay_ns),
                  H5T_NATIVE_DOUBLE);
        H5Tinsert(telescope_trigger_type, "coincidence_time_ns",
                  HOFFSET(TelescopeTriggerRow, coincidence_time_ns),
                  H5T_NATIVE_DOUBLE);
        writeCompound1D(trigger_group, "telescope", telescope_trigger_type,
                        telescope_trigger_rows);
        H5Tclose(telescope_trigger_type);

        hid_t array_trigger_type = H5Tcreate(H5T_COMPOUND, sizeof(ArrayTriggerRow));
        H5Tinsert(array_trigger_type, "event_id",
                  HOFFSET(ArrayTriggerRow, event_id), H5T_NATIVE_INT64);
        H5Tinsert(array_trigger_type, "array_triggered",
                  HOFFSET(ArrayTriggerRow, array_triggered), H5T_NATIVE_SCHAR);
        H5Tinsert(array_trigger_type, "n_triggered_telescopes",
                  HOFFSET(ArrayTriggerRow, n_triggered_telescopes), H5T_NATIVE_INT32);
        writeCompound1D(trigger_group, "array", array_trigger_type, array_trigger_rows);
        H5Tclose(array_trigger_type);
        H5Gclose(trigger_group);

        if (!whiteboard_hits.empty()) {
            hid_t whiteboard_group = H5Gcreate2(file, "whiteboard",
                                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            hid_t hit_type = H5Tcreate(H5T_COMPOUND, sizeof(WhiteboardHdf5Row));
            H5Tinsert(hit_type, "event_id",
                      HOFFSET(WhiteboardHdf5Row, event_id), H5T_NATIVE_INT64);
            H5Tinsert(hit_type, "telescope_id",
                      HOFFSET(WhiteboardHdf5Row, telescope_id), H5T_NATIVE_INT32);
            H5Tinsert(hit_type, "photon_index",
                      HOFFSET(WhiteboardHdf5Row, photon_index), H5T_NATIVE_INT64);
            H5Tinsert(hit_type, "mirror_id",
                      HOFFSET(WhiteboardHdf5Row, mirror_id), H5T_NATIVE_INT32);
            H5Tinsert(hit_type, "surface_x_m",
                      HOFFSET(WhiteboardHdf5Row, surface_x_m), H5T_NATIVE_FLOAT);
            H5Tinsert(hit_type, "surface_y_m",
                      HOFFSET(WhiteboardHdf5Row, surface_y_m), H5T_NATIVE_FLOAT);
            H5Tinsert(hit_type, "surface_z_m",
                      HOFFSET(WhiteboardHdf5Row, surface_z_m), H5T_NATIVE_FLOAT);
            H5Tinsert(hit_type, "u_m", HOFFSET(WhiteboardHdf5Row, u_m), H5T_NATIVE_FLOAT);
            H5Tinsert(hit_type, "v_m", HOFFSET(WhiteboardHdf5Row, v_m), H5T_NATIVE_FLOAT);
            H5Tinsert(hit_type, "dir_x",
                      HOFFSET(WhiteboardHdf5Row, dir_x), H5T_NATIVE_FLOAT);
            H5Tinsert(hit_type, "dir_y",
                      HOFFSET(WhiteboardHdf5Row, dir_y), H5T_NATIVE_FLOAT);
            H5Tinsert(hit_type, "dir_z",
                      HOFFSET(WhiteboardHdf5Row, dir_z), H5T_NATIVE_FLOAT);
            H5Tinsert(hit_type, "time_ns",
                      HOFFSET(WhiteboardHdf5Row, time_ns), H5T_NATIVE_FLOAT);
            H5Tinsert(hit_type, "wavelength_nm",
                      HOFFSET(WhiteboardHdf5Row, wavelength_nm), H5T_NATIVE_FLOAT);
            H5Tinsert(hit_type, "weight",
                      HOFFSET(WhiteboardHdf5Row, weight), H5T_NATIVE_FLOAT);
            H5Tinsert(hit_type, "relative_efficiency",
                      HOFFSET(WhiteboardHdf5Row, relative_efficiency), H5T_NATIVE_FLOAT);
            H5Tinsert(hit_type, "signal_weight",
                      HOFFSET(WhiteboardHdf5Row, signal_weight), H5T_NATIVE_FLOAT);
            if (output_cfg.whiteboard_emitter_info) {
                H5Tinsert(hit_type, "has_emitter",
                          HOFFSET(WhiteboardHdf5Row, has_emitter), H5T_NATIVE_UINT8);
                H5Tinsert(hit_type, "emitter_mass_gev",
                          HOFFSET(WhiteboardHdf5Row, emitter_mass_gev), H5T_NATIVE_FLOAT);
                H5Tinsert(hit_type, "emitter_charge",
                          HOFFSET(WhiteboardHdf5Row, emitter_charge), H5T_NATIVE_FLOAT);
                H5Tinsert(hit_type, "emitter_energy_gev",
                          HOFFSET(WhiteboardHdf5Row, emitter_energy_gev), H5T_NATIVE_FLOAT);
                H5Tinsert(hit_type, "emitter_time_ns",
                          HOFFSET(WhiteboardHdf5Row, emitter_time_ns), H5T_NATIVE_FLOAT);
            }
            writeCompound1D(whiteboard_group, "hits", hit_type, whiteboard_hits);
            H5Tclose(hit_type);
            H5Gclose(whiteboard_group);
        }

        H5Fclose(file);
    } catch (...) {
        H5Fclose(file);
        throw;
    }
}
#endif

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
    printField("eventio_reference_z_m",
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
        printField("mode", "whiteboard only");
    }

    printSection("SiPM");
    printField("size_m", doubleToString(sipm_cfg.size_m));
    printField("pde", factorDescription(efficiency_cfg.sipm_pde));

    printSection("Electronics");
    printField("pipeline_enabled", detector_cfg.enabled ? "true" : "false");
    printField("pde_stage", "upstream sipm.pde, applied exactly once");
    printField("microcell_enabled",
               detector_cfg.microcell.enabled ? "true" : "false");
    printField("microcell_model", detector_cfg.microcell.model);
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
        missing = "proxy waveform output, " + missing;
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

} // namespace

int main(int argc, char** argv) {
    if (argc < 2 || argc > 3) {
        std::cerr << "usage: run_corsika_trace <config.txt> [corsika_eventio_file]\n";
        return 2;
    }

    try {
#ifndef LACT_HAS_HESSIO
        throw std::runtime_error(
            "run_corsika_trace requires libhessio. Build external/hessioxxx/source "
            "and reconfigure LACT_sim.");
#else
        std::cout.setf(std::ios::unitbuf);
        std::cout << std::fixed << std::setprecision(6);
        const auto t_start = std::chrono::steady_clock::now();
        auto main_cfg = readKeyValueConfig(argv[1]);
        ComponentConfigPaths component_paths;
        auto cfg = expandConfig(main_cfg, argv[1], component_paths);

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
        if (argc == 3 && !photon_csv_mode) {
            source_runtime_cfg.eventio_path = argv[2];
            cfg["source.eventio_path"] = argv[2];
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
        NsbConfig nsb_cfg = buildNsbConfig(cfg);
        TriggerConfig trigger_cfg = buildTriggerConfig(cfg);
        CameraGeometry camera = buildCameraGeometry(camera_cfg);
        auto light_collector = buildLightCollector(camera_cfg, camera);
        TelescopeOpticsCache telescope_optics(nominal_facets, error_cfg);
        const MirrorLayout& mirrors = telescope_optics.layoutFor(telescope_cfg.id);
        OpticalEfficiencyConfig efficiency_cfg = buildEfficiencyConfig(cfg);
        resolveNsbSpectralRate(nsb_cfg, efficiency_cfg, camera, telescope_cfg);
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
                "electronics.pipeline.enabled=true");
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
            if (waveform_cfg.enabled &&
                waveform_cfg.source == "electronics" &&
                output_cfg.hdf5_write_waveforms) {
                throw std::runtime_error(
                    "electronics waveform serialization is currently "
                    "defined by the lact_event ROOT schema; set "
                    "output.hdf5_write_waveforms=false or use ROOT output");
            }
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
        printField("main", argv[1]);
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
                argv[1],
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

        auto writeRootEventIfReady = [&](int event_id) {
            if (event_id < 0) {
                return;
            }
#ifdef LACT_HAS_ROOT
            if (!lact_root_stream_writer) {
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
            lact_root_stream_writer->writeEvent(event_summaries,
                                                event_pixels,
                                                event_waveforms,
                                                event_raw_waveform_hits);
#endif
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
            writeRootEventIfReady(event_id);
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
                Photon& photon = candidate.photon;
                photon.normalizeDirection();

                double atmosphere_t = 1.0;
                if (atmosphere.enabled() &&
                    !photon.optical_efficiency_preapplied) {
                    atmosphere_t = atmosphere.transmission(
                        photon.wavelength_nm,
                        bunch.emission_altitude_km,
                        global_dir);
                }
                const auto pre_geometry = response_sampler.applyPreGeometry(
                    candidate,
                    atmosphere_t,
                    eff.preGeometryDetectionProbability(photon.wavelength_nm));
                accumulateAtmosphereHistogram(
                    atmosphere_histogram,
                    atmosphere_histogram_cfg,
                    bunch.emission_altitude_km,
                    pre_geometry.expected_weight_before_atmosphere,
                    pre_geometry.expected_weight_after_atmosphere,
                    pre_geometry.expected_weight_after_atmosphere);
                if (!pre_geometry.survives) {
                    return;
                }

            if (profile_cfg.enabled) {
                t_step = std::chrono::steady_clock::now();
            }
            const bool backproject_this_bunch = bunch.eventio_2d && eventio_2d_backproject;
            OpticalSurfaceHit hit = backproject_this_bunch
                                        ? tracer.traceBackprojectedToPlane(
                                              photon, telescope_mirrors, plane, eff)
                                        : tracer.traceToPlane(
                                              photon, telescope_mirrors, plane, eff);
            if (profile_cfg.enabled) {
                addElapsed(profile_stats, &ProfileStats::trace_to_plane_s, t_step);
            }
            if (hit.hit_mirror) {
                summary.hit_mirror_before_obstruction += 1;
                if (hit.hit_surface) {
                    summary.hit_output_before_obstruction += 1;
                }
                if (profile_cfg.enabled) {
                    t_step = std::chrono::steady_clock::now();
                }
                // A 2D EventIO position is only an anchor on the incoming
                // line. Obstruction is checked on the physical upstream ray,
                // independent of which side of the mirror contains the
                // record plane.
                if (incomingRayBlockedByObstruction(hit.mirror_point,
                                                    photon.dir,
                                                    obstruction, nullptr)) {
                    if (profile_cfg.enabled) {
                        addElapsed(profile_stats, &ProfileStats::obstruction_s, t_step);
                    }
                    summary.blocked_by_obstruction += 1;
                    summary.blocked_incoming += 1;
                    if (mirror_diagnostic_out) {
                        writeCorsikaMirrorDiagnosticHit(
                            mirror_diagnostic_out, bunch, photon_index, hit,
                            "blocked_incoming");
                    }
                    return;
                }
                if (profile_cfg.enabled) {
                    addElapsed(profile_stats, &ProfileStats::obstruction_s, t_step);
                }
                summary.hit_mirror += 1;
            }
            if (!hit.hit_surface) {
                if (hit.hit_mirror && mirror_diagnostic_out) {
                    writeCorsikaMirrorDiagnosticHit(
                        mirror_diagnostic_out, bunch, photon_index, hit,
                        "reflected_missed_output");
                }
                return;
            }
            if (profile_cfg.enabled) {
                t_step = std::chrono::steady_clock::now();
            }
            if (segmentBlockedByObstruction(hit.mirror_point, hit.surface_point,
                                            obstruction, nullptr)) {
                if (profile_cfg.enabled) {
                    addElapsed(profile_stats, &ProfileStats::obstruction_s, t_step);
                }
                summary.blocked_by_obstruction += 1;
                summary.blocked_reflected += 1;
                if (mirror_diagnostic_out) {
                    writeCorsikaMirrorDiagnosticHit(
                        mirror_diagnostic_out, bunch, photon_index, hit,
                        "blocked_reflected");
                }
                return;
            }
            if (profile_cfg.enabled) {
                addElapsed(profile_stats, &ProfileStats::obstruction_s, t_step);
            }

            summary.hit_output_plane += 1;
            if (mirror_diagnostic_out) {
                writeCorsikaMirrorDiagnosticHit(
                    mirror_diagnostic_out, bunch, photon_index, hit,
                    "reflected_to_output");
            }

            if (!camera_cfg.enabled) {
                if (candidate.stochastic &&
                    !response_sampler.acceptPostGeometry(candidate, hit)) {
                    return;
                }
                const double base_signal = hit.weight * hit.relative_efficiency;
                if (profile_cfg.enabled) {
                    t_step = std::chrono::steady_clock::now();
                }
                if (save_hdf5) {
                    whiteboard_hits.push_back(makeWhiteboardHdf5Row(bunch, photon_index, hit));
                }
                if (save_csv) {
                    writeCorsikaWhiteboardHit(whiteboard_out, bunch, photon_index, hit,
                                              output_cfg.whiteboard_emitter_info);
                }
                summary.weighted_signal += base_signal;
                summary.weighted_time_sum += base_signal * hit.time_ns;
                summary.weighted_time2_sum += base_signal * hit.time_ns * hit.time_ns;
                if (profile_cfg.enabled) {
                    addElapsed(profile_stats, &ProfileStats::whiteboard_accumulate_s, t_step);
                }
                return;
            }

            if (profile_cfg.enabled) {
                t_step = std::chrono::steady_clock::now();
            }
            applyCameraResponse(camera, light_collector.get(), plane, sipm_cfg,
                                electronics, hit,
                                propagation_cfg.speed_of_light_m_per_ns);
            if (profile_cfg.enabled) {
                addElapsed(profile_stats, &ProfileStats::camera_response_s, t_step);
            }
            appendCollectorDebugPhoton(collector_debug_rows,
                                       collector_debug_cfg,
                                       bunch.event_id,
                                       bunch.telescope_id,
                                       hit);
            if (!hit.hit_camera) {
                summary.lost_between_pixels += 1;
                return;
            }
            summary.hit_camera += 1;
            summary.unique_pixels.insert(hit.pixel_id);
            if (candidate.stochastic &&
                !response_sampler.acceptPostGeometry(candidate, hit)) {
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
                t_step = std::chrono::steady_clock::now();
            }
            accumulatePixelHit(pixels, bunch.event_id, bunch.telescope_id, hit);
            accumulateWaveformHit(waveforms,
                                  raw_waveform_hits,
                                  waveform_cfg,
                                  detector_pipeline_cfg.enabled,
                                  bunch.event_id,
                                  bunch.telescope_id,
                                  hit);
            if (profile_cfg.enabled) {
                addElapsed(profile_stats, &ProfileStats::camera_accumulate_s, t_step);
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
            writePixelCsv(output_cfg.pixel_csv, pixels);
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
                                 argv[1],
                                 cfg,
                                 component_paths,
                                 source_runtime_cfg,
                                 telescope_cfg,
                                 metadata,
                                 camera,
                                 nominal_facets,
                                 sipm_cfg,
                                 electronics_cfg,
                                 efficiency_cfg,
                                 nsb_cfg,
                                 trigger_cfg,
                                 summaries,
                                 pixels,
                                 waveforms,
                                 raw_waveform_hits,
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
            writeSummaryCsv(output_cfg.summary_csv, summaries);
        }
        const auto t_done = std::chrono::steady_clock::now();

        printSection("Input");
        printField("config", argv[1]);
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
