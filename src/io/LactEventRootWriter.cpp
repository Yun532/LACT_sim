#include "io/LactEventRootWriter.hpp"

#include "app/TriggerResponse.hpp"
#include "io/EventIOArrayTiming.hpp"

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace lact {
namespace {

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

void configureRootTreeAutoFlush(TTree* tree, double auto_flush_mb)
{
    if (!tree || auto_flush_mb <= 0.0) {
        return;
    }
    const auto bytes = static_cast<Long64_t>(auto_flush_mb * 1024.0 * 1024.0);
    if (bytes <= 0) {
        return;
    }
    tree->SetAutoFlush(-bytes);
    tree->SetAutoSave(-bytes);
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

bool waveformUsesImageReference(const WaveformOutputConfig& cfg)
{
    return cfg.enabled &&
        (cfg.time_reference == "image_mean" || cfg.time_reference == "image_first");
}

int lactRootPixelShapeCode(PixelShape shape)
{
    if (shape == PixelShape::Square) return 1;
    if (shape == PixelShape::Hexagonal) return 2;
    if (shape == PixelShape::Circular) return 3;
    return 0;
}

std::string lactRootReadTextIfExists(const std::string& path)
{
    if (path.empty() || !std::filesystem::exists(path)) {
        return "";
    }
    std::ifstream ifs(path);
    std::ostringstream ss;
    ss << ifs.rdbuf();
    return ss.str();
}

double mirrorFacetArea(const MirrorFacet& facet)
{
    if (facet.aperture_shape == ApertureShape::Square) {
        return facet.size1 * facet.size1;
    }
    if (facet.aperture_shape == ApertureShape::Hexagon) {
        return 0.5 * std::sqrt(3.0) * facet.size1 * facet.size1;
    }
    return 3.14159265358979323846 * facet.size1 * facet.size1;
}

struct LactRootObservation {
    long long event_id = 0;
    int telescope_id = 0;
    bool triggered = false;
    int n_pixels_camera = 0;
    int n_pixels_saved = 0;
    std::vector<int> pixel_id;
    std::vector<float> image_pe;
    std::vector<float> image_cherenkov_pe;
    std::vector<float> image_fired_cherenkov_pe;
    std::vector<float> image_fired_pe;
    std::vector<float> image_nsb_pe;
    std::vector<float> image_time_mean_ns;
    std::vector<float> image_time_rms_ns;
    std::vector<float> image_time_peak_ns;
    double total_pe = 0.0;
    double time_first_ns = std::numeric_limits<double>::quiet_NaN();
    double time_mean_ns = std::numeric_limits<double>::quiet_NaN();
    double time_rms_ns = std::numeric_limits<double>::quiet_NaN();
    double time_peak_ns = std::numeric_limits<double>::quiet_NaN();
    double impact_parameter_m = std::numeric_limits<double>::quiet_NaN();
    int n_pixels_above_threshold = 0;
    double trigger_time_ns = std::numeric_limits<double>::quiet_NaN();
    double trigger_first_time_ns = std::numeric_limits<double>::quiet_NaN();
    double trigger_max_multiplicity_time_ns =
        std::numeric_limits<double>::quiet_NaN();
    double geometric_delay_ns = std::numeric_limits<double>::quiet_NaN();
    double coincidence_time_ns = std::numeric_limits<double>::quiet_NaN();
};

struct LactRootWaveform {
    long long event_id = 0;
    int telescope_id = 0;
    int n_pixels_camera = 0;
    int n_time_bins = 0;
    std::vector<int> pixel_id;
    std::vector<unsigned short> time_bin;
    std::vector<float> sample_value;
};

struct LactRootPrimaryHit {
    long long event_id = 0;
    int telescope_id = 0;
    int pixel_id = -1;
    double time_ns = 0.0;
    double sensor_x_m = 0.0;
    double sensor_y_m = 0.0;
    double wavelength_nm = 0.0;
    double primary_pe = 0.0;
    int origin = 0;
};

struct LactRootFiredHit {
    long long event_id = 0;
    int telescope_id = 0;
    int pixel_id = -1;
    double time_ns = 0.0;
    int channel_id = -1;
    int microcell_id = -1;
    double fired_pe = 0.0;
    int origin = 0;
};

struct LactRootMicrocellDecision {
    long long event_id = 0;
    int telescope_id = 0;
    int pixel_id = -1;
    double time_ns = 0.0;
    double sensor_x_m = 0.0;
    double sensor_y_m = 0.0;
    int grid_column = -1;
    int grid_row = -1;
    int channel_id = -1;
    int microcell_id = -1;
    bool fired = false;
    int origin = 0;
};

struct LactRootPreparedData {
    std::vector<int> pixel_axis;
    std::vector<double> time_edges_ns;
    std::vector<double> time_centers_ns;
    std::vector<LactRootObservation> observations;
    std::vector<LactRootWaveform> waveforms;
    std::vector<LactRootPrimaryHit> primary_hits;
    std::vector<LactRootFiredHit> fired_hits;
    std::vector<LactRootMicrocellDecision> microcell_decisions;
};

LactRootPreparedData prepareLactRootObservations(
    const CorsikaTraceOutputConfig& output_cfg,
    const WaveformOutputConfig& waveform_cfg,
    const electronics::DetectorPipelineConfig& detector_cfg,
    const NsbConfig& nsb_cfg,
    const TriggerConfig& trigger_cfg,
    const SourceRuntimeConfig& source_runtime_cfg,
    const TelescopeConfig& telescope_cfg,
    const EventIOMetadata& metadata,
    const CameraGeometry& camera,
    const std::map<SummaryKey, TraceSummary>& summaries,
    const std::map<PixelKey, PixelAccumulator>& pixels,
    const std::map<WaveformKey, WaveformPixelAccumulator>& waveforms,
    const std::vector<RawWaveformHit>& raw_waveform_hits)
{
    LactRootPreparedData prepared;
    prepared.pixel_axis.reserve(camera.size());
    std::map<int, std::size_t> pixel_to_col;
    for (std::size_t i = 0; i < camera.pixels().size(); ++i) {
        const int pixel_id = camera.pixels()[i].id;
        prepared.pixel_axis.push_back(pixel_id);
        pixel_to_col[pixel_id] = i;
    }

    std::set<SummaryKey> image_keys;
    for (const auto& kv : summaries) image_keys.insert(kv.first);
    for (const auto& kv : pixels) {
        image_keys.insert({std::get<0>(kv.first), std::get<1>(kv.first)});
    }
    for (const auto& kv : waveforms) {
        image_keys.insert({std::get<0>(kv.first), std::get<1>(kv.first)});
    }
    for (const auto& hit : raw_waveform_hits) {
        image_keys.insert({hit.event_id, hit.telescope_id});
    }

    const bool detector_pipeline_available =
        detector_cfg.enabled && !prepared.pixel_axis.empty();
    const bool waveform_pe_available =
        waveform_cfg.enabled &&
        (waveform_cfg.source == "pe" ||
         waveform_cfg.source == "electronics") &&
        !prepared.pixel_axis.empty();
    const bool write_time_series = waveform_pe_available &&
        (output_cfg.lact_profile == "timeseries_pe" ||
         output_cfg.lact_profile == "debug_full");
    const bool evaluate_time_series =
        detector_pipeline_available ||
        (waveform_pe_available && (write_time_series || trigger_cfg.enabled));
    const std::size_t n_pixels = prepared.pixel_axis.size();
    const std::size_t n_bins =
        evaluate_time_series
            ? (detector_pipeline_available
                   ? static_cast<std::size_t>(std::ceil(
                         (detector_cfg.sampling.end_ns -
                          detector_cfg.sampling.start_ns) /
                         detector_cfg.sampling.width_ns))
                   : waveformBinCount(waveform_cfg))
            : 0;
    if (evaluate_time_series) {
        prepared.time_edges_ns.resize(n_bins + 1);
        prepared.time_centers_ns.resize(n_bins);
        for (std::size_t i = 0; i <= n_bins; ++i) {
            prepared.time_edges_ns[i] =
                (detector_pipeline_available
                     ? detector_cfg.sampling.start_ns
                     : waveform_cfg.time_window_start_ns) +
                static_cast<double>(i) *
                    (detector_pipeline_available
                         ? detector_cfg.sampling.width_ns
                         : waveform_cfg.time_bin_width_ns);
        }
        for (std::size_t i = 0; i < n_bins; ++i) {
            prepared.time_centers_ns[i] =
                0.5 * (prepared.time_edges_ns[i] + prepared.time_edges_ns[i + 1]);
        }
    }

    struct PreparedObservation {
        LactRootObservation observation;
        LactRootWaveform waveform;
        bool has_waveform = false;
    };
    std::vector<PreparedObservation> candidates;
    std::map<int, std::vector<TelescopeTriggerTime>>
        telescope_trigger_times_by_event;
    candidates.reserve(image_keys.size());
    prepared.observations.reserve(image_keys.size());
    prepared.waveforms.reserve(write_time_series ? image_keys.size() : 0);
    for (const auto& key : image_keys) {
        const int event_id = key.first;
        const int telescope_id = key.second;
        LactRootObservation obs;
        obs.event_id = event_id;
        obs.telescope_id = telescope_id;
        obs.n_pixels_camera = static_cast<int>(n_pixels);

        std::vector<double> image_pe_by_col(n_pixels, 0.0);
        std::vector<double> image_cherenkov_pe_by_col(n_pixels, 0.0);
        std::vector<double> image_fired_cherenkov_pe_by_col(n_pixels, 0.0);
        std::vector<double> image_fired_pe_by_col(n_pixels, 0.0);
        std::vector<double> image_nsb_pe_by_col(n_pixels, 0.0);
        std::vector<double> time_sum_by_col(n_pixels, 0.0);
        std::vector<double> time2_sum_by_col(n_pixels, 0.0);
        std::vector<double> time_weight_by_col(n_pixels, 0.0);
        std::unordered_map<std::size_t, double> waveform_pe;
        std::vector<double> waveform_peak_by_col;
        std::vector<std::size_t> waveform_peak_bin_by_col;
        double reference_time_ns = 0.0;
        double trigger_time_ns = std::numeric_limits<double>::quiet_NaN();
        bool detector_camera_triggered = false;

        auto summary_it = summaries.find(key);
        if (summary_it != summaries.end()) {
            const auto& s = summary_it->second;
            if (std::isfinite(s.first_cherenkov_time_ns)) {
                obs.time_first_ns = s.first_cherenkov_time_ns;
            }
            if (s.weighted_signal > 0.0) {
                const double mean = s.weighted_time_sum / s.weighted_signal;
                const double var =
                    std::max(0.0, s.weighted_time2_sum / s.weighted_signal - mean * mean);
                obs.time_mean_ns = mean;
                obs.time_rms_ns = std::sqrt(var);
            }
        }

        // Pixel timing moments describe the Cherenkov signal itself. They are
        // independent of waveform binning and must not be diluted by NSB.
        const PixelKey pixel_begin_key{
            event_id, telescope_id, std::numeric_limits<int>::min()};
        const PixelKey pixel_end_key{
            event_id, telescope_id, std::numeric_limits<int>::max()};
        for (auto it = pixels.lower_bound(pixel_begin_key);
             it != pixels.end() && it->first <= pixel_end_key;
             ++it) {
            const auto& p = it->second;
            const auto col_it = pixel_to_col.find(p.pixel_id);
            if (col_it == pixel_to_col.end()) continue;
            const std::size_t col = col_it->second;
            time_weight_by_col[col] = p.signal;
            time_sum_by_col[col] = p.time_sum;
            time2_sum_by_col[col] = p.time2_sum;
        }

        if (detector_pipeline_available) {
            if (waveform_cfg.time_reference == "image_first" &&
                std::isfinite(obs.time_first_ns)) {
                reference_time_ns = obs.time_first_ns;
            } else if (waveform_cfg.time_reference == "image_mean" &&
                       std::isfinite(obs.time_mean_ns)) {
                reference_time_ns = obs.time_mean_ns;
            }

            std::vector<electronics::PrimaryPeHit> primary_hits;
            for (const auto& hit : raw_waveform_hits) {
                if (hit.event_id != event_id ||
                    hit.telescope_id != telescope_id) {
                    continue;
                }
                const auto col_it = pixel_to_col.find(hit.pixel_id);
                if (col_it == pixel_to_col.end()) continue;
                primary_hits.push_back({
                    event_id,
                    telescope_id,
                    static_cast<int>(col_it->second),
                    hit.time_ns - reference_time_ns,
                    hit.sensor_x_m,
                    hit.sensor_y_m,
                    hit.wavelength_nm,
                    hit.pe,
                    hit.origin,
                });
            }
            if (nsb_cfg.enabled &&
                nsb_cfg.rate_pe_per_ns_per_pixel > 0.0) {
                auto nsb_hits = electronics::generateUniformNsbPrimaryHits(
                    event_id,
                    telescope_id,
                    n_pixels,
                    nsb_cfg.rate_pe_per_ns_per_pixel,
                    detector_cfg.sampling.start_ns,
                    detector_cfg.sampling.end_ns,
                    detector_cfg.microcell.sensor_size_x_m,
                    detector_cfg.microcell.sensor_size_y_m,
                    nsb_cfg.seed);
                primary_hits.insert(primary_hits.end(),
                                    nsb_hits.begin(), nsb_hits.end());
            }

            const auto detector_result =
                electronics::runDetectorPipeline(
                    detector_cfg, n_pixels, primary_hits);
            std::vector<double> camera_time_series(n_bins, 0.0);
            waveform_peak_by_col.assign(n_pixels, -1.0);
            waveform_peak_bin_by_col.assign(n_pixels, 0);
            for (std::size_t col = 0; col < n_pixels; ++col) {
                const auto& pixel = detector_result.pixels[col];
                image_cherenkov_pe_by_col[col] =
                    pixel.primary_cherenkov_pe;
                image_nsb_pe_by_col[col] = pixel.primary_nsb_pe;
                image_fired_cherenkov_pe_by_col[col] =
                    pixel.fired_cherenkov_pe;
                image_fired_pe_by_col[col] =
                    pixel.fired_cherenkov_pe +
                    pixel.fired_nsb_pe +
                    pixel.fired_dark_pe;
                image_pe_by_col[col] = image_fired_pe_by_col[col];
            }
            if (detector_cfg.save_primary_sequence) {
                for (const auto& hit : detector_result.primary_hits) {
                    prepared.primary_hits.push_back({
                        hit.event_id,
                        hit.telescope_id,
                        prepared.pixel_axis.at(
                            static_cast<std::size_t>(hit.pixel_id)),
                        hit.time_ns,
                        hit.sensor_x_m,
                        hit.sensor_y_m,
                        hit.wavelength_nm,
                        hit.primary_pe,
                        static_cast<int>(hit.origin),
                    });
                }
            }
            if (detector_cfg.save_fired_sequence) {
                for (const auto& hit : detector_result.fired_hits) {
                    prepared.fired_hits.push_back({
                        hit.event_id,
                        hit.telescope_id,
                        prepared.pixel_axis.at(
                            static_cast<std::size_t>(hit.pixel_id)),
                        hit.time_ns,
                        hit.channel_id,
                        hit.microcell_id,
                        hit.fired_pe,
                        static_cast<int>(hit.origin),
                    });
                }
            }
            if (detector_cfg.save_microcell_decisions) {
                for (const auto& item :
                     detector_result.microcell_decisions) {
                    prepared.microcell_decisions.push_back({
                        item.event_id,
                        item.telescope_id,
                        prepared.pixel_axis.at(
                            static_cast<std::size_t>(item.pixel_id)),
                        item.time_ns,
                        item.sensor_x_m,
                        item.sensor_y_m,
                        item.grid_column,
                        item.grid_row,
                        item.channel_id,
                        item.microcell_id,
                        item.fired,
                        static_cast<int>(item.origin),
                    });
                }
            }

            if (detector_cfg.single_pe.enabled) {
                for (std::size_t col = 0; col < n_pixels; ++col) {
                    for (std::size_t bin = 0; bin < n_bins; ++bin) {
                        const double value =
                            detector_result.waveform[col * n_bins + bin];
                        if (value == 0.0) continue;
                        waveform_pe[col * n_bins + bin] = value;
                        camera_time_series[bin] += value;
                    }
                }
            } else {
                for (const auto& hit : detector_result.fired_hits) {
                    const int bin = static_cast<int>(std::floor(
                        (hit.time_ns - detector_cfg.sampling.start_ns) /
                        detector_cfg.sampling.width_ns));
                    if (bin < 0 || static_cast<std::size_t>(bin) >= n_bins) {
                        continue;
                    }
                    const std::size_t col =
                        static_cast<std::size_t>(hit.pixel_id);
                    waveform_pe[
                        col * n_bins + static_cast<std::size_t>(bin)] +=
                        hit.fired_pe;
                    camera_time_series[static_cast<std::size_t>(bin)] +=
                        hit.fired_pe;
                }
            }
            for (const auto& entry : waveform_pe) {
                const std::size_t col = entry.first / n_bins;
                const std::size_t bin = entry.first % n_bins;
                if (entry.second > waveform_peak_by_col[col]) {
                    waveform_peak_by_col[col] = entry.second;
                    waveform_peak_bin_by_col[col] = bin;
                }
            }
            const auto camera_peak = std::max_element(
                camera_time_series.begin(), camera_time_series.end());
            if (camera_peak != camera_time_series.end() &&
                *camera_peak > 0.0) {
                const std::size_t bin = static_cast<std::size_t>(
                    std::distance(camera_time_series.begin(), camera_peak));
                obs.time_peak_ns =
                    reference_time_ns + prepared.time_centers_ns[bin];
            }
            obs.n_pixels_above_threshold =
                detector_result.camera_trigger.max_pixels_above_threshold;
            detector_camera_triggered =
                detector_result.camera_trigger.triggered;
            if (detector_result.camera_trigger.triggered) {
                trigger_time_ns =
                    reference_time_ns +
                    detector_result.camera_trigger.trigger_time_ns;
                obs.trigger_first_time_ns = trigger_time_ns;
                obs.trigger_max_multiplicity_time_ns = trigger_time_ns;
            }
        } else if (evaluate_time_series) {
            if (waveform_cfg.time_reference == "image_first" &&
                std::isfinite(obs.time_first_ns)) {
                reference_time_ns = obs.time_first_ns;
            } else if (waveform_cfg.time_reference == "image_mean" &&
                       std::isfinite(obs.time_mean_ns)) {
                reference_time_ns = obs.time_mean_ns;
            }

            std::vector<double> camera_time_series(n_bins, 0.0);
            auto add_waveform_pe = [&](std::size_t col,
                                       std::size_t bin,
                                       double pe,
                                       double cherenkov_pe,
                                       double nsb_pe) {
                if (col >= n_pixels || bin >= n_bins || pe == 0.0) {
                    return;
                }
                const std::size_t index = col * n_bins + bin;
                waveform_pe[index] += pe;
                image_pe_by_col[col] += pe;
                image_cherenkov_pe_by_col[col] += cherenkov_pe;
                image_nsb_pe_by_col[col] += nsb_pe;
                camera_time_series[bin] += pe;
            };

            if (waveformUsesImageReference(waveform_cfg)) {
                for (const auto& hit : raw_waveform_hits) {
                    if (hit.event_id != event_id || hit.telescope_id != telescope_id) continue;
                    const auto col_it = pixel_to_col.find(hit.pixel_id);
                    if (col_it == pixel_to_col.end()) continue;
                    const int bin = waveformBinForTime(
                        waveform_cfg, hit.time_ns - reference_time_ns);
                    if (bin < 0) continue;
                    add_waveform_pe(col_it->second,
                                    static_cast<std::size_t>(bin),
                                    hit.pe,
                                    hit.pe,
                                    0.0);
                }
            } else {
                const WaveformKey begin_key{
                    event_id, telescope_id, std::numeric_limits<int>::min(),
                    std::numeric_limits<int>::min()};
                const WaveformKey end_key{
                    event_id, telescope_id, std::numeric_limits<int>::max(),
                    std::numeric_limits<int>::max()};
                for (auto it = waveforms.lower_bound(begin_key);
                     it != waveforms.end() && it->first <= end_key;
                     ++it) {
                    const auto& w = it->second;
                    const auto col_it = pixel_to_col.find(w.pixel_id);
                    if (col_it == pixel_to_col.end() || w.time_bin < 0 ||
                        static_cast<std::size_t>(w.time_bin) >= n_bins) {
                        continue;
                    }
                    add_waveform_pe(col_it->second,
                                    static_cast<std::size_t>(w.time_bin),
                                    w.pe,
                                    w.pe,
                                    0.0);
                }
            }

            const auto waveform_pe_at = [&](std::size_t col, std::size_t bin) {
                const auto it = waveform_pe.find(col * n_bins + bin);
                return it == waveform_pe.end() ? 0.0 : it->second;
            };

            // Triggering, the integrated image and serialized waveforms must
            // all see exactly the same deterministic NSB realization.
            generateTimeBinnedNsbPe(
                nsb_cfg,
                waveform_cfg,
                event_id,
                telescope_id,
                n_pixels,
                n_bins,
                [&](std::size_t col, std::size_t bin, float nsb_pe) {
                    add_waveform_pe(col, bin, nsb_pe, 0.0, nsb_pe);
                });

            waveform_peak_by_col.assign(n_pixels, -1.0);
            waveform_peak_bin_by_col.assign(n_pixels, 0);
            for (const auto& entry : waveform_pe) {
                const std::size_t col = entry.first / n_bins;
                const std::size_t bin = entry.first % n_bins;
                const double pe = entry.second;
                if (pe > waveform_peak_by_col[col] ||
                    (pe == waveform_peak_by_col[col] &&
                     bin < waveform_peak_bin_by_col[col])) {
                    waveform_peak_by_col[col] = pe;
                    waveform_peak_bin_by_col[col] = bin;
                }
            }

            double camera_peak = -1.0;
            std::size_t camera_peak_bin = 0;
            for (std::size_t bin = 0; bin < n_bins; ++bin) {
                if (camera_time_series[bin] > camera_peak) {
                    camera_peak = camera_time_series[bin];
                    camera_peak_bin = bin;
                }
            }
            if (camera_peak > 0.0) {
                obs.time_peak_ns =
                    reference_time_ns + prepared.time_centers_ns[camera_peak_bin];
            }
            const auto final_trigger = evaluateBinnedPeTrigger(
                n_pixels,
                n_bins,
                waveform_cfg.time_bin_width_ns,
                reference_time_ns + prepared.time_centers_ns.front(),
                trigger_cfg,
                waveform_pe_at);
            obs.n_pixels_above_threshold =
                final_trigger.n_pixels_above_threshold;
            trigger_time_ns = final_trigger.trigger_time_ns;
            obs.trigger_first_time_ns =
                final_trigger.first_trigger_time_ns;
            obs.trigger_max_multiplicity_time_ns =
                final_trigger.max_multiplicity_time_ns;
        } else {
            for (auto it = pixels.lower_bound(pixel_begin_key);
                 it != pixels.end() && it->first <= pixel_end_key;
                 ++it) {
                const auto& p = it->second;
                const auto col_it = pixel_to_col.find(p.pixel_id);
                if (col_it == pixel_to_col.end()) continue;
                const std::size_t col = col_it->second;
                image_pe_by_col[col] = p.pe;
                image_cherenkov_pe_by_col[col] = p.pe;
            }
            generateIntegratedNsbPe(
                nsb_cfg,
                event_id,
                telescope_id,
                n_pixels,
                nsb_cfg.window_ns,
                [&](std::size_t col, float nsb_pe) {
                    image_pe_by_col[col] += nsb_pe;
                    image_nsb_pe_by_col[col] += nsb_pe;
                });
        }

        for (std::size_t col = 0; col < n_pixels; ++col) {
            const double pe = image_pe_by_col[col];
            const double primary_pe =
                image_cherenkov_pe_by_col[col] + image_nsb_pe_by_col[col];
            if (pe <= 0.0 && primary_pe <= 0.0) continue;
            obs.pixel_id.push_back(prepared.pixel_axis[col]);
            obs.image_pe.push_back(static_cast<float>(pe));
            obs.image_cherenkov_pe.push_back(
                static_cast<float>(image_cherenkov_pe_by_col[col]));
            if (detector_pipeline_available) {
                obs.image_fired_cherenkov_pe.push_back(
                    static_cast<float>(
                        image_fired_cherenkov_pe_by_col[col]));
                obs.image_fired_pe.push_back(
                    static_cast<float>(image_fired_pe_by_col[col]));
            }
            if (output_cfg.lact_root_write_components) {
                obs.image_nsb_pe.push_back(static_cast<float>(image_nsb_pe_by_col[col]));
            }
            const double time_weight = time_weight_by_col[col];
            const double mean = time_weight > 0.0
                ? time_sum_by_col[col] / time_weight
                : std::numeric_limits<double>::quiet_NaN();
            double rms = std::numeric_limits<double>::quiet_NaN();
            if (time_weight > 0.0 && std::isfinite(mean)) {
                const double var = std::max(
                    0.0, time2_sum_by_col[col] / time_weight - mean * mean);
                rms = std::sqrt(var);
            }
            obs.image_time_mean_ns.push_back(static_cast<float>(mean));
            obs.image_time_rms_ns.push_back(static_cast<float>(rms));
            if (evaluate_time_series && col < waveform_peak_by_col.size()) {
                const double peak = waveform_peak_by_col[col];
                const std::size_t peak_bin = waveform_peak_bin_by_col[col];
                obs.image_time_peak_ns.push_back(static_cast<float>(
                    peak > 0.0
                        ? reference_time_ns + prepared.time_centers_ns[peak_bin]
                        : std::numeric_limits<double>::quiet_NaN()));
            } else {
                obs.image_time_peak_ns.push_back(std::numeric_limits<float>::quiet_NaN());
            }
            obs.total_pe += pe;
            if (!evaluate_time_series && pe >= trigger_cfg.pixel_threshold_pe) {
                ++obs.n_pixels_above_threshold;
            }
        }

        obs.n_pixels_saved = static_cast<int>(obs.pixel_id.size());
        obs.triggered = trigger_cfg.enabled &&
            (detector_pipeline_available
                 ? detector_camera_triggered
                 : obs.n_pixels_above_threshold >=
                       trigger_cfg.camera_multiplicity);
        if (obs.triggered) {
            obs.trigger_time_ns =
                std::isfinite(trigger_time_ns) ? trigger_time_ns :
                (std::isfinite(obs.time_peak_ns) ? obs.time_peak_ns :
                 (std::isfinite(obs.time_mean_ns) ? obs.time_mean_ns : 0.0));
            if (!std::isfinite(obs.trigger_first_time_ns)) {
                obs.trigger_first_time_ns = obs.trigger_time_ns;
            }
            if (!std::isfinite(obs.trigger_max_multiplicity_time_ns)) {
                obs.trigger_max_multiplicity_time_ns = obs.trigger_time_ns;
            }
            telescope_trigger_times_by_event[event_id].push_back(
                TelescopeTriggerTime{telescope_id, obs.trigger_time_ns});
        }

        PreparedObservation candidate;
        if (write_time_series && !waveform_pe.empty()) {
            LactRootWaveform wf;
            wf.event_id = event_id;
            wf.telescope_id = telescope_id;
            wf.n_pixels_camera = static_cast<int>(n_pixels);
            wf.n_time_bins = static_cast<int>(n_bins);
            std::vector<std::size_t> waveform_indices;
            waveform_indices.reserve(waveform_pe.size());
            for (const auto& entry : waveform_pe) {
                waveform_indices.push_back(entry.first);
            }
            std::sort(waveform_indices.begin(), waveform_indices.end());
            for (const std::size_t index : waveform_indices) {
                const double pe = waveform_pe[index];
                if (pe <= 0.0) continue;
                const std::size_t col = index / n_bins;
                const std::size_t bin = index % n_bins;
                wf.pixel_id.push_back(prepared.pixel_axis[col]);
                wf.time_bin.push_back(static_cast<unsigned short>(bin));
                wf.sample_value.push_back(static_cast<float>(pe));
            }
            candidate.waveform = std::move(wf);
            candidate.has_waveform = true;
        }

        candidate.observation = std::move(obs);
        candidates.push_back(std::move(candidate));
    }

    std::map<int, ArrayTriggerDecision> array_trigger_decisions;
    std::map<SummaryKey, TelescopeTriggerTime> array_trigger_times;
    for (auto& item : telescope_trigger_times_by_event) {
        applyEventIOArrayTimingCorrection(
            item.second, item.first, source_runtime_cfg.event_id_mode, trigger_cfg,
            telescope_cfg, metadata);
        for (const auto& trigger_time : item.second) {
            array_trigger_times[{item.first, trigger_time.telescope_id}] =
                trigger_time;
        }
        array_trigger_decisions[item.first] =
            evaluateArrayTrigger(item.second, trigger_cfg);
    }

    for (auto& candidate : candidates) {
        const auto corrected_time = array_trigger_times.find({
            static_cast<int>(candidate.observation.event_id),
            candidate.observation.telescope_id});
        if (corrected_time != array_trigger_times.end()) {
            candidate.observation.geometric_delay_ns =
                corrected_time->second.geometric_delay_ns;
            candidate.observation.coincidence_time_ns =
                std::isfinite(corrected_time->second.coincidence_time_ns)
                    ? corrected_time->second.coincidence_time_ns
                    : corrected_time->second.trigger_time_ns;
        }
        if (output_cfg.save_only_triggered && trigger_cfg.enabled) {
            const auto& array_decision =
                array_trigger_decisions[candidate.observation.event_id];
            const bool telescope_is_coincident = std::binary_search(
                array_decision.coincident_telescope_ids.begin(),
                array_decision.coincident_telescope_ids.end(),
                candidate.observation.telescope_id);
            if (!array_decision.triggered ||
                !candidate.observation.triggered ||
                !telescope_is_coincident) {
                continue;
            }
        }
        if (candidate.has_waveform) {
            prepared.waveforms.push_back(std::move(candidate.waveform));
        }
        prepared.observations.push_back(std::move(candidate.observation));
    }

    return prepared;
}

} // namespace

struct LactEventRootStreamWriter::Impl {
    CorsikaTraceOutputConfig output_cfg;
    WaveformOutputConfig waveform_cfg;
    SourceRuntimeConfig source_runtime_cfg;
    TelescopeConfig telescope_cfg;
    EventIOMetadata metadata;
    CameraGeometry camera;
    electronics::DetectorPipelineConfig detector_cfg;
    NsbConfig nsb_cfg;
    TriggerConfig trigger_cfg;
    std::unique_ptr<TFile> file;
    std::unique_ptr<TTree> corsika_tree;
    std::unique_ptr<TTree> observation_tree;
    std::unique_ptr<TTree> waveform_config_tree;
    std::unique_ptr<TTree> waveform_tree;
    std::unique_ptr<TTree> primary_hit_tree;
    std::unique_ptr<TTree> fired_hit_tree;
    std::unique_ptr<TTree> microcell_decision_tree;
    std::unique_ptr<TTree> trace_tree;
    bool finished = false;
    bool waveform_config_written = false;
    int events_since_flush = 0;
    std::set<long long> written_events;

    int run_id = 0;

    long long root_event_id = 0;
    int shower_event_id = 0, array_id = 0, primary_type = 0;
    bool has_simtel_mc_shower = false;
    double energy_gev = 0.0, theta_deg = 0.0, phi_deg = 0.0;
    double azimuth_north_to_east_deg = 0.0, altitude_deg = 0.0;
    double core_x_north_m = 0.0, core_y_west_m = 0.0, array_rotation_deg = 0.0;
    double array_time_offset_ns = 0.0, area_weight_m2 = 0.0;
    bool has_explicit_area_weight = false;
    double h_first_int_m = std::numeric_limits<double>::quiet_NaN();
    double x_max_g_cm2 = std::numeric_limits<double>::quiet_NaN();
    double h_max_m = std::numeric_limits<double>::quiet_NaN();
    double starting_grammage_g_cm2 = std::numeric_limits<double>::quiet_NaN();
    double ground_gammas = std::numeric_limits<double>::quiet_NaN();
    double ground_electrons = std::numeric_limits<double>::quiet_NaN();
    double ground_hadrons = std::numeric_limits<double>::quiet_NaN();
    double ground_muons = std::numeric_limits<double>::quiet_NaN();

    long long obs_event_id = 0;
    int obs_telescope_id = 0;
    bool triggered = false;
    int n_pixels_camera = 0, n_pixels_saved = 0;
    std::vector<int> obs_pixel_id;
    std::vector<float> image_pe;
    std::vector<float> image_cherenkov_pe;
    std::vector<float> image_fired_cherenkov_pe;
    std::vector<float> image_fired_pe;
    std::vector<float> image_nsb_pe;
    std::vector<float> image_time_mean_ns;
    std::vector<float> image_time_rms_ns;
    std::vector<float> image_time_peak_ns;
    double total_pe = 0.0, time_first_ns = 0.0, time_mean_ns = 0.0;
    double time_rms_ns = 0.0, time_peak_ns = 0.0, impact_parameter_m = 0.0;
    int n_pixels_above_threshold = 0;
    double trigger_time_ns = 0.0;
    double trigger_first_time_ns = 0.0;
    double trigger_max_multiplicity_time_ns = 0.0;
    double geometric_delay_ns = 0.0, coincidence_time_ns = 0.0;

    bool waveform_enabled = true;
    std::string waveform_source;
    std::string waveform_sample_unit;
    std::string time_reference;
    double time_bin_width_ns = 0.0;
    double time_window_start_ns = 0.0;
    double time_window_end_ns = 0.0;
    int n_time_bins = 0;
    std::vector<double> time_edges_ns;
    std::vector<double> time_centers_ns;

    long long wf_event_id = 0;
    int wf_telescope_id = 0, wf_n_pixels_camera = 0, wf_n_time_bins = 0;
    std::vector<int> wf_pixel_id;
    std::vector<unsigned short> wf_time_bin;
    std::vector<float> wf_sample_value;

    LactRootPrimaryHit primary_hit_row;
    LactRootFiredHit fired_hit_row;
    LactRootMicrocellDecision microcell_decision_row;

    long long trace_event_id = 0;
    int trace_telescope_id = 0;
    unsigned long long input_bunches = 0;
    double input_photons = 0.0;
    unsigned long long blocked_by_obstruction = 0, blocked_incoming = 0;
    unsigned long long blocked_reflected = 0, hit_mirror = 0, hit_output_plane = 0;
    unsigned long long hit_camera = 0, accepted_camera = 0, lost_between_pixels = 0;
    int unique_hit_pixels = 0;
    double signal_pe = 0.0, trace_time_mean_ns = 0.0, trace_time_rms_ns = 0.0;

    Impl(const CorsikaTraceOutputConfig &output_cfg_in,
         const WaveformOutputConfig &waveform_cfg_in,
         const std::string &main_config_path,
         const std::map<std::string, std::string> &cfg,
         const SourceRuntimeConfig &source_runtime_cfg_in,
         const TelescopeConfig &telescope_cfg_in,
         const EventIOMetadata &metadata_in, const CameraGeometry &camera_in,
         const std::vector<MirrorFacet> &facets, const NsbConfig &nsb_cfg_in,
         const TriggerConfig &trigger_cfg_in)
        : output_cfg(output_cfg_in), waveform_cfg(waveform_cfg_in),
          source_runtime_cfg(source_runtime_cfg_in),
          telescope_cfg(telescope_cfg_in), metadata(metadata_in),
          camera(camera_in), detector_cfg(buildDetectorPipelineConfig(cfg)),
          nsb_cfg(nsb_cfg_in), trigger_cfg(trigger_cfg_in) {
      const std::filesystem::path out_path(output_cfg.lact_root_path);
      if (out_path.has_parent_path()) {
        std::filesystem::create_directories(out_path.parent_path());
      }

      file.reset(TFile::Open(output_cfg.lact_root_path.c_str(), "RECREATE"));
      if (!file || file->IsZombie()) {
        throw std::runtime_error("failed to create lact_event ROOT file: " +
                                 output_cfg.lact_root_path);
      }

      std::string schema_name = "lact_event_root";
      int schema_version = 1;
      std::string profile = output_cfg.lact_profile;
      std::string producer = "LACT_sim";
      std::string producer_version = LACT_PRODUCER_VERSION;
      std::string source_kind =
          source_runtime_cfg.use_photon_csv ? "PhotonCsv" : "EventIO";
      std::string source_path = source_runtime_cfg.use_photon_csv
                                    ? source_runtime_cfg.csv_path
                                    : source_runtime_cfg.eventio_path;
      const auto source_hash_it = cfg.find("provenance.source_sha256");
      std::string source_sha256 =
          source_hash_it == cfg.end() ? std::string{} : source_hash_it->second;
      run_id = 0;
      std::string coordinate_convention =
          "CORSIKA IACT NWU; azimuth North-to-East";
      std::string event_id_mode = source_runtime_cfg.event_id_mode;
      std::string config_text = lactRootReadTextIfExists(main_config_path);
      std::ostringstream expanded;
      for (const auto &kv : cfg)
        expanded << kv.first << '=' << kv.second << '\n';
      std::string expanded_config_text = expanded.str();

      TTree config_tree("config", "LACT_sim lact_event ROOT configuration");
      config_tree.Branch("schema_name", &schema_name);
      config_tree.Branch("schema_version", &schema_version);
      config_tree.Branch("profile", &profile);
      config_tree.Branch("producer", &producer);
      config_tree.Branch("producer_version", &producer_version);
      config_tree.Branch("source_kind", &source_kind);
      config_tree.Branch("source_path", &source_path);
      config_tree.Branch("source_sha256", &source_sha256);
      config_tree.Branch("run_id", &run_id);
      config_tree.Branch("coordinate_convention", &coordinate_convention);
      config_tree.Branch("event_id_mode", &event_id_mode);
      config_tree.Branch("config_text", &config_text);
      config_tree.Branch("expanded_config_text", &expanded_config_text);
      config_tree.Fill();
      config_tree.Write();

      int camera_id = 0;
      TTree camera_tree("camera_pixels", "Camera pixel geometry");
      int pixel_id = 0;
      double pixel_x_m = 0.0, pixel_y_m = 0.0, pixel_size_m = 0.0;
      int pixel_shape_code = 0;
      camera_tree.Branch("camera_id", &camera_id);
      camera_tree.Branch("pixel_id", &pixel_id);
      camera_tree.Branch("x_m", &pixel_x_m);
      camera_tree.Branch("y_m", &pixel_y_m);
      camera_tree.Branch("size_m", &pixel_size_m);
      camera_tree.Branch("shape_code", &pixel_shape_code);
      for (const auto &pixel : camera.pixels()) {
        pixel_id = pixel.id;
        pixel_x_m = pixel.center.x;
        pixel_y_m = pixel.center.y;
        pixel_size_m = pixel.size;
        pixel_shape_code = lactRootPixelShapeCode(pixel.shape);
        camera_tree.Fill();
      }
      camera_tree.Write();

      TTree optics_tree("optics", "Optics descriptions");
      int optics_id = 0;
      std::string optics_name =
          telescope_cfg.name.empty() ? "LACT" : telescope_cfg.name;
      int num_mirrors = static_cast<int>(facets.size());
      double mirror_area_m2 = 0.0;
      for (const auto &facet : facets)
        mirror_area_m2 += mirrorFacetArea(facet);
      double equivalent_focal_length_m = telescope_cfg.focal_length_m;
      double effective_focal_length_m = telescope_cfg.focal_length_m;
      optics_tree.Branch("optics_id", &optics_id);
      optics_tree.Branch("name", &optics_name);
      optics_tree.Branch("num_mirrors", &num_mirrors);
      optics_tree.Branch("mirror_area_m2", &mirror_area_m2);
      optics_tree.Branch("equivalent_focal_length_m",
                         &equivalent_focal_length_m);
      optics_tree.Branch("effective_focal_length_m", &effective_focal_length_m);
      optics_tree.Fill();
      optics_tree.Write();

      TTree telescope_tree("telescopes", "Telescope positions and pointing");
      int telescope_id = 0;
      std::string telescope_name = telescope_cfg.name;
      double array_x_north_m = 0.0, array_y_west_m = 0.0, array_z_up_m = 0.0;
      double radius_m = 0.0;
      double pointing_az_deg = telescope_cfg.pointing_az_deg;
      double pointing_el_deg = telescope_cfg.pointing_el_deg;
      telescope_tree.Branch("telescope_id", &telescope_id);
      telescope_tree.Branch("name", &telescope_name);
      telescope_tree.Branch("array_x_north_m", &array_x_north_m);
      telescope_tree.Branch("array_y_west_m", &array_y_west_m);
      telescope_tree.Branch("array_z_up_m", &array_z_up_m);
      telescope_tree.Branch("radius_m", &radius_m);
      telescope_tree.Branch("pointing_az_deg", &pointing_az_deg);
      telescope_tree.Branch("pointing_el_deg", &pointing_el_deg);
      telescope_tree.Branch("camera_id", &camera_id);
      telescope_tree.Branch("optics_id", &optics_id);
      if (!metadata.telescopes.empty()) {
        for (const auto &tel : metadata.telescopes) {
          telescope_id = tel.telescope_id;
          array_x_north_m = tel.x_m;
          array_y_west_m = tel.y_m;
          array_z_up_m = tel.z_m;
          radius_m = tel.radius_m;
          telescope_tree.Fill();
        }
      } else {
        telescope_id = telescope_cfg.id;
        array_x_north_m = telescope_cfg.position_m.x;
        array_y_west_m = telescope_cfg.position_m.y;
        array_z_up_m = telescope_cfg.position_m.z;
        radius_m = 0.0;
        telescope_tree.Fill();
      }
      telescope_tree.Write();

      corsika_tree = std::make_unique<TTree>("corsika_events",
                                             "CORSIKA event truth/provenance");
      corsika_tree->SetDirectory(file.get());
      corsika_tree->Branch("event_id", &root_event_id);
      corsika_tree->Branch("shower_event_id", &shower_event_id);
      corsika_tree->Branch("array_id", &array_id);
      corsika_tree->Branch("array_time_offset_ns", &array_time_offset_ns);
      corsika_tree->Branch("area_weight_m2", &area_weight_m2);
      corsika_tree->Branch("has_explicit_area_weight", &has_explicit_area_weight);
      corsika_tree->Branch("run_id", &run_id);
      corsika_tree->Branch("primary_type", &primary_type);
      corsika_tree->Branch("energy_gev", &energy_gev);
      corsika_tree->Branch("theta_deg", &theta_deg);
      corsika_tree->Branch("phi_deg", &phi_deg);
      corsika_tree->Branch("azimuth_north_to_east_deg",
                           &azimuth_north_to_east_deg);
      corsika_tree->Branch("altitude_deg", &altitude_deg);
      corsika_tree->Branch("core_x_north_m", &core_x_north_m);
      corsika_tree->Branch("core_y_west_m", &core_y_west_m);
      corsika_tree->Branch("array_rotation_deg", &array_rotation_deg);
      corsika_tree->Branch("h_first_int_m", &h_first_int_m);
      corsika_tree->Branch("x_max_g_cm2", &x_max_g_cm2);
      corsika_tree->Branch("h_max_m", &h_max_m);
      corsika_tree->Branch("starting_grammage_g_cm2", &starting_grammage_g_cm2);
      corsika_tree->Branch("ground_gammas", &ground_gammas);
      corsika_tree->Branch("ground_electrons", &ground_electrons);
      corsika_tree->Branch("ground_hadrons", &ground_hadrons);
      corsika_tree->Branch("ground_muons", &ground_muons);
      corsika_tree->Branch("has_simtel_mc_shower", &has_simtel_mc_shower);
      configureRootTreeAutoFlush(corsika_tree.get(),
                                 output_cfg.lact_root_auto_flush_mb);

      observation_tree = std::make_unique<TTree>(
          "observations", "Event-telescope integrated p.e. observations");
      observation_tree->SetDirectory(file.get());
      observation_tree->Branch("event_id", &obs_event_id);
      observation_tree->Branch("telescope_id", &obs_telescope_id);
      observation_tree->Branch("triggered", &triggered);
      observation_tree->Branch("n_pixels_camera", &n_pixels_camera);
      observation_tree->Branch("n_pixels_saved", &n_pixels_saved);
      observation_tree->Branch("pixel_id", &obs_pixel_id);
      observation_tree->Branch("image_pe", &image_pe);
      observation_tree->Branch("image_cherenkov_pe", &image_cherenkov_pe);
      if (detector_cfg.enabled) {
        observation_tree->Branch("image_fired_cherenkov_pe",
                                 &image_fired_cherenkov_pe);
        observation_tree->Branch("image_fired_pe", &image_fired_pe);
      }
      if (output_cfg.lact_root_write_components) {
        observation_tree->Branch("image_nsb_pe", &image_nsb_pe);
      }
      observation_tree->Branch("image_time_mean_ns", &image_time_mean_ns);
      observation_tree->Branch("image_time_rms_ns", &image_time_rms_ns);
      observation_tree->Branch("image_time_peak_ns", &image_time_peak_ns);
      observation_tree->Branch("total_pe", &total_pe);
      observation_tree->Branch("time_first_ns", &time_first_ns);
      observation_tree->Branch("time_mean_ns", &time_mean_ns);
      observation_tree->Branch("time_rms_ns", &time_rms_ns);
      observation_tree->Branch("time_peak_ns", &time_peak_ns);
      observation_tree->Branch("impact_parameter_m", &impact_parameter_m);
      observation_tree->Branch("n_pixels_above_threshold",
                               &n_pixels_above_threshold);
      observation_tree->Branch("trigger_time_ns", &trigger_time_ns);
      observation_tree->Branch("trigger_first_time_ns",
                               &trigger_first_time_ns);
      observation_tree->Branch("trigger_max_multiplicity_time_ns",
                               &trigger_max_multiplicity_time_ns);
      observation_tree->Branch("geometric_delay_ns", &geometric_delay_ns);
      observation_tree->Branch("coincidence_time_ns", &coincidence_time_ns);
      configureRootTreeAutoFlush(observation_tree.get(),
                                 output_cfg.lact_root_auto_flush_mb);

      const bool write_time_series =
          (output_cfg.lact_profile == "timeseries_pe" ||
           output_cfg.lact_profile == "debug_full") &&
          waveform_cfg.enabled && waveform_cfg.source == "pe";
      const bool write_detector_time_series =
          (output_cfg.lact_profile == "timeseries_pe" ||
           output_cfg.lact_profile == "debug_full") &&
          waveform_cfg.enabled && detector_cfg.enabled;
      if (write_time_series || write_detector_time_series) {
        waveform_tree = std::make_unique<TTree>(
            "waveforms", "Sparse p.e. waveform COO rows");
        waveform_tree->SetDirectory(file.get());
        waveform_tree->Branch("event_id", &wf_event_id);
        waveform_tree->Branch("telescope_id", &wf_telescope_id);
        waveform_tree->Branch("n_pixels_camera", &wf_n_pixels_camera);
        waveform_tree->Branch("n_time_bins", &wf_n_time_bins);
        waveform_tree->Branch("pixel_id", &wf_pixel_id);
        waveform_tree->Branch("time_bin", &wf_time_bin);
        // sample_value is the canonical branch. The old "pe" name is kept
        // as a compatibility alias for readers of pre-electronics files.
        waveform_tree->Branch("sample_value", &wf_sample_value);
        waveform_tree->Branch("pe", &wf_sample_value);
        configureRootTreeAutoFlush(waveform_tree.get(),
                                   output_cfg.lact_root_auto_flush_mb);
      }

      if (detector_cfg.enabled && detector_cfg.save_primary_sequence) {
        primary_hit_tree = std::make_unique<TTree>(
            "primary_pe_hits",
            "Detected primary p.e. after optics, collector and PDE");
        primary_hit_tree->SetDirectory(file.get());
        primary_hit_tree->Branch("event_id", &primary_hit_row.event_id);
        primary_hit_tree->Branch("telescope_id",
                                 &primary_hit_row.telescope_id);
        primary_hit_tree->Branch("pixel_id", &primary_hit_row.pixel_id);
        primary_hit_tree->Branch("time_ns", &primary_hit_row.time_ns);
        primary_hit_tree->Branch("sensor_x_m",
                                 &primary_hit_row.sensor_x_m);
        primary_hit_tree->Branch("sensor_y_m",
                                 &primary_hit_row.sensor_y_m);
        primary_hit_tree->Branch("wavelength_nm",
                                 &primary_hit_row.wavelength_nm);
        primary_hit_tree->Branch("primary_pe",
                                 &primary_hit_row.primary_pe);
        primary_hit_tree->Branch("origin", &primary_hit_row.origin);
      }
      if (detector_cfg.enabled && detector_cfg.save_fired_sequence) {
        fired_hit_tree = std::make_unique<TTree>(
            "fired_pe_hits",
            "SiPM microcell firings after explicit saturation");
        fired_hit_tree->SetDirectory(file.get());
        fired_hit_tree->Branch("event_id", &fired_hit_row.event_id);
        fired_hit_tree->Branch("telescope_id",
                               &fired_hit_row.telescope_id);
        fired_hit_tree->Branch("pixel_id", &fired_hit_row.pixel_id);
        fired_hit_tree->Branch("time_ns", &fired_hit_row.time_ns);
        fired_hit_tree->Branch("channel_id", &fired_hit_row.channel_id);
        fired_hit_tree->Branch("microcell_id",
                               &fired_hit_row.microcell_id);
        fired_hit_tree->Branch("fired_pe", &fired_hit_row.fired_pe);
        fired_hit_tree->Branch("origin", &fired_hit_row.origin);
      }
      if (detector_cfg.enabled &&
          detector_cfg.save_microcell_decisions) {
        microcell_decision_tree = std::make_unique<TTree>(
            "microcell_decisions",
            "Explicit no-recovery microcell decisions");
        microcell_decision_tree->SetDirectory(file.get());
        microcell_decision_tree->Branch(
            "event_id", &microcell_decision_row.event_id);
        microcell_decision_tree->Branch(
            "telescope_id", &microcell_decision_row.telescope_id);
        microcell_decision_tree->Branch(
            "pixel_id", &microcell_decision_row.pixel_id);
        microcell_decision_tree->Branch(
            "time_ns", &microcell_decision_row.time_ns);
        microcell_decision_tree->Branch(
            "sensor_x_m", &microcell_decision_row.sensor_x_m);
        microcell_decision_tree->Branch(
            "sensor_y_m", &microcell_decision_row.sensor_y_m);
        microcell_decision_tree->Branch(
            "grid_column", &microcell_decision_row.grid_column);
        microcell_decision_tree->Branch(
            "grid_row", &microcell_decision_row.grid_row);
        microcell_decision_tree->Branch(
            "channel_id", &microcell_decision_row.channel_id);
        microcell_decision_tree->Branch(
            "microcell_id", &microcell_decision_row.microcell_id);
        microcell_decision_tree->Branch(
            "fired", &microcell_decision_row.fired);
        microcell_decision_tree->Branch(
            "origin", &microcell_decision_row.origin);
      }

      trace_tree = std::make_unique<TTree>("trace_summary",
                                           "Event-telescope trace summary");
      trace_tree->SetDirectory(file.get());
      trace_tree->Branch("event_id", &trace_event_id);
      trace_tree->Branch("telescope_id", &trace_telescope_id);
      trace_tree->Branch("input_bunches", &input_bunches);
      trace_tree->Branch("input_photons", &input_photons);
      trace_tree->Branch("blocked_by_obstruction", &blocked_by_obstruction);
      trace_tree->Branch("blocked_incoming", &blocked_incoming);
      trace_tree->Branch("blocked_reflected", &blocked_reflected);
      trace_tree->Branch("hit_mirror", &hit_mirror);
      trace_tree->Branch("hit_output_plane", &hit_output_plane);
      trace_tree->Branch("hit_camera", &hit_camera);
      trace_tree->Branch("accepted_camera", &accepted_camera);
      trace_tree->Branch("lost_between_pixels", &lost_between_pixels);
      trace_tree->Branch("unique_hit_pixels", &unique_hit_pixels);
      trace_tree->Branch("signal_pe", &signal_pe);
      trace_tree->Branch("time_mean_ns", &trace_time_mean_ns);
      trace_tree->Branch("time_rms_ns", &trace_time_rms_ns);
      configureRootTreeAutoFlush(trace_tree.get(),
                                 output_cfg.lact_root_auto_flush_mb);
    }

    void writeCorsikaEvent(long long event_id) {
      if (!written_events.insert(event_id).second)
        return;
      const auto event_meta =
          outputEventMetadata(static_cast<int>(event_id),
                              source_runtime_cfg.event_id_mode, metadata);
      shower_event_id = event_meta.shower_event;
      array_id = event_meta.array_id;
      array_time_offset_ns = event_meta.array_time_offset_ns;
      area_weight_m2 = event_meta.area_weight_m2;
      has_explicit_area_weight = event_meta.has_explicit_area_weight;
      root_event_id = event_id;
      primary_type = 0;
      energy_gev = event_meta.energy_gev;
      theta_deg = std::numeric_limits<double>::quiet_NaN();
      phi_deg = std::numeric_limits<double>::quiet_NaN();
      azimuth_north_to_east_deg = event_meta.azimuth_north_to_east_deg;
      core_x_north_m = event_meta.core_x_north_m;
      core_y_west_m = event_meta.core_y_west_m;
      array_rotation_deg = std::numeric_limits<double>::quiet_NaN();
      altitude_deg = std::numeric_limits<double>::quiet_NaN();
      h_first_int_m = std::numeric_limits<double>::quiet_NaN();
      x_max_g_cm2 = std::numeric_limits<double>::quiet_NaN();
      h_max_m = std::numeric_limits<double>::quiet_NaN();
      starting_grammage_g_cm2 = std::numeric_limits<double>::quiet_NaN();
      ground_gammas = std::numeric_limits<double>::quiet_NaN();
      ground_electrons = std::numeric_limits<double>::quiet_NaN();
      ground_hadrons = std::numeric_limits<double>::quiet_NaN();
      ground_muons = std::numeric_limits<double>::quiet_NaN();
      has_simtel_mc_shower = false;
      const int target_shower_event = shower_event_id;
      auto event_it =
          std::find_if(metadata.events.begin(), metadata.events.end(),
                       [target_shower_event](const EventIOEventHeader &event) {
                         return event.shower_event_id == target_shower_event;
                       });
      if (event_it != metadata.events.end()) {
        primary_type = event_it->primary_type;
        theta_deg = event_it->theta_deg;
        phi_deg = event_it->phi_deg;
        altitude_deg = std::isfinite(event_it->altitude_deg)
                           ? event_it->altitude_deg
                           : 90.0 - event_it->theta_deg;
        array_rotation_deg = event_it->array_rotation_deg;
        h_first_int_m = event_it->h_first_int_m;
        x_max_g_cm2 = event_it->x_max_g_cm2;
        h_max_m = event_it->h_max_m;
        starting_grammage_g_cm2 = event_it->starting_grammage_g_cm2;
        ground_gammas = event_it->ground_gammas;
        ground_electrons = event_it->ground_electrons;
        ground_hadrons = event_it->ground_hadrons;
        ground_muons = event_it->ground_muons;
        has_simtel_mc_shower = event_it->has_simtel_mc_shower;
      }
      if (corsika_tree->Fill() < 0) {
        throw std::runtime_error("failed to fill ROOT corsika_events tree");
      }
    }

    void writeObservation(const LactRootObservation& obs)
    {
        obs_event_id = obs.event_id;
        obs_telescope_id = obs.telescope_id;
        triggered = obs.triggered;
        n_pixels_camera = obs.n_pixels_camera;
        n_pixels_saved = obs.n_pixels_saved;
        obs_pixel_id = obs.pixel_id;
        image_pe = obs.image_pe;
        image_cherenkov_pe = obs.image_cherenkov_pe;
        image_fired_cherenkov_pe = obs.image_fired_cherenkov_pe;
        image_fired_pe = obs.image_fired_pe;
        image_nsb_pe = obs.image_nsb_pe;
        image_time_mean_ns = obs.image_time_mean_ns;
        image_time_rms_ns = obs.image_time_rms_ns;
        image_time_peak_ns = obs.image_time_peak_ns;
        total_pe = obs.total_pe;
        time_first_ns = obs.time_first_ns;
        time_mean_ns = obs.time_mean_ns;
        time_rms_ns = obs.time_rms_ns;
        time_peak_ns = obs.time_peak_ns;
        impact_parameter_m = obs.impact_parameter_m;
        n_pixels_above_threshold = obs.n_pixels_above_threshold;
        trigger_time_ns = obs.trigger_time_ns;
        trigger_first_time_ns = obs.trigger_first_time_ns;
        trigger_max_multiplicity_time_ns =
            obs.trigger_max_multiplicity_time_ns;
        geometric_delay_ns = obs.geometric_delay_ns;
        coincidence_time_ns = obs.coincidence_time_ns;
        if (observation_tree->Fill() < 0) {
            throw std::runtime_error("failed to fill ROOT observations tree");
        }
    }

    void writeWaveformConfig(const LactRootPreparedData& prepared)
    {
        if (waveform_config_written || prepared.time_centers_ns.empty()) return;
        waveform_config_tree = std::make_unique<TTree>("waveform_config", "p.e. waveform metadata");
        waveform_config_tree->SetDirectory(nullptr);
        waveform_enabled = true;
        waveform_source = waveform_cfg.source;
        waveform_sample_unit =
            detector_cfg.enabled
                ? (detector_cfg.single_pe.enabled
                       ? (detector_cfg.single_pe.unit == "mv"
                              ? "mV"
                              : "pe_charge_per_sample")
                       : "fired_pe_per_sample")
                : "pe_per_sample";
        time_reference = waveform_cfg.time_reference;
        time_bin_width_ns = detector_cfg.enabled
            ? detector_cfg.sampling.width_ns
            : waveform_cfg.time_bin_width_ns;
        time_window_start_ns = detector_cfg.enabled
            ? detector_cfg.sampling.start_ns
            : waveform_cfg.time_window_start_ns;
        time_window_end_ns = detector_cfg.enabled
            ? detector_cfg.sampling.end_ns
            : waveform_cfg.time_window_end_ns;
        n_time_bins = static_cast<int>(prepared.time_centers_ns.size());
        time_edges_ns = prepared.time_edges_ns;
        time_centers_ns = prepared.time_centers_ns;
        waveform_config_tree->Branch("waveform_enabled", &waveform_enabled);
        waveform_config_tree->Branch("waveform_source", &waveform_source);
        waveform_config_tree->Branch("sample_unit", &waveform_sample_unit);
        waveform_config_tree->Branch("time_reference", &time_reference);
        waveform_config_tree->Branch("time_bin_width_ns", &time_bin_width_ns);
        waveform_config_tree->Branch("time_window_start_ns", &time_window_start_ns);
        waveform_config_tree->Branch("time_window_end_ns", &time_window_end_ns);
        waveform_config_tree->Branch("n_time_bins", &n_time_bins);
        waveform_config_tree->Branch("time_edges_ns", &time_edges_ns);
        waveform_config_tree->Branch("time_centers_ns", &time_centers_ns);
        waveform_config_tree->Fill();
        waveform_config_tree->Write();
        waveform_config_written = true;
    }

    void writeWaveform(const LactRootWaveform& wf)
    {
        if (!waveform_tree) return;
        wf_event_id = wf.event_id;
        wf_telescope_id = wf.telescope_id;
        wf_n_pixels_camera = wf.n_pixels_camera;
        wf_n_time_bins = wf.n_time_bins;
        wf_pixel_id = wf.pixel_id;
        wf_time_bin = wf.time_bin;
        wf_sample_value = wf.sample_value;
        if (waveform_tree->Fill() < 0) {
            throw std::runtime_error("failed to fill ROOT waveforms tree");
        }
    }

    void writePrimaryHit(const LactRootPrimaryHit& hit)
    {
        if (!primary_hit_tree) return;
        primary_hit_row = hit;
        if (primary_hit_tree->Fill() < 0) {
            throw std::runtime_error(
                "failed to fill ROOT primary_pe_hits tree");
        }
    }

    void writeFiredHit(const LactRootFiredHit& hit)
    {
        if (!fired_hit_tree) return;
        fired_hit_row = hit;
        if (fired_hit_tree->Fill() < 0) {
            throw std::runtime_error(
                "failed to fill ROOT fired_pe_hits tree");
        }
    }

    void writeMicrocellDecision(
        const LactRootMicrocellDecision& item)
    {
        if (!microcell_decision_tree) return;
        microcell_decision_row = item;
        if (microcell_decision_tree->Fill() < 0) {
            throw std::runtime_error(
                "failed to fill ROOT microcell_decisions tree");
        }
    }

    void
    writeTraceSummary(const std::map<SummaryKey, TraceSummary> &summaries) {
      for (const auto &kv : summaries) {
        const auto &s = kv.second;
        trace_event_id = s.event_id;
        trace_telescope_id = s.telescope_id;
        input_bunches = s.input_bunches;
        input_photons = s.input_photons;
        blocked_by_obstruction = s.blocked_by_obstruction;
        blocked_incoming = s.blocked_incoming;
        blocked_reflected = s.blocked_reflected;
        hit_mirror = s.hit_mirror;
        hit_output_plane = s.hit_output_plane;
        hit_camera = s.hit_camera;
        accepted_camera = s.accepted_camera;
        lost_between_pixels = s.lost_between_pixels;
        unique_hit_pixels = static_cast<int>(s.unique_pixels.size());
        signal_pe = s.weighted_signal;
        trace_time_mean_ns = std::numeric_limits<double>::quiet_NaN();
        trace_time_rms_ns = std::numeric_limits<double>::quiet_NaN();
        if (s.weighted_signal > 0.0) {
          trace_time_mean_ns = s.weighted_time_sum / s.weighted_signal;
          trace_time_rms_ns = std::sqrt(
              std::max(0.0, s.weighted_time2_sum / s.weighted_signal -
                                trace_time_mean_ns * trace_time_mean_ns));
        }
        if (trace_tree->Fill() < 0) {
          throw std::runtime_error("failed to fill ROOT trace_summary tree");
        }
      }
    }

    void flushBufferedTrees()
    {
        auto flush_tree = [](TTree* tree) {
            if (!tree) return;
            if (tree->FlushBaskets() < 0) {
                throw std::runtime_error(
                    std::string("failed to flush ROOT tree ") + tree->GetName());
            }
            if (tree->AutoSave("SaveSelf") < 0) {
                throw std::runtime_error(
                    std::string("failed to autosave ROOT tree ") + tree->GetName());
            }
        };
        flush_tree(corsika_tree.get());
        flush_tree(observation_tree.get());
        flush_tree(waveform_tree.get());
        flush_tree(primary_hit_tree.get());
        flush_tree(fired_hit_tree.get());
        flush_tree(microcell_decision_tree.get());
        flush_tree(trace_tree.get());
        if (file) file->Flush();
    }

    void writeEvent(const std::map<SummaryKey, TraceSummary>& summaries,
                    const std::map<PixelKey, PixelAccumulator>& pixels,
                    const std::map<WaveformKey, WaveformPixelAccumulator>& waveforms,
                    const std::vector<RawWaveformHit>& raw_waveform_hits)
    {
        LactRootPreparedData prepared = prepareLactRootObservations(
            output_cfg, waveform_cfg, detector_cfg, nsb_cfg, trigger_cfg,
            source_runtime_cfg, telescope_cfg, metadata, camera,
            summaries, pixels, waveforms, raw_waveform_hits);
        writeWaveformConfig(prepared);
        for (const auto& obs : prepared.observations) {
            writeCorsikaEvent(obs.event_id);
            writeObservation(obs);
        }
        for (const auto& wf : prepared.waveforms) {
            writeWaveform(wf);
        }
        for (const auto& hit : prepared.primary_hits) {
            writePrimaryHit(hit);
        }
        for (const auto& hit : prepared.fired_hits) {
            writeFiredHit(hit);
        }
        for (const auto& item : prepared.microcell_decisions) {
            writeMicrocellDecision(item);
        }
        writeTraceSummary(summaries);
        if (output_cfg.lact_root_flush_events > 0 &&
            ++events_since_flush >= output_cfg.lact_root_flush_events) {
            flushBufferedTrees();
            events_since_flush = 0;
        }
    }

    void finish()
    {
        if (finished) return;
        file->cd();
        flushBufferedTrees();
        if (file->Write("", TObject::kOverwrite) < 0) {
            throw std::runtime_error("failed to finalize lact_event ROOT file");
        }
        // Streaming trees must be attached while filling so ROOT can flush
        // their baskets.  Detach them before closing the file because the
        // unique_ptr members, rather than TFile, own their lifetime.
        if (corsika_tree) corsika_tree->SetDirectory(nullptr);
        if (observation_tree) observation_tree->SetDirectory(nullptr);
        if (waveform_tree) waveform_tree->SetDirectory(nullptr);
        if (primary_hit_tree) primary_hit_tree->SetDirectory(nullptr);
        if (fired_hit_tree) fired_hit_tree->SetDirectory(nullptr);
        if (microcell_decision_tree) {
            microcell_decision_tree->SetDirectory(nullptr);
        }
        if (trace_tree) trace_tree->SetDirectory(nullptr);
        file->Close();
        finished = true;
    }
};

LactEventRootStreamWriter::LactEventRootStreamWriter(
    const CorsikaTraceOutputConfig& output_cfg,
    const WaveformOutputConfig& waveform_cfg,
    const std::string& main_config_path,
    const std::map<std::string, std::string>& cfg,
    const SourceRuntimeConfig& source_runtime_cfg,
    const TelescopeConfig& telescope_cfg,
    const EventIOMetadata& metadata,
    const CameraGeometry& camera,
    const std::vector<MirrorFacet>& facets,
    const NsbConfig& nsb_cfg,
    const TriggerConfig& trigger_cfg)
    : impl_(std::make_unique<Impl>(output_cfg, waveform_cfg, main_config_path, cfg,
                                   source_runtime_cfg, telescope_cfg, metadata,
                                   camera, facets, nsb_cfg, trigger_cfg))
{
}

LactEventRootStreamWriter::~LactEventRootStreamWriter()
{
    if (impl_) {
        impl_->finish();
    }
}

void LactEventRootStreamWriter::writeEvent(
    const std::map<SummaryKey, TraceSummary>& summaries,
    const std::map<PixelKey, PixelAccumulator>& pixels,
    const std::map<WaveformKey, WaveformPixelAccumulator>& waveforms,
    const std::vector<RawWaveformHit>& raw_waveform_hits)
{
    impl_->writeEvent(summaries, pixels, waveforms, raw_waveform_hits);
}

void LactEventRootStreamWriter::finish()
{
    impl_->finish();
}

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
                        const std::vector<RawWaveformHit>& raw_waveform_hits)
{
    LactEventRootStreamWriter writer(output_cfg, waveform_cfg, main_config_path, cfg,
                                     source_runtime_cfg, telescope_cfg, metadata, camera,
                                     facets, nsb_cfg, trigger_cfg);
    writer.writeEvent(summaries, pixels, waveforms, raw_waveform_hits);
    writer.finish();
}

} // namespace lact
