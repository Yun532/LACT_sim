#include "io/LactEventRootWriter.hpp"

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
    out.shower_event = showerEventFromOutputEvent(event_id, event_id_mode);
    out.array_id = arrayIdFromOutputEvent(event_id, event_id_mode);

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
        const std::size_t offset_index = static_cast<std::size_t>(out.array_id);
        if (out.array_id >= 0 && offset_index < offsets->x_m.size() &&
            offset_index < offsets->y_m.size()) {
            out.core_x_north_m = -offsets->x_m[offset_index];
            out.core_y_west_m = -offsets->y_m[offset_index];
            out.used_array_offset = true;
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
};

struct LactRootWaveform {
    long long event_id = 0;
    int telescope_id = 0;
    int n_pixels_camera = 0;
    int n_time_bins = 0;
    std::vector<int> pixel_id;
    std::vector<unsigned short> time_bin;
    std::vector<float> pe;
};

struct LactRootPreparedData {
    std::vector<int> pixel_axis;
    std::vector<double> time_edges_ns;
    std::vector<double> time_centers_ns;
    std::vector<LactRootObservation> observations;
    std::vector<LactRootWaveform> waveforms;
};

LactRootPreparedData prepareLactRootObservations(
    const CorsikaTraceOutputConfig& output_cfg,
    const WaveformOutputConfig& waveform_cfg,
    const NsbConfig& nsb_cfg,
    const TriggerConfig& trigger_cfg,
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

    const bool write_time_series =
        (output_cfg.lact_profile == "timeseries_pe" ||
         output_cfg.lact_profile == "debug_full") &&
        waveform_cfg.enabled && waveform_cfg.source == "pe" && !prepared.pixel_axis.empty();
    const std::size_t n_pixels = prepared.pixel_axis.size();
    const std::size_t n_bins = write_time_series ? waveformBinCount(waveform_cfg) : 0;
    if (write_time_series) {
        prepared.time_edges_ns.resize(n_bins + 1);
        prepared.time_centers_ns.resize(n_bins);
        for (std::size_t i = 0; i <= n_bins; ++i) {
            prepared.time_edges_ns[i] =
                waveform_cfg.time_window_start_ns +
                static_cast<double>(i) * waveform_cfg.time_bin_width_ns;
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
    std::map<int, int> triggered_telescopes_by_event;
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
        std::vector<double> image_nsb_pe_by_col(n_pixels, 0.0);
        std::vector<double> time_sum_by_col(n_pixels, 0.0);
        std::vector<double> time2_sum_by_col(n_pixels, 0.0);
        std::unordered_map<std::size_t, double> waveform_pe;
        std::vector<double> waveform_peak_by_col;
        std::vector<std::size_t> waveform_peak_bin_by_col;
        double reference_time_ns = 0.0;
        double trigger_time_ns = std::numeric_limits<double>::quiet_NaN();

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

        if (write_time_series) {
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

            const std::size_t trigger_window_bins =
                trigger_cfg.coincidence_window_ns > 0.0
                    ? std::min(n_bins, std::max<std::size_t>(
                                      1,
                                      static_cast<std::size_t>(std::ceil(
                                          trigger_cfg.coincidence_window_ns /
                                          waveform_cfg.time_bin_width_ns))))
                    : n_bins;
            const auto waveform_pe_at = [&](std::size_t col, std::size_t bin) {
                const auto it = waveform_pe.find(col * n_bins + bin);
                return it == waveform_pe.end() ? 0.0 : it->second;
            };
            const auto find_best_trigger_window =
                [&](std::size_t first_window,
                    std::size_t last_window) -> std::pair<int, std::size_t> {
                    if (n_bins == 0 || first_window >= n_bins) {
                        return {0, 0};
                    }
                    last_window = std::min(last_window, n_bins - 1);
                    std::vector<int> pixels_above_threshold_by_window(n_bins, 0);
                    for (std::size_t col = 0; col < n_pixels; ++col) {
                        double window_pe = 0.0;
                        const std::size_t first_end =
                            std::min(n_bins, first_window + trigger_window_bins);
                        for (std::size_t bin = first_window; bin < first_end; ++bin) {
                            window_pe += waveform_pe_at(col, bin);
                        }
                        if (window_pe >= trigger_cfg.pixel_threshold_pe) {
                            ++pixels_above_threshold_by_window[first_window];
                        }
                        for (std::size_t window_start = first_window + 1;
                             window_start <= last_window;
                             ++window_start) {
                            window_pe -= waveform_pe_at(col, window_start - 1);
                            const std::size_t add_bin =
                                window_start + trigger_window_bins - 1;
                            if (add_bin < n_bins) {
                                window_pe += waveform_pe_at(col, add_bin);
                            }
                            if (window_pe >= trigger_cfg.pixel_threshold_pe) {
                                ++pixels_above_threshold_by_window[window_start];
                            }
                        }
                    }

                    int best_pixels_above_threshold = 0;
                    std::size_t best_trigger_window_start = first_window;
                    for (std::size_t window_start = first_window;
                         window_start <= last_window;
                         ++window_start) {
                        const int pixels_above_threshold =
                            pixels_above_threshold_by_window[window_start];
                        if (pixels_above_threshold > best_pixels_above_threshold) {
                            best_pixels_above_threshold = pixels_above_threshold;
                            best_trigger_window_start = window_start;
                        }
                    }
                    return {best_pixels_above_threshold, best_trigger_window_start};
                };

            double camera_peak = -1.0;
            std::size_t camera_peak_bin = 0;
            for (std::size_t bin = 0; bin < n_bins; ++bin) {
                if (camera_time_series[bin] > camera_peak) {
                    camera_peak = camera_time_series[bin];
                    camera_peak_bin = bin;
                }
            }

            const auto cherenkov_trigger =
                find_best_trigger_window(0, n_bins - 1);
            bool should_add_full_nsb =
                !output_cfg.save_only_triggered || !trigger_cfg.enabled ||
                cherenkov_trigger.first >= trigger_cfg.camera_multiplicity;
            if (!should_add_full_nsb && nsb_cfg.enabled && camera_peak > 0.0) {
                const std::size_t first_window =
                    camera_peak_bin >= trigger_window_bins
                        ? camera_peak_bin - trigger_window_bins + 1
                        : 0;
                const std::size_t last_window = std::min(camera_peak_bin, n_bins - 1);
                const std::size_t first_light_bin = first_window;
                const std::size_t last_light_bin = std::min(
                    n_bins - 1,
                    last_window + trigger_window_bins - 1);
                const std::size_t light_bin_count =
                    last_light_bin - first_light_bin + 1;
                std::vector<float> light_nsb_pe(n_pixels * light_bin_count, 0.0f);
                for (std::size_t col = 0; col < n_pixels; ++col) {
                    for (std::size_t offset = 0; offset < light_bin_count; ++offset) {
                        const std::size_t bin = first_light_bin + offset;
                        light_nsb_pe[col * light_bin_count + offset] =
                            sampleTimeBinnedNsbPeCell(
                                nsb_cfg,
                                waveform_cfg,
                                event_id,
                                telescope_id,
                                n_pixels,
                                n_bins,
                                col,
                                bin);
                    }
                }
                const auto light_nsb_at = [&](std::size_t col, std::size_t bin) {
                    if (bin < first_light_bin || bin > last_light_bin) {
                        return 0.0f;
                    }
                    return light_nsb_pe[col * light_bin_count + (bin - first_light_bin)];
                };
                int best_pixels_with_light_nsb = 0;
                std::size_t best_light_window = first_window;
                for (std::size_t window_start = first_window;
                     window_start <= last_window;
                     ++window_start) {
                    const std::size_t window_end =
                        std::min(n_bins, window_start + trigger_window_bins);
                    int pixels_above_threshold = 0;
                    for (std::size_t col = 0; col < n_pixels; ++col) {
                        double window_pe = 0.0;
                        for (std::size_t bin = window_start; bin < window_end; ++bin) {
                            window_pe += waveform_pe_at(col, bin);
                            window_pe += light_nsb_at(col, bin);
                        }
                        if (window_pe >= trigger_cfg.pixel_threshold_pe) {
                            ++pixels_above_threshold;
                        }
                    }
                    if (pixels_above_threshold > best_pixels_with_light_nsb) {
                        best_pixels_with_light_nsb = pixels_above_threshold;
                        best_light_window = window_start;
                    }
                }
                should_add_full_nsb =
                    best_pixels_with_light_nsb >= trigger_cfg.camera_multiplicity;
                if (should_add_full_nsb) {
                    trigger_time_ns =
                        reference_time_ns +
                        prepared.time_centers_ns[std::min(
                            n_bins - 1,
                            best_light_window + trigger_window_bins / 2)];
                }
            }

            if (should_add_full_nsb) {
                for (std::size_t col = 0; col < n_pixels; ++col) {
                    for (std::size_t bin = 0; bin < n_bins; ++bin) {
                        const float nsb_pe = sampleTimeBinnedNsbPeCell(
                            nsb_cfg,
                            waveform_cfg,
                            event_id,
                            telescope_id,
                            n_pixels,
                            n_bins,
                            col,
                            bin);
                        if (nsb_pe > 0.0f) {
                            add_waveform_pe(col, bin, nsb_pe, 0.0, nsb_pe);
                        }
                    }
                }
            }

            waveform_peak_by_col.assign(n_pixels, -1.0);
            waveform_peak_bin_by_col.assign(n_pixels, 0);
            for (const auto& entry : waveform_pe) {
                const std::size_t col = entry.first / n_bins;
                const std::size_t bin = entry.first % n_bins;
                const double pe = entry.second;
                if (pe > waveform_peak_by_col[col]) {
                    waveform_peak_by_col[col] = pe;
                    waveform_peak_bin_by_col[col] = bin;
                }
            }
            for (std::size_t col = 0; col < n_pixels; ++col) {
                const double total = image_pe_by_col[col];
                const double peak = waveform_peak_by_col[col];
                const std::size_t peak_bin = waveform_peak_bin_by_col[col];
                if (total > 0.0 && peak > 0.0) {
                    const double peak_time =
                        reference_time_ns + prepared.time_centers_ns[peak_bin];
                    time_sum_by_col[col] = total * peak_time;
                    time2_sum_by_col[col] = total * peak_time * peak_time;
                }
            }

            camera_peak = -1.0;
            camera_peak_bin = 0;
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
            const auto final_trigger = find_best_trigger_window(0, n_bins - 1);
            obs.n_pixels_above_threshold = final_trigger.first;
            if (!std::isfinite(trigger_time_ns) && final_trigger.first > 0) {
                const std::size_t trigger_bin = std::min(
                    n_bins - 1,
                    final_trigger.second + trigger_window_bins / 2);
                trigger_time_ns = reference_time_ns + prepared.time_centers_ns[trigger_bin];
            }
        } else {
            const PixelKey begin_key{event_id, telescope_id, std::numeric_limits<int>::min()};
            const PixelKey end_key{event_id, telescope_id, std::numeric_limits<int>::max()};
            for (auto it = pixels.lower_bound(begin_key);
                 it != pixels.end() && it->first <= end_key;
                 ++it) {
                const auto& p = it->second;
                const auto col_it = pixel_to_col.find(p.pixel_id);
                if (col_it == pixel_to_col.end()) continue;
                const std::size_t col = col_it->second;
                image_pe_by_col[col] = p.pe;
                image_cherenkov_pe_by_col[col] = p.pe;
                time_sum_by_col[col] = p.time_sum;
                time2_sum_by_col[col] = p.time2_sum;
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
            if (pe <= 0.0) continue;
            obs.pixel_id.push_back(prepared.pixel_axis[col]);
            obs.image_pe.push_back(static_cast<float>(pe));
            obs.image_cherenkov_pe.push_back(
                static_cast<float>(image_cherenkov_pe_by_col[col]));
            if (output_cfg.lact_root_write_components) {
                obs.image_nsb_pe.push_back(static_cast<float>(image_nsb_pe_by_col[col]));
            }
            const double mean = time_sum_by_col[col] > 0.0
                ? time_sum_by_col[col] / pe
                : std::numeric_limits<double>::quiet_NaN();
            double rms = std::numeric_limits<double>::quiet_NaN();
            if (time2_sum_by_col[col] > 0.0 && std::isfinite(mean)) {
                const double var = std::max(0.0, time2_sum_by_col[col] / pe - mean * mean);
                rms = std::sqrt(var);
            }
            obs.image_time_mean_ns.push_back(static_cast<float>(mean));
            obs.image_time_rms_ns.push_back(static_cast<float>(rms));
            if (write_time_series && col < waveform_peak_by_col.size()) {
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
            if (!write_time_series && pe >= trigger_cfg.pixel_threshold_pe) {
                ++obs.n_pixels_above_threshold;
            }
        }

        obs.n_pixels_saved = static_cast<int>(obs.pixel_id.size());
        obs.triggered = trigger_cfg.enabled &&
            obs.n_pixels_above_threshold >= trigger_cfg.camera_multiplicity;
        if (obs.triggered) {
            obs.trigger_time_ns =
                std::isfinite(trigger_time_ns) ? trigger_time_ns :
                (std::isfinite(obs.time_peak_ns) ? obs.time_peak_ns : obs.time_mean_ns);
            triggered_telescopes_by_event[event_id] += 1;
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
                wf.pe.push_back(static_cast<float>(pe));
            }
            candidate.waveform = std::move(wf);
            candidate.has_waveform = true;
        }

        candidate.observation = std::move(obs);
        candidates.push_back(std::move(candidate));
    }

    for (auto& candidate : candidates) {
        if (output_cfg.save_only_triggered && trigger_cfg.enabled) {
            const int n_triggered =
                triggered_telescopes_by_event[candidate.observation.event_id];
            if (n_triggered < trigger_cfg.array_multiplicity ||
                !candidate.observation.triggered) {
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
    NsbConfig nsb_cfg;
    TriggerConfig trigger_cfg;
    std::unique_ptr<TFile> file;
    std::unique_ptr<TTree> corsika_tree;
    std::unique_ptr<TTree> observation_tree;
    std::unique_ptr<TTree> waveform_config_tree;
    std::unique_ptr<TTree> waveform_tree;
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
    std::vector<float> image_nsb_pe;
    std::vector<float> image_time_mean_ns;
    std::vector<float> image_time_rms_ns;
    std::vector<float> image_time_peak_ns;
    double total_pe = 0.0, time_first_ns = 0.0, time_mean_ns = 0.0;
    double time_rms_ns = 0.0, time_peak_ns = 0.0, impact_parameter_m = 0.0;
    int n_pixels_above_threshold = 0;
    double trigger_time_ns = 0.0;

    bool waveform_enabled = true;
    std::string waveform_source;
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
    std::vector<float> wf_pe;

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
          camera(camera_in), nsb_cfg(nsb_cfg_in), trigger_cfg(trigger_cfg_in) {
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
      std::string producer_version = "unknown";
      std::string source_kind = "EventIO";
      std::string source_path = source_runtime_cfg.eventio_path;
      std::string source_sha256;
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
      corsika_tree->SetDirectory(nullptr);
      corsika_tree->Branch("event_id", &root_event_id);
      corsika_tree->Branch("shower_event_id", &shower_event_id);
      corsika_tree->Branch("array_id", &array_id);
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
      observation_tree->SetDirectory(nullptr);
      observation_tree->Branch("event_id", &obs_event_id);
      observation_tree->Branch("telescope_id", &obs_telescope_id);
      observation_tree->Branch("triggered", &triggered);
      observation_tree->Branch("n_pixels_camera", &n_pixels_camera);
      observation_tree->Branch("n_pixels_saved", &n_pixels_saved);
      observation_tree->Branch("pixel_id", &obs_pixel_id);
      observation_tree->Branch("image_pe", &image_pe);
      observation_tree->Branch("image_cherenkov_pe", &image_cherenkov_pe);
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
      configureRootTreeAutoFlush(observation_tree.get(),
                                 output_cfg.lact_root_auto_flush_mb);

      const bool write_time_series =
          (output_cfg.lact_profile == "timeseries_pe" ||
           output_cfg.lact_profile == "debug_full") &&
          waveform_cfg.enabled && waveform_cfg.source == "pe";
      if (write_time_series) {
        waveform_tree = std::make_unique<TTree>(
            "waveforms", "Sparse p.e. waveform COO rows");
        waveform_tree->SetDirectory(nullptr);
        waveform_tree->Branch("event_id", &wf_event_id);
        waveform_tree->Branch("telescope_id", &wf_telescope_id);
        waveform_tree->Branch("n_pixels_camera", &wf_n_pixels_camera);
        waveform_tree->Branch("n_time_bins", &wf_n_time_bins);
        waveform_tree->Branch("pixel_id", &wf_pixel_id);
        waveform_tree->Branch("time_bin", &wf_time_bin);
        waveform_tree->Branch("pe", &wf_pe);
        configureRootTreeAutoFlush(waveform_tree.get(),
                                   output_cfg.lact_root_auto_flush_mb);
      }

      trace_tree = std::make_unique<TTree>("trace_summary",
                                           "Event-telescope trace summary");
      trace_tree->SetDirectory(nullptr);
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
      corsika_tree->Fill();
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
        observation_tree->Fill();
    }

    void writeWaveformConfig(const LactRootPreparedData& prepared)
    {
        if (waveform_config_written || prepared.time_centers_ns.empty()) return;
        waveform_config_tree = std::make_unique<TTree>("waveform_config", "p.e. waveform metadata");
        waveform_config_tree->SetDirectory(nullptr);
        waveform_enabled = true;
        waveform_source = waveform_cfg.source;
        time_reference = waveform_cfg.time_reference;
        time_bin_width_ns = waveform_cfg.time_bin_width_ns;
        time_window_start_ns = waveform_cfg.time_window_start_ns;
        time_window_end_ns = waveform_cfg.time_window_end_ns;
        n_time_bins = static_cast<int>(prepared.time_centers_ns.size());
        time_edges_ns = prepared.time_edges_ns;
        time_centers_ns = prepared.time_centers_ns;
        waveform_config_tree->Branch("waveform_enabled", &waveform_enabled);
        waveform_config_tree->Branch("waveform_source", &waveform_source);
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
        wf_pe = wf.pe;
        waveform_tree->Fill();
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
        trace_tree->Fill();
      }
    }

    void flushBufferedTrees()
    {
        if (corsika_tree) corsika_tree->FlushBaskets();
        if (observation_tree) observation_tree->FlushBaskets();
        if (waveform_tree) waveform_tree->FlushBaskets();
        if (trace_tree) trace_tree->FlushBaskets();
        if (file) file->Flush();
    }

    void writeEvent(const std::map<SummaryKey, TraceSummary>& summaries,
                    const std::map<PixelKey, PixelAccumulator>& pixels,
                    const std::map<WaveformKey, WaveformPixelAccumulator>& waveforms,
                    const std::vector<RawWaveformHit>& raw_waveform_hits)
    {
        LactRootPreparedData prepared = prepareLactRootObservations(
            output_cfg, waveform_cfg, nsb_cfg, trigger_cfg, camera,
            summaries, pixels, waveforms, raw_waveform_hits);
        writeWaveformConfig(prepared);
        for (const auto& obs : prepared.observations) {
            writeCorsikaEvent(obs.event_id);
            writeObservation(obs);
        }
        for (const auto& wf : prepared.waveforms) {
            writeWaveform(wf);
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
        if (corsika_tree) corsika_tree->Write();
        if (observation_tree) observation_tree->Write();
        if (waveform_tree) waveform_tree->Write();
        if (trace_tree) trace_tree->Write();
        file->Write();
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
