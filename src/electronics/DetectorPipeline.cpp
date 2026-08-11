#include "electronics/DetectorPipeline.hpp"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace lact::electronics {
namespace {

struct PulsePoint {
    double time_ns = 0.0;
    double amplitude = 0.0;
};

std::string lower(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

void addPrimary(PixelSummary& pixel, HitOrigin origin, double pe)
{
    if (origin == HitOrigin::Cherenkov) pixel.primary_cherenkov_pe += pe;
    else if (origin == HitOrigin::Nsb) pixel.primary_nsb_pe += pe;
    else pixel.primary_dark_pe += pe;
}

void addFired(PixelSummary& pixel, HitOrigin origin, double pe)
{
    if (origin == HitOrigin::Cherenkov) pixel.fired_cherenkov_pe += pe;
    else if (origin == HitOrigin::Nsb) pixel.fired_nsb_pe += pe;
    else pixel.fired_dark_pe += pe;
}

double trapezoidArea(const std::vector<PulsePoint>& pulse)
{
    double area = 0.0;
    for (std::size_t i = 1; i < pulse.size(); ++i) {
        area += 0.5 * (pulse[i - 1].amplitude + pulse[i].amplitude) *
                (pulse[i].time_ns - pulse[i - 1].time_ns);
    }
    return area;
}

std::vector<PulsePoint> loadPulseCsv(const std::string& path)
{
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("failed to open single-p.e. CSV: " + path);
    }
    std::vector<PulsePoint> pulse;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::replace(line.begin(), line.end(), ';', ',');
        std::istringstream row(line);
        std::string first;
        std::string second;
        if (!std::getline(row, first, ',') ||
            !std::getline(row, second, ',')) {
            continue;
        }
        try {
            const double time_ns = std::stod(first);
            const double amplitude = std::stod(second);
            if (std::isfinite(time_ns) && std::isfinite(amplitude)) {
                pulse.push_back({time_ns, amplitude});
            }
        } catch (...) {
            // Header row.
        }
    }
    if (pulse.size() < 2) {
        throw std::runtime_error(
            "single-p.e. CSV must contain at least two numeric rows");
    }
    std::sort(pulse.begin(), pulse.end(),
              [](const PulsePoint& a, const PulsePoint& b) {
                  return a.time_ns < b.time_ns;
              });
    const std::size_t baseline_count =
        std::max<std::size_t>(1, std::min<std::size_t>(100, pulse.size() / 5));
    double baseline = 0.0;
    for (std::size_t i = 0; i < baseline_count; ++i) {
        baseline += pulse[i].amplitude;
    }
    baseline /= static_cast<double>(baseline_count);
    for (auto& point : pulse) point.amplitude -= baseline;

    auto peak = std::max_element(
        pulse.begin(), pulse.end(),
        [](const PulsePoint& a, const PulsePoint& b) {
            return std::abs(a.amplitude) < std::abs(b.amplitude);
        });
    if (peak != pulse.end() && peak->amplitude < 0.0) {
        for (auto& point : pulse) point.amplitude = -point.amplitude;
    }
    return pulse;
}

std::vector<PulsePoint> makePulse(const SinglePeConfig& config,
                                  double sample_width_ns)
{
    std::vector<PulsePoint> pulse;
    if (config.model == "measured_csv") {
        pulse = loadPulseCsv(config.csv_path);
        if (config.template_time_reference == "onset") {
            const double onset = pulse.front().time_ns;
            for (auto& point : pulse) point.time_ns -= onset;
        }
    } else if (config.model == "analytic") {
        const double step_ns = std::min(0.05, sample_width_ns / 32.0);
        const std::size_t count = static_cast<std::size_t>(
            std::ceil(config.support_ns / step_ns));
        pulse.reserve(count + 1);
        for (std::size_t i = 0; i <= count; ++i) {
            const double time_ns =
                std::min(config.support_ns, static_cast<double>(i) * step_ns);
            const double amplitude =
                std::exp(-time_ns / config.fall_ns) -
                std::exp(-time_ns / config.rise_ns);
            pulse.push_back({time_ns, amplitude});
        }
    } else if (config.model == "impulse") {
        pulse = {{0.0, 2.0 / sample_width_ns},
                 {0.5 * sample_width_ns, 2.0 / sample_width_ns},
                 {sample_width_ns, 0.0}};
    }
    if (config.model != "measured_csv") {
        pulse.erase(
            std::remove_if(pulse.begin(), pulse.end(),
                           [&](const PulsePoint& point) {
                               return point.time_ns > config.support_ns;
                           }),
            pulse.end());
    }
    if (pulse.size() < 2) {
        throw std::runtime_error(
            "single-p.e. pulse has fewer than two points inside support");
    }
    if (config.model != "measured_csv" &&
        pulse.back().time_ns < config.support_ns) {
        pulse.push_back({config.support_ns, 0.0});
    }
    for (auto& point : pulse) point.amplitude *= config.amplitude_scale;
    return pulse;
}

std::vector<double> cumulativePulse(const std::vector<PulsePoint>& pulse)
{
    std::vector<double> cumulative(pulse.size(), 0.0);
    for (std::size_t i = 1; i < pulse.size(); ++i) {
        cumulative[i] = cumulative[i - 1] +
            0.5 * (pulse[i - 1].amplitude + pulse[i].amplitude) *
                (pulse[i].time_ns - pulse[i - 1].time_ns);
    }
    return cumulative;
}

double cumulativeAt(const std::vector<PulsePoint>& pulse,
                    const std::vector<double>& cumulative,
                    double time_ns)
{
    if (time_ns <= pulse.front().time_ns) return 0.0;
    if (time_ns >= pulse.back().time_ns) return cumulative.back();
    const auto upper = std::upper_bound(
        pulse.begin(), pulse.end(), time_ns,
        [](double value, const PulsePoint& point) {
            return value < point.time_ns;
        });
    const std::size_t right =
        static_cast<std::size_t>(std::distance(pulse.begin(), upper));
    const std::size_t left = right - 1;
    const double fraction =
        (time_ns - pulse[left].time_ns) /
        (pulse[right].time_ns - pulse[left].time_ns);
    const double amplitude =
        pulse[left].amplitude +
        fraction * (pulse[right].amplitude - pulse[left].amplitude);
    const double dt = time_ns - pulse[left].time_ns;
    return cumulative[left] +
        0.5 * (pulse[left].amplitude + amplitude) * dt;
}

int sampleBin(const SamplingConfig& sampling, double time_ns)
{
    if (time_ns < sampling.start_ns || time_ns >= sampling.end_ns) return -1;
    return static_cast<int>(
        std::floor((time_ns - sampling.start_ns) / sampling.width_ns));
}

CameraTriggerDecision evaluatePeCountTrigger(
    const DetectorPipelineConfig& config,
    std::size_t n_pixels,
    std::size_t n_samples,
    const std::vector<FiredCellHit>& fired_hits)
{
    CameraTriggerDecision out;
    out.pixels_above_threshold.assign(n_samples, 0);
    if (!config.camera_trigger.enabled || n_samples == 0) return out;

    std::vector<double> counts(n_pixels * n_samples, 0.0);
    for (const auto& hit : fired_hits) {
        const int bin = sampleBin(config.sampling, hit.time_ns);
        if (bin < 0 || hit.pixel_id < 0 ||
            static_cast<std::size_t>(hit.pixel_id) >= n_pixels) {
            continue;
        }
        counts[static_cast<std::size_t>(hit.pixel_id) * n_samples +
               static_cast<std::size_t>(bin)] += hit.fired_pe;
    }
    const std::size_t window_bins = std::max<std::size_t>(
        1, static_cast<std::size_t>(std::ceil(
               config.camera_trigger.coincidence_window_ns /
               config.sampling.width_ns)));
    for (std::size_t pixel = 0; pixel < n_pixels; ++pixel) {
        for (std::size_t end = 0; end < n_samples; ++end) {
            double pe = 0.0;
            const std::size_t start =
                end + 1 > window_bins ? end + 1 - window_bins : 0;
            for (std::size_t bin = start; bin <= end; ++bin) {
                pe += counts[pixel * n_samples + bin];
            }
            if (pe >= config.camera_trigger.pixel_threshold_pe) {
                ++out.pixels_above_threshold[end];
            }
        }
    }
    for (std::size_t bin = 0; bin < n_samples; ++bin) {
        out.max_pixels_above_threshold =
            std::max(out.max_pixels_above_threshold,
                     out.pixels_above_threshold[bin]);
        if (out.first_trigger_bin < 0 &&
            out.pixels_above_threshold[bin] >=
                config.camera_trigger.multiplicity) {
            out.first_trigger_bin = static_cast<int>(bin);
        }
    }
    out.triggered = out.first_trigger_bin >= 0;
    if (out.triggered) {
        out.trigger_time_ns =
            config.sampling.start_ns +
            (static_cast<double>(out.first_trigger_bin) + 0.5) *
                config.sampling.width_ns;
    }
    return out;
}

CameraTriggerDecision evaluateVoltageTrigger(
    const DetectorPipelineConfig& config,
    std::size_t n_pixels,
    std::size_t n_samples,
    const std::vector<double>& waveform)
{
    CameraTriggerDecision out;
    out.pixels_above_threshold.assign(n_samples, 0);
    if (!config.camera_trigger.enabled || n_samples == 0) return out;
    const std::size_t window_bins = std::max<std::size_t>(
        1, static_cast<std::size_t>(std::ceil(
               config.camera_trigger.coincidence_window_ns /
               config.sampling.width_ns)));
    for (std::size_t pixel = 0; pixel < n_pixels; ++pixel) {
        for (std::size_t end = 0; end < n_samples; ++end) {
            const std::size_t start =
                end + 1 > window_bins ? end + 1 - window_bins : 0;
            bool crossed = false;
            for (std::size_t bin = start; bin <= end; ++bin) {
                if (waveform[pixel * n_samples + bin] >=
                    config.camera_trigger.pixel_threshold_mv) {
                    crossed = true;
                    break;
                }
            }
            if (crossed) ++out.pixels_above_threshold[end];
        }
    }
    for (std::size_t bin = 0; bin < n_samples; ++bin) {
        out.max_pixels_above_threshold =
            std::max(out.max_pixels_above_threshold,
                     out.pixels_above_threshold[bin]);
        if (out.first_trigger_bin < 0 &&
            out.pixels_above_threshold[bin] >=
                config.camera_trigger.multiplicity) {
            out.first_trigger_bin = static_cast<int>(bin);
        }
    }
    out.triggered = out.first_trigger_bin >= 0;
    if (out.triggered) {
        out.trigger_time_ns =
            config.sampling.start_ns +
            (static_cast<double>(out.first_trigger_bin) + 0.5) *
                config.sampling.width_ns;
    }
    return out;
}

} // namespace

double interChannelActiveFraction(const MicrocellConfig& config)
{
    if (config.layout != "s17351_tiled_2x4") {
        return 1.0;
    }
    const double sensor_area =
        config.sensor_size_x_m * config.sensor_size_y_m;
    const double channel_area =
        static_cast<double>(config.channel_columns * config.channel_rows) *
        config.channel_size_x_m * config.channel_size_y_m;
    if (!(sensor_area > 0.0) || !(channel_area > 0.0) ||
        channel_area > sensor_area) {
        throw std::runtime_error(
            "invalid S17351 area for inter-channel gap fraction");
    }
    return channel_area / sensor_area;
}

MicrocellAddress mapMicrocellPosition(
    const MicrocellConfig& config,
    double sensor_x_m,
    double sensor_y_m)
{
    MicrocellAddress out;
    if (!std::isfinite(sensor_x_m) || !std::isfinite(sensor_y_m) ||
        !(config.sensor_size_x_m > 0.0) ||
        !(config.sensor_size_y_m > 0.0)) {
        return out;
    }
    const double half_x = 0.5 * config.sensor_size_x_m;
    const double half_y = 0.5 * config.sensor_size_y_m;
    if (sensor_x_m < -half_x || sensor_x_m > half_x ||
        sensor_y_m < -half_y || sensor_y_m > half_y) {
        return out;
    }
    out.inside_sensor = true;

    if (config.layout == "uniform_interleaved") {
        const double x_fraction = std::clamp(
            (sensor_x_m + half_x) / config.sensor_size_x_m,
            0.0, std::nextafter(1.0, 0.0));
        const double y_fraction = std::clamp(
            (sensor_y_m + half_y) / config.sensor_size_y_m,
            0.0, std::nextafter(1.0, 0.0));
        out.grid_column = static_cast<int>(
            x_fraction * config.grid_columns);
        out.grid_row = static_cast<int>(
            y_fraction * config.grid_rows);
        const std::uint64_t global_cell =
            static_cast<std::uint64_t>(out.grid_row) *
                static_cast<std::uint64_t>(config.grid_columns) +
            static_cast<std::uint64_t>(out.grid_column);
        out.channel_id = static_cast<int>(
            global_cell %
            static_cast<std::uint64_t>(config.channels_per_pixel));
        out.microcell_id = static_cast<int>(
            global_cell /
            static_cast<std::uint64_t>(config.channels_per_pixel));
        out.inside_channel = true;
        return out;
    }

    if (config.layout != "s17351_tiled_2x4") {
        return out;
    }

    const double local_x = sensor_x_m + half_x;
    const double local_y = sensor_y_m + half_y;
    const double tile_pitch_x =
        config.channel_size_x_m + config.channel_gap_x_m;
    const double tile_pitch_y =
        config.channel_size_y_m + config.channel_gap_y_m;
    int channel_column = static_cast<int>(std::floor(local_x / tile_pitch_x));
    int channel_row_from_bottom =
        static_cast<int>(std::floor(local_y / tile_pitch_y));
    channel_column = std::clamp(
        channel_column, 0, config.channel_columns - 1);
    channel_row_from_bottom = std::clamp(
        channel_row_from_bottom, 0, config.channel_rows - 1);
    const double channel_x =
        local_x - static_cast<double>(channel_column) * tile_pitch_x;
    const double channel_y =
        local_y - static_cast<double>(channel_row_from_bottom) * tile_pitch_y;
    if (channel_x < 0.0 || channel_x >= config.channel_size_x_m ||
        channel_y < 0.0 || channel_y >= config.channel_size_y_m) {
        out.channel_gap = true;
        return out;
    }
    out.inside_channel = true;

    const double pitch_x =
        config.channel_size_x_m /
        static_cast<double>(config.microcell_columns_per_channel);
    const double pitch_y =
        config.channel_size_y_m /
        static_cast<double>(config.microcell_rows_per_channel);
    const int cell_column = std::min(
        config.microcell_columns_per_channel - 1,
        static_cast<int>(std::floor(channel_x / pitch_x)));
    const int cell_row = std::min(
        config.microcell_rows_per_channel - 1,
        static_cast<int>(std::floor(channel_y / pitch_y)));
    out.grid_column =
        channel_column * config.microcell_columns_per_channel + cell_column;
    out.grid_row =
        channel_row_from_bottom * config.microcell_rows_per_channel + cell_row;

    // Front-side S17351 drawing: B-1..B-4 are the left column and
    // A-1..A-4 are the right column, numbered from top to bottom.
    const int row_from_top =
        config.channel_rows - 1 - channel_row_from_bottom;
    out.channel_id =
        channel_column == 1 ? row_from_top
                            : config.channel_rows + row_from_top;
    out.microcell_id =
        cell_row * config.microcell_columns_per_channel + cell_column;

    return out;
}

const char* hitOriginName(HitOrigin origin)
{
    if (origin == HitOrigin::Cherenkov) return "cherenkov";
    if (origin == HitOrigin::Nsb) return "nsb";
    return "dark";
}

HitOrigin parseHitOrigin(const std::string& text)
{
    const std::string value = lower(text);
    if (value == "cherenkov" || value == "signal") {
        return HitOrigin::Cherenkov;
    }
    if (value == "nsb") return HitOrigin::Nsb;
    if (value == "dark") return HitOrigin::Dark;
    throw std::runtime_error("hit origin must be cherenkov, nsb, or dark");
}

void validateDetectorPipelineConfig(const DetectorPipelineConfig& config)
{
    if (!(config.sampling.width_ns > 0.0) ||
        !(config.sampling.end_ns > config.sampling.start_ns)) {
        throw std::runtime_error("invalid detector sampling time grid");
    }
    if (config.microcell.model != "explicit_no_recovery") {
        throw std::runtime_error(
            "only microcell.model=explicit_no_recovery is supported");
    }
    if (!(config.microcell.layout == "uniform_interleaved" ||
          config.microcell.layout == "s17351_tiled_2x4")) {
        throw std::runtime_error(
            "microcell.layout must be uniform_interleaved or "
            "s17351_tiled_2x4");
    }
    if (!(config.microcell.sensor_size_x_m > 0.0) ||
        !(config.microcell.sensor_size_y_m > 0.0) ||
        config.microcell.grid_columns <= 0 ||
        config.microcell.grid_rows <= 0 ||
        config.microcell.channels_per_pixel <= 0 ||
        config.microcell.microcells_per_channel <= 0) {
        throw std::runtime_error("invalid microcell geometry");
    }
    const long long grid_cells =
        static_cast<long long>(config.microcell.grid_columns) *
        static_cast<long long>(config.microcell.grid_rows);
    const long long configured_cells =
        static_cast<long long>(config.microcell.channels_per_pixel) *
        static_cast<long long>(config.microcell.microcells_per_channel);
    if (grid_cells != configured_cells) {
        throw std::runtime_error(
            "microcell grid must equal channels_per_pixel * "
            "microcells_per_channel");
    }
    if (config.microcell.layout == "s17351_tiled_2x4") {
        const auto& microcell = config.microcell;
        if (microcell.channel_columns != 2 ||
            microcell.channel_rows != 4 ||
            microcell.channels_per_pixel !=
                microcell.channel_columns * microcell.channel_rows ||
            microcell.microcell_columns_per_channel <= 0 ||
            microcell.microcell_rows_per_channel <= 0 ||
            microcell.microcells_per_channel !=
                microcell.microcell_columns_per_channel *
                    microcell.microcell_rows_per_channel ||
            microcell.grid_columns !=
                microcell.channel_columns *
                    microcell.microcell_columns_per_channel ||
            microcell.grid_rows !=
                microcell.channel_rows *
                    microcell.microcell_rows_per_channel ||
            !(microcell.channel_size_x_m > 0.0) ||
            !(microcell.channel_size_y_m > 0.0) ||
            microcell.channel_gap_x_m < 0.0 ||
            microcell.channel_gap_y_m < 0.0) {
            throw std::runtime_error("invalid S17351 tiled microcell geometry");
        }
        const double expected_x =
            microcell.channel_columns * microcell.channel_size_x_m +
            (microcell.channel_columns - 1) * microcell.channel_gap_x_m;
        const double expected_y =
            microcell.channel_rows * microcell.channel_size_y_m +
            (microcell.channel_rows - 1) * microcell.channel_gap_y_m;
        constexpr double geometry_tolerance_m = 1.0e-9;
        if (std::abs(expected_x - microcell.sensor_size_x_m) >
                geometry_tolerance_m ||
            std::abs(expected_y - microcell.sensor_size_y_m) >
                geometry_tolerance_m) {
            throw std::runtime_error(
                "S17351 sensor size must include channel tiles and gaps");
        }
        const double active_fraction = interChannelActiveFraction(microcell);
        if (!(active_fraction > 0.0) || active_fraction > 1.0) {
            throw std::runtime_error(
                "invalid S17351 inter-channel active fraction");
        }
    }
    if (!(config.single_pe.model == "analytic" ||
          config.single_pe.model == "measured_csv" ||
          config.single_pe.model == "impulse")) {
        throw std::runtime_error(
            "single_pe.model must be analytic, measured_csv, or impulse");
    }
    if (!(config.single_pe.unit == "pe_charge" ||
          config.single_pe.unit == "mv")) {
        throw std::runtime_error("single_pe.unit must be pe_charge or mV");
    }
    if (config.single_pe.model == "measured_csv" &&
        config.single_pe.csv_path.empty()) {
        throw std::runtime_error(
            "single_pe.csv_path is required for measured_csv");
    }
    if (!(config.single_pe.template_time_reference == "onset" ||
          config.single_pe.template_time_reference == "peak")) {
        throw std::runtime_error(
            "single_pe.template_time_reference must be onset or peak");
    }
    const auto& charge = config.single_pe.charge_fluctuation;
    if (charge.enabled) {
        if (!config.single_pe.enabled) {
            throw std::runtime_error(
                "single-p.e. charge fluctuation requires single_pe.enabled=true");
        }
        if (charge.model != "empirical") {
            throw std::runtime_error(
                "single_pe.charge_fluctuation.model must be empirical");
        }
        if (charge.empirical_samples.empty()) {
            throw std::runtime_error(
                "single-p.e. empirical charge samples are empty");
        }
        for (const double value : charge.empirical_samples) {
            if (!std::isfinite(value) || value <= 0.0) {
                throw std::runtime_error(
                    "single-p.e. empirical charge samples must be finite and > 0");
            }
        }
    }
    const auto& jitter = config.single_pe.time_jitter;
    if (jitter.enabled) {
        if (!config.single_pe.enabled) {
            throw std::runtime_error(
                "single-p.e. time jitter requires single_pe.enabled=true");
        }
        if (jitter.model != "gaussian") {
            throw std::runtime_error(
                "single_pe.time_jitter.model must be gaussian");
        }
        if (!std::isfinite(jitter.sigma_ns) || jitter.sigma_ns < 0.0) {
            throw std::runtime_error(
                "single_pe.time_jitter.sigma_ns must be finite and >= 0");
        }
    }
    if (!(config.single_pe.rise_ns > 0.0) ||
        !(config.single_pe.fall_ns > config.single_pe.rise_ns) ||
        !(config.single_pe.support_ns > 0.0)) {
        throw std::runtime_error("invalid single-p.e. pulse parameters");
    }
    if (!(config.camera_trigger.mode == "pe_count" ||
          config.camera_trigger.mode == "voltage")) {
        throw std::runtime_error(
            "camera_trigger.mode must be pe_count or voltage");
    }
    if (config.camera_trigger.multiplicity <= 0 ||
        config.camera_trigger.coincidence_window_ns < 0.0) {
        throw std::runtime_error("invalid camera trigger configuration");
    }
    if (config.camera_trigger.enabled &&
        config.camera_trigger.mode == "voltage" &&
        (!config.single_pe.enabled || config.single_pe.unit != "mv")) {
        throw std::runtime_error(
            "voltage trigger requires single_pe.enabled=true and "
            "single_pe.unit=mV");
    }
}

DetectorPipelineResult runDetectorPipeline(
    const DetectorPipelineConfig& config,
    std::size_t n_pixels,
    const std::vector<PrimaryPeHit>& input_hits)
{
    validateDetectorPipelineConfig(config);
    DetectorPipelineResult out;
    out.n_pixels = n_pixels;
    out.primary_hits = input_hits;
    out.pixels.resize(n_pixels);
    out.n_samples = static_cast<std::size_t>(std::ceil(
        (config.sampling.end_ns - config.sampling.start_ns) /
        config.sampling.width_ns));
    out.time_edges_ns.resize(out.n_samples + 1);
    out.time_centers_ns.resize(out.n_samples);
    for (std::size_t i = 0; i <= out.n_samples; ++i) {
        out.time_edges_ns[i] =
            config.sampling.start_ns +
            static_cast<double>(i) * config.sampling.width_ns;
        if (i < out.n_samples) {
            out.time_centers_ns[i] =
                out.time_edges_ns[i] + 0.5 * config.sampling.width_ns;
        }
    }
    if (!input_hits.empty()) {
        out.event_id = input_hits.front().event_id;
        out.telescope_id = input_hits.front().telescope_id;
    }

    std::vector<std::size_t> order(input_hits.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t a, std::size_t b) {
                         const auto& left = input_hits[a];
                         const auto& right = input_hits[b];
                         if (left.time_ns != right.time_ns) {
                             return left.time_ns < right.time_ns;
                         }
                         if (left.pixel_id != right.pixel_id) {
                             return left.pixel_id < right.pixel_id;
                         }
                         return a < b;
                     });

    const std::uint64_t cells_per_pixel =
        static_cast<std::uint64_t>(config.microcell.grid_columns) *
        static_cast<std::uint64_t>(config.microcell.grid_rows);
    std::unordered_set<std::uint64_t> occupied;
    occupied.reserve(input_hits.size());

    for (std::size_t ordered_index : order) {
        const auto& hit = input_hits[ordered_index];
        if (hit.event_id != out.event_id ||
            hit.telescope_id != out.telescope_id) {
            throw std::runtime_error(
                "runDetectorPipeline accepts one event/telescope at a time");
        }
        if (hit.pixel_id < 0 ||
            static_cast<std::size_t>(hit.pixel_id) >= n_pixels ||
            !std::isfinite(hit.time_ns) ||
            !std::isfinite(hit.primary_pe) || hit.primary_pe < 0.0) {
            throw std::runtime_error("invalid primary p.e. hit");
        }
        addPrimary(out.pixels[static_cast<std::size_t>(hit.pixel_id)],
                   hit.origin, hit.primary_pe);
        if (!config.microcell.enabled) {
            if (hit.primary_pe > 0.0) {
                out.fired_hits.push_back({
                    hit.event_id, hit.telescope_id, hit.pixel_id, hit.time_ns,
                    -1, -1, hit.primary_pe, hit.origin});
                addFired(out.pixels[static_cast<std::size_t>(hit.pixel_id)],
                         hit.origin, hit.primary_pe);
            }
            continue;
        }

        const double rounded = std::round(hit.primary_pe);
        if (std::abs(hit.primary_pe - rounded) > 1.0e-9) {
            throw std::runtime_error(
                "explicit microcell tracking requires integral primary_pe; "
                "use stochastic_pe input");
        }
        const int repetitions = static_cast<int>(rounded);
        const auto address = mapMicrocellPosition(
            config.microcell, hit.sensor_x_m, hit.sensor_y_m);
        if (!address.inside_sensor) {
            throw std::runtime_error(
                "primary p.e. position lies outside configured SiPM area");
        }
        if (!address.inside_channel) {
            for (int repetition = 0; repetition < repetitions; ++repetition) {
                if (config.save_microcell_decisions) {
                    out.microcell_decisions.push_back({
                        ordered_index,
                        hit.event_id,
                        hit.telescope_id,
                        hit.pixel_id,
                        hit.time_ns,
                        hit.sensor_x_m,
                        hit.sensor_y_m,
                        address.grid_column,
                        address.grid_row,
                        address.channel_id,
                        address.microcell_id,
                        false,
                        true,
                        false,
                        hit.origin,
                    });
                }
                out.pixels[static_cast<std::size_t>(hit.pixel_id)]
                    .gap_lost_pe += 1.0;
            }
            continue;
        }
        const std::uint64_t cell_in_pixel =
            static_cast<std::uint64_t>(address.channel_id) *
                static_cast<std::uint64_t>(
                    config.microcell.microcells_per_channel) +
            static_cast<std::uint64_t>(address.microcell_id);
        const std::uint64_t occupancy_key =
            static_cast<std::uint64_t>(hit.pixel_id) * cells_per_pixel +
            cell_in_pixel;
        for (int repetition = 0; repetition < repetitions; ++repetition) {
            const bool fired =
                !config.microcell.saturation_enabled ||
                occupied.insert(occupancy_key).second;
            if (config.save_microcell_decisions) {
                out.microcell_decisions.push_back({
                    ordered_index,
                    hit.event_id,
                    hit.telescope_id,
                    hit.pixel_id,
                    hit.time_ns,
                    hit.sensor_x_m,
                    hit.sensor_y_m,
                    address.grid_column,
                    address.grid_row,
                    address.channel_id,
                    address.microcell_id,
                    fired,
                    false,
                    !fired,
                    hit.origin,
                });
            }
            if (fired) {
                out.fired_hits.push_back({
                    hit.event_id, hit.telescope_id, hit.pixel_id, hit.time_ns,
                    address.channel_id, address.microcell_id, 1.0, hit.origin});
                addFired(out.pixels[static_cast<std::size_t>(hit.pixel_id)],
                         hit.origin, 1.0);
            } else {
                out.pixels[static_cast<std::size_t>(hit.pixel_id)]
                    .saturation_lost_pe += 1.0;
            }
        }
    }

    if (config.single_pe.enabled) {
        out.charge_fluctuation_enabled =
            config.single_pe.charge_fluctuation.enabled;
        out.time_jitter_enabled = config.single_pe.time_jitter.enabled;
        out.template_time_reference =
            config.single_pe.template_time_reference;
        std::seed_seq seed{
            static_cast<std::uint32_t>(config.random_seed),
            static_cast<std::uint32_t>(config.random_seed >> 32U),
            static_cast<std::uint32_t>(out.event_id),
            static_cast<std::uint32_t>(out.telescope_id),
            0x5349504dU,
        };
        std::mt19937_64 rng(seed);
        if (config.single_pe.charge_fluctuation.enabled) {
            const auto& samples =
                config.single_pe.charge_fluctuation.empirical_samples;
            std::uniform_int_distribution<std::size_t> choose(
                0, samples.size() - 1);
            for (auto& hit : out.fired_hits) {
                hit.charge_factor = samples[choose(rng)];
            }
        }
        if (config.single_pe.time_jitter.enabled &&
            config.single_pe.time_jitter.sigma_ns > 0.0) {
            std::normal_distribution<double> jitter(
                0.0, config.single_pe.time_jitter.sigma_ns);
            for (auto& hit : out.fired_hits) {
                hit.time_jitter_ns = jitter(rng);
            }
        }
        out.sample_unit =
            config.single_pe.unit == "mv" ? "mV" : "pe_charge_per_sample";
        out.waveform.assign(n_pixels * out.n_samples, 0.0);
        if (config.save_channel_waveforms) {
            out.channel_waveform.assign(
                n_pixels *
                    static_cast<std::size_t>(
                        config.microcell.channels_per_pixel) *
                    out.n_samples,
                0.0);
        }
        const auto pulse = makePulse(config.single_pe,
                                     config.sampling.width_ns);
        out.reference_pulse_time_ns.reserve(pulse.size());
        out.reference_pulse_amplitude.reserve(pulse.size());
        for (const auto& point : pulse) {
            out.reference_pulse_time_ns.push_back(point.time_ns);
            out.reference_pulse_amplitude.push_back(point.amplitude);
        }
        const auto cumulative = cumulativePulse(pulse);
        const double raw_area = trapezoidArea(pulse);
        if (!(raw_area > 0.0)) {
            throw std::runtime_error(
                "single-p.e. pulse must have positive area");
        }
        if (config.single_pe.unit == "mv") {
            out.single_pe_area_mv_ns = raw_area;
        }
        for (const auto& hit : out.fired_hits) {
            const double waveform_time_ns =
                hit.time_ns + hit.time_jitter_ns;
            for (std::size_t bin = 0; bin < out.n_samples; ++bin) {
                const double elapsed_start =
                    out.time_edges_ns[bin] - waveform_time_ns;
                const double elapsed_end =
                    out.time_edges_ns[bin + 1] - waveform_time_ns;
                double sample =
                    cumulativeAt(pulse, cumulative, elapsed_end) -
                    cumulativeAt(pulse, cumulative, elapsed_start);
                if (config.single_pe.unit == "pe_charge") {
                    sample /= raw_area;
                } else {
                    sample /= config.sampling.width_ns;
                }
                sample *= hit.fired_pe * hit.charge_factor;
                if (sample == 0.0) continue;
                const std::size_t pixel =
                    static_cast<std::size_t>(hit.pixel_id);
                out.waveform[pixel * out.n_samples + bin] += sample;
                if (!out.channel_waveform.empty()) {
                    const int channel = hit.channel_id >= 0
                        ? hit.channel_id
                        : 0;
                    const std::size_t flat =
                        (pixel *
                             static_cast<std::size_t>(
                                 config.microcell.channels_per_pixel) +
                         static_cast<std::size_t>(channel)) *
                            out.n_samples +
                        bin;
                    out.channel_waveform[flat] += sample;
                }
            }
        }
    }

    if (config.camera_trigger.mode == "voltage") {
        out.camera_trigger = evaluateVoltageTrigger(
            config, n_pixels, out.n_samples, out.waveform);
    } else {
        out.camera_trigger = evaluatePeCountTrigger(
            config, n_pixels, out.n_samples, out.fired_hits);
    }
    return out;
}

std::vector<PrimaryPeHit> generateUniformNsbPrimaryHits(
    int event_id,
    int telescope_id,
    std::size_t n_pixels,
    double rate_pe_per_ns_per_pixel,
    double start_ns,
    double end_ns,
    const MicrocellConfig& microcell,
    std::uint64_t seed)
{
    if (!(rate_pe_per_ns_per_pixel >= 0.0) ||
        !(end_ns > start_ns) ||
        !(microcell.sensor_size_x_m > 0.0) ||
        !(microcell.sensor_size_y_m > 0.0)) {
        throw std::runtime_error("invalid NSB primary-hit generation settings");
    }
    std::seed_seq sequence{
        static_cast<std::uint32_t>(seed),
        static_cast<std::uint32_t>(seed >> 32U),
        static_cast<std::uint32_t>(event_id),
        static_cast<std::uint32_t>(telescope_id)};
    std::mt19937_64 rng(sequence);
    std::poisson_distribution<int> count_distribution(
        rate_pe_per_ns_per_pixel * (end_ns - start_ns));
    std::uniform_real_distribution<double> time_distribution(start_ns, end_ns);
    std::uniform_real_distribution<double> x_distribution(
        -0.5 * microcell.sensor_size_x_m,
        0.5 * microcell.sensor_size_x_m);
    std::uniform_real_distribution<double> y_distribution(
        -0.5 * microcell.sensor_size_y_m,
        0.5 * microcell.sensor_size_y_m);
    std::uniform_int_distribution<int> channel_distribution(
        0, microcell.channels_per_pixel - 1);
    std::uniform_int_distribution<int> microcell_distribution(
        0, microcell.microcells_per_channel - 1);
    std::uniform_real_distribution<double> centered_distribution(-0.5, 0.5);

    auto active_position = [&]() {
        if (microcell.layout != "s17351_tiled_2x4") {
            return std::pair<double, double>{
                x_distribution(rng), y_distribution(rng)};
        }
        const int channel = channel_distribution(rng);
        const int row_from_top =
            channel < microcell.channel_rows
                ? channel
                : channel - microcell.channel_rows;
        const int channel_column =
            channel < microcell.channel_rows ? 1 : 0;
        const int channel_row_from_bottom =
            microcell.channel_rows - 1 - row_from_top;
        const int cell = microcell_distribution(rng);
        const int cell_column =
            cell % microcell.microcell_columns_per_channel;
        const int cell_row =
            cell / microcell.microcell_columns_per_channel;
        const double pitch_x =
            microcell.channel_size_x_m /
            static_cast<double>(microcell.microcell_columns_per_channel);
        const double pitch_y =
            microcell.channel_size_y_m /
            static_cast<double>(microcell.microcell_rows_per_channel);
        const double x =
            -0.5 * microcell.sensor_size_x_m +
            channel_column *
                (microcell.channel_size_x_m + microcell.channel_gap_x_m) +
            (static_cast<double>(cell_column) + 0.5) * pitch_x +
            centered_distribution(rng) * pitch_x;
        const double y =
            -0.5 * microcell.sensor_size_y_m +
            channel_row_from_bottom *
                (microcell.channel_size_y_m + microcell.channel_gap_y_m) +
            (static_cast<double>(cell_row) + 0.5) * pitch_y +
            centered_distribution(rng) * pitch_y;
        return std::pair<double, double>{x, y};
    };
    std::vector<PrimaryPeHit> hits;
    const double expected =
        static_cast<double>(n_pixels) *
        rate_pe_per_ns_per_pixel * (end_ns - start_ns);
    hits.reserve(static_cast<std::size_t>(std::ceil(expected * 1.1)));
    for (std::size_t pixel = 0; pixel < n_pixels; ++pixel) {
        const int count = count_distribution(rng);
        for (int i = 0; i < count; ++i) {
            const auto position = active_position();
            hits.push_back({
                event_id,
                telescope_id,
                static_cast<int>(pixel),
                time_distribution(rng),
                position.first,
                position.second,
                0.0,
                1.0,
                HitOrigin::Nsb,
            });
        }
    }
    return hits;
}

} // namespace lact::electronics
