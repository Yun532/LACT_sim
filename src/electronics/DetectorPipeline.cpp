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
    const double onset = pulse.front().time_ns;
    for (auto& point : pulse) point.time_ns -= onset;
    return pulse;
}

std::vector<PulsePoint> makePulse(const SinglePeConfig& config,
                                  double sample_width_ns)
{
    std::vector<PulsePoint> pulse;
    if (config.model == "measured_csv") {
        pulse = loadPulseCsv(config.csv_path);
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
    pulse.erase(
        std::remove_if(pulse.begin(), pulse.end(),
                       [&](const PulsePoint& point) {
                           return point.time_ns > config.support_ns;
                       }),
        pulse.end());
    if (pulse.size() < 2) {
        throw std::runtime_error(
            "single-p.e. pulse has fewer than two points inside support");
    }
    if (pulse.back().time_ns < config.support_ns) {
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
        for (std::size_t start = 0; start < n_samples; ++start) {
            double pe = 0.0;
            const std::size_t end = std::min(n_samples, start + window_bins);
            for (std::size_t bin = start; bin < end; ++bin) {
                pe += counts[pixel * n_samples + bin];
            }
            if (pe >= config.camera_trigger.pixel_threshold_pe) {
                ++out.pixels_above_threshold[start];
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
        for (std::size_t start = 0; start < n_samples; ++start) {
            const std::size_t end = std::min(n_samples, start + window_bins);
            bool crossed = false;
            for (std::size_t bin = start; bin < end; ++bin) {
                if (waveform[pixel * n_samples + bin] >=
                    config.camera_trigger.pixel_threshold_mv) {
                    crossed = true;
                    break;
                }
            }
            if (crossed) ++out.pixels_above_threshold[start];
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
    if (config.microcell.layout != "uniform_interleaved") {
        throw std::runtime_error(
            "only microcell.layout=uniform_interleaved is currently supported");
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
        const double half_x = 0.5 * config.microcell.sensor_size_x_m;
        const double half_y = 0.5 * config.microcell.sensor_size_y_m;
        if (hit.sensor_x_m < -half_x || hit.sensor_x_m > half_x ||
            hit.sensor_y_m < -half_y || hit.sensor_y_m > half_y) {
            throw std::runtime_error(
                "primary p.e. position lies outside configured SiPM area");
        }
        const double x_fraction =
            std::clamp((hit.sensor_x_m + half_x) /
                           config.microcell.sensor_size_x_m,
                       0.0, std::nextafter(1.0, 0.0));
        const double y_fraction =
            std::clamp((hit.sensor_y_m + half_y) /
                           config.microcell.sensor_size_y_m,
                       0.0, std::nextafter(1.0, 0.0));
        const int column = static_cast<int>(
            x_fraction * config.microcell.grid_columns);
        const int row = static_cast<int>(
            y_fraction * config.microcell.grid_rows);
        const std::uint64_t global_cell =
            static_cast<std::uint64_t>(row) *
                static_cast<std::uint64_t>(config.microcell.grid_columns) +
            static_cast<std::uint64_t>(column);
        const int channel = static_cast<int>(
            global_cell %
            static_cast<std::uint64_t>(config.microcell.channels_per_pixel));
        const int microcell = static_cast<int>(
            global_cell /
            static_cast<std::uint64_t>(config.microcell.channels_per_pixel));
        const std::uint64_t occupancy_key =
            static_cast<std::uint64_t>(hit.pixel_id) * cells_per_pixel +
            global_cell;
        for (int repetition = 0; repetition < repetitions; ++repetition) {
            const bool fired = occupied.insert(occupancy_key).second;
            out.microcell_decisions.push_back({
                ordered_index,
                hit.event_id,
                hit.telescope_id,
                hit.pixel_id,
                hit.time_ns,
                hit.sensor_x_m,
                hit.sensor_y_m,
                column,
                row,
                channel,
                microcell,
                fired,
                hit.origin,
            });
            if (fired) {
                out.fired_hits.push_back({
                    hit.event_id, hit.telescope_id, hit.pixel_id, hit.time_ns,
                    channel, microcell, 1.0, hit.origin});
                addFired(out.pixels[static_cast<std::size_t>(hit.pixel_id)],
                         hit.origin, 1.0);
            } else {
                out.pixels[static_cast<std::size_t>(hit.pixel_id)]
                    .saturation_lost_pe += 1.0;
            }
        }
    }

    if (config.single_pe.enabled) {
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
        const auto cumulative = cumulativePulse(pulse);
        const double raw_area = trapezoidArea(pulse);
        if (!(raw_area > 0.0)) {
            throw std::runtime_error(
                "single-p.e. pulse must have positive area");
        }
        for (const auto& hit : out.fired_hits) {
            for (std::size_t bin = 0; bin < out.n_samples; ++bin) {
                const double elapsed_start =
                    out.time_edges_ns[bin] - hit.time_ns;
                const double elapsed_end =
                    out.time_edges_ns[bin + 1] - hit.time_ns;
                double sample =
                    cumulativeAt(pulse, cumulative, elapsed_end) -
                    cumulativeAt(pulse, cumulative, elapsed_start);
                if (config.single_pe.unit == "pe_charge") {
                    sample /= raw_area;
                } else {
                    sample /= config.sampling.width_ns;
                }
                sample *= hit.fired_pe;
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
    double sensor_size_x_m,
    double sensor_size_y_m,
    std::uint64_t seed)
{
    if (!(rate_pe_per_ns_per_pixel >= 0.0) ||
        !(end_ns > start_ns) ||
        !(sensor_size_x_m > 0.0) ||
        !(sensor_size_y_m > 0.0)) {
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
        -0.5 * sensor_size_x_m, 0.5 * sensor_size_x_m);
    std::uniform_real_distribution<double> y_distribution(
        -0.5 * sensor_size_y_m, 0.5 * sensor_size_y_m);
    std::vector<PrimaryPeHit> hits;
    const double expected =
        static_cast<double>(n_pixels) *
        rate_pe_per_ns_per_pixel * (end_ns - start_ns);
    hits.reserve(static_cast<std::size_t>(std::ceil(expected * 1.1)));
    for (std::size_t pixel = 0; pixel < n_pixels; ++pixel) {
        const int count = count_distribution(rng);
        for (int i = 0; i < count; ++i) {
            hits.push_back({
                event_id,
                telescope_id,
                static_cast<int>(pixel),
                time_distribution(rng),
                x_distribution(rng),
                y_distribution(rng),
                0.0,
                1.0,
                HitOrigin::Nsb,
            });
        }
    }
    return hits;
}

} // namespace lact::electronics
