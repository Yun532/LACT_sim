#include "app/OpticalSimCommon.hpp"
#include "electronics/DetectorPipeline.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using EventKey = std::pair<int, int>;
using lact::electronics::DetectorPipelineResult;
using lact::electronics::PrimaryPeHit;

std::map<EventKey, std::vector<PrimaryPeHit>> readPrimaryHits(
    const std::string& path,
    int pixel_id_base)
{
    std::map<EventKey, std::vector<PrimaryPeHit>> events;
    if (lact::lowerCopy(path) == "none") return events;
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("failed to open primary p.e. CSV: " + path);
    }
    std::string line;
    while (std::getline(input, line)) {
        line = lact::trim(line);
        if (line.empty() || line[0] == '#') continue;
        const auto cells = lact::splitCsvCells(line);
        if (cells.size() < 7) continue;
        try {
            std::size_t used = 0;
            const int event_id = std::stoi(cells[0], &used);
            if (used != cells[0].size()) continue;
            PrimaryPeHit hit;
            hit.event_id = event_id;
            hit.telescope_id = std::stoi(cells[1]);
            hit.pixel_id = std::stoi(cells[2]) - pixel_id_base;
            hit.time_ns = std::stod(cells[3]);
            hit.sensor_x_m = std::stod(cells[4]);
            hit.sensor_y_m = std::stod(cells[5]);
            hit.primary_pe = std::stod(cells[6]);
            hit.wavelength_nm =
                cells.size() > 7 ? std::stod(cells[7]) : 0.0;
            hit.origin = cells.size() > 8
                ? lact::electronics::parseHitOrigin(cells[8])
                : lact::electronics::HitOrigin::Cherenkov;
            events[{hit.event_id, hit.telescope_id}].push_back(hit);
        } catch (const std::invalid_argument&) {
            // Header row.
        }
    }
    return events;
}

void appendNsb(std::map<EventKey, std::vector<PrimaryPeHit>>& events,
               const std::map<std::string, std::string>& config,
               const lact::electronics::DetectorPipelineConfig& detector,
               std::size_t n_pixels)
{
    const bool enabled = lact::getBool(
        config, "electronics.nsb.enabled",
        lact::getBool(config, "nsb.enabled", false));
    if (!enabled) return;
    const double rate = lact::getDouble(
        config, "electronics.nsb.rate_pe_per_ns_per_pixel",
        lact::getDouble(config, "nsb.rate_pe_per_ns_per_pixel", 0.0));
    const auto automatic_window =
        lact::electronics::waveformContributingPrimaryWindow(detector);
    const double start_ns = lact::getDouble(
        config, "electronics.nsb.start_ns", automatic_window.start_ns);
    const double end_ns = lact::getDouble(
        config, "electronics.nsb.end_ns", automatic_window.end_ns);
    const std::uint64_t seed = lact::getUInt64(
        config, "electronics.nsb.seed",
        lact::getUInt64(config, "nsb.seed", 12345ULL));
    if (!(rate >= 0.0) || !(end_ns > start_ns)) {
        throw std::runtime_error("invalid standalone NSB configuration");
    }
    if (events.empty()) {
        const int event_id =
            lact::getInt(config, "electronics.event_id", 1);
        const int telescope_id =
            lact::getInt(config, "electronics.telescope_id", 0);
        events[{event_id, telescope_id}] = {};
    }
    for (auto& event : events) {
        auto nsb_hits = lact::electronics::generateUniformNsbPrimaryHits(
            event.first.first,
            event.first.second,
            n_pixels,
            rate,
            start_ns,
            end_ns,
            detector.microcell,
            seed);
        for (auto& hit : nsb_hits) {
            hit.count_in_integrated_image =
                hit.time_ns >= detector.sampling.start_ns &&
                hit.time_ns < detector.sampling.end_ns;
        }
        event.second.insert(
            event.second.end(), nsb_hits.begin(), nsb_hits.end());
    }
}

std::ofstream openOutput(const std::filesystem::path& directory,
                         const std::string& name)
{
    std::ofstream output(directory / name);
    if (!output) {
        throw std::runtime_error(
            "failed to create output file: " + (directory / name).string());
    }
    output << std::setprecision(12);
    return output;
}

void writeResult(const std::filesystem::path& directory,
                 const DetectorPipelineResult& result,
                 const lact::electronics::DetectorPipelineConfig& config,
                 int pixel_id_base)
{
    {
        auto output = openOutput(directory, "primary_hits.csv");
        output << "event_id,telescope_id,pixel_id,time_ns,sensor_x_m,"
                  "sensor_y_m,primary_pe,wavelength_nm,origin\n";
        for (const auto& hit : result.primary_hits) {
            output << hit.event_id << ',' << hit.telescope_id << ','
                   << hit.pixel_id + pixel_id_base << ',' << hit.time_ns << ','
                   << hit.sensor_x_m << ',' << hit.sensor_y_m << ','
                   << hit.primary_pe << ',' << hit.wavelength_nm << ','
                   << lact::electronics::hitOriginName(hit.origin) << '\n';
        }
    }
    {
        auto output = openOutput(directory, "microcell_decisions.csv");
        output << "primary_index,event_id,telescope_id,pixel_id,time_ns,"
                  "sensor_x_m,sensor_y_m,grid_column,grid_row,channel_id,"
                  "microcell_id,fired,gap_rejected,saturation_rejected,"
                  "origin\n";
        for (const auto& item : result.microcell_decisions) {
            output << item.primary_index << ',' << item.event_id << ','
                   << item.telescope_id << ','
                   << item.pixel_id + pixel_id_base << ',' << item.time_ns << ','
                   << item.sensor_x_m << ',' << item.sensor_y_m << ','
                   << item.grid_column << ',' << item.grid_row << ','
                   << item.channel_id << ',' << item.microcell_id << ','
                   << (item.fired ? 1 : 0) << ','
                   << (item.gap_rejected ? 1 : 0) << ','
                   << (item.saturation_rejected ? 1 : 0) << ','
                   << lact::electronics::hitOriginName(item.origin) << '\n';
        }
    }
    {
        auto output = openOutput(directory, "fired_hits.csv");
        output << "event_id,telescope_id,pixel_id,time_ns,channel_id,"
                  "microcell_id,fired_pe,origin,charge_factor,"
                  "time_jitter_ns\n";
        for (const auto& hit : result.fired_hits) {
            output << hit.event_id << ',' << hit.telescope_id << ','
                   << hit.pixel_id + pixel_id_base << ',' << hit.time_ns << ','
                   << hit.channel_id << ',' << hit.microcell_id << ','
                   << hit.fired_pe << ','
                   << lact::electronics::hitOriginName(hit.origin) << ','
                   << hit.charge_factor << ',' << hit.time_jitter_ns << '\n';
        }
    }
    {
        auto output = openOutput(directory, "images.csv");
        output << "event_id,telescope_id,pixel_id,primary_cherenkov_pe,"
                  "primary_nsb_pe,primary_dark_pe,primary_total_pe,"
                  "fired_cherenkov_pe,fired_nsb_pe,fired_dark_pe,"
                  "fired_total_pe,gap_lost_pe,saturation_lost_pe\n";
        for (std::size_t pixel = 0; pixel < result.pixels.size(); ++pixel) {
            const auto& item = result.pixels[pixel];
            const double primary =
                item.primary_cherenkov_pe + item.primary_nsb_pe +
                item.primary_dark_pe;
            const double fired =
                item.fired_cherenkov_pe + item.fired_nsb_pe +
                item.fired_dark_pe;
            output << result.event_id << ',' << result.telescope_id << ','
                   << static_cast<int>(pixel) + pixel_id_base << ','
                   << item.primary_cherenkov_pe << ','
                   << item.primary_nsb_pe << ',' << item.primary_dark_pe << ','
                   << primary << ',' << item.fired_cherenkov_pe << ','
                   << item.fired_nsb_pe << ',' << item.fired_dark_pe << ','
                   << fired << ',' << item.gap_lost_pe << ','
                   << item.saturation_lost_pe << '\n';
        }
    }
    {
        auto output = openOutput(directory, "waveform.csv");
        output << "event_id,telescope_id,pixel_id,time_bin,time_center_ns,"
                  "sample_value,sample_unit\n";
        for (std::size_t pixel = 0; pixel < result.n_pixels; ++pixel) {
            for (std::size_t bin = 0; bin < result.n_samples; ++bin) {
                const double value = result.waveform.empty()
                    ? 0.0
                    : result.waveform[pixel * result.n_samples + bin];
                if (value == 0.0) continue;
                output << result.event_id << ',' << result.telescope_id << ','
                       << static_cast<int>(pixel) + pixel_id_base << ',' << bin
                       << ',' << result.time_centers_ns[bin] << ',' << value
                       << ',' << result.sample_unit << '\n';
            }
        }
    }
    if (!result.channel_waveform.empty()) {
        auto output = openOutput(directory, "channel_waveform.csv");
        output << "event_id,telescope_id,pixel_id,channel_id,time_bin,"
                  "time_center_ns,sample_value,sample_unit\n";
        const std::size_t channels = static_cast<std::size_t>(
            config.microcell.channels_per_pixel);
        for (std::size_t pixel = 0; pixel < result.n_pixels; ++pixel) {
            for (std::size_t channel = 0; channel < channels; ++channel) {
                for (std::size_t bin = 0; bin < result.n_samples; ++bin) {
                    const std::size_t flat =
                        (pixel * channels + channel) * result.n_samples + bin;
                    const double value = result.channel_waveform[flat];
                    if (value == 0.0) continue;
                    output << result.event_id << ',' << result.telescope_id
                           << ',' << static_cast<int>(pixel) + pixel_id_base
                           << ',' << channel << ',' << bin << ','
                           << result.time_centers_ns[bin] << ',' << value << ','
                           << result.sample_unit << '\n';
                }
            }
        }
    }
    {
        auto output = openOutput(directory, "trigger.csv");
        output << "event_id,telescope_id,enabled,mode,triggered,"
                  "max_pixels_above_threshold,first_trigger_bin,"
                  "trigger_time_ns\n";
        output << result.event_id << ',' << result.telescope_id << ','
               << (config.camera_trigger.enabled ? 1 : 0) << ','
               << config.camera_trigger.mode << ','
               << (result.camera_trigger.triggered ? 1 : 0) << ','
               << result.camera_trigger.max_pixels_above_threshold << ','
               << result.camera_trigger.first_trigger_bin << ','
               << result.camera_trigger.trigger_time_ns << '\n';
        auto multiplicity =
            openOutput(directory, "trigger_multiplicity.csv");
        multiplicity << "event_id,telescope_id,time_bin,time_center_ns,"
                        "pixels_above_threshold\n";
        for (std::size_t bin = 0;
             bin < result.camera_trigger.pixels_above_threshold.size();
             ++bin) {
            multiplicity << result.event_id << ',' << result.telescope_id
                         << ',' << bin << ',' << result.time_centers_ns[bin]
                         << ','
                         << result.camera_trigger.pixels_above_threshold[bin]
                         << '\n';
        }
    }
    {
        auto output = openOutput(directory, "metadata.json");
        output << "{\n"
               << "  \"event_id\": " << result.event_id << ",\n"
               << "  \"telescope_id\": " << result.telescope_id << ",\n"
               << "  \"n_pixels\": " << result.n_pixels << ",\n"
               << "  \"n_primary_hits\": " << result.primary_hits.size()
               << ",\n"
               << "  \"n_fired_hits\": " << result.fired_hits.size() << ",\n"
               << "  \"microcell_enabled\": "
               << (config.microcell.enabled ? "true" : "false") << ",\n"
               << "  \"microcell_saturation_enabled\": "
               << (config.microcell.saturation_enabled ? "true" : "false")
               << ",\n"
               << "  \"microcell_layout\": \""
               << config.microcell.layout << "\",\n"
               << "  \"inter_channel_gap_geometry\": "
               << "\"part of tiled layout\",\n"
               << "  \"pde_includes_inter_channel_gaps\": "
               << (config.microcell.pde_includes_inter_channel_gaps
                       ? "true"
                       : "false")
               << ",\n"
               << "  \"inter_channel_active_fraction\": "
               << ::lact::electronics::interChannelActiveFraction(
                      config.microcell)
               << ",\n"
               << "  \"single_pe_enabled\": "
               << (config.single_pe.enabled ? "true" : "false") << ",\n"
               << "  \"single_pe_area_mv_ns\": "
               << result.single_pe_area_mv_ns << ",\n"
               << "  \"single_pe_template_time_reference\": \""
               << result.template_time_reference << "\",\n"
               << "  \"single_pe_charge_fluctuation_enabled\": "
               << (result.charge_fluctuation_enabled ? "true" : "false")
               << ",\n"
               << "  \"single_pe_time_jitter_enabled\": "
               << (result.time_jitter_enabled ? "true" : "false")
               << ",\n"
               << "  \"sample_unit\": \"" << result.sample_unit << "\",\n"
               << "  \"sample_width_ns\": " << config.sampling.width_ns
               << ",\n"
               << "  \"camera_trigger_enabled\": "
               << (config.camera_trigger.enabled ? "true" : "false") << ",\n"
               << "  \"camera_trigger_mode\": \""
               << config.camera_trigger.mode << "\",\n"
               << "  \"camera_triggered\": "
               << (result.camera_trigger.triggered ? "true" : "false")
               << "\n}\n";
    }
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const auto command = lact::parseConfigCommandLine(argc, argv);
        if (command.help || command.positional.size() != 3) {
            std::cerr
                << "usage: run_camera_electronics CONFIG "
                   "PRIMARY_HITS.csv|none OUTPUT_DIR [-C key=value ...]\n";
            return command.help ? 0 : 2;
        }
        const std::string config_path = command.positional[0];
        const std::string primary_path = command.positional[1];
        const std::string output_path = command.positional[2];
        lact::ComponentConfigPaths paths;
        auto main_config = lact::readKeyValueConfig(config_path);
        lact::applyConfigOverrides(main_config, command.overrides);
        auto config = lact::expandConfig(main_config, config_path, paths);
        lact::applyConfigOverrides(config, command.overrides);
        const auto detector = lact::buildDetectorPipelineConfig(config);
        const std::size_t n_pixels = static_cast<std::size_t>(
            lact::getInt(config, "electronics.n_pixels", 1664));
        const int pixel_id_base =
            lact::getInt(config, "electronics.input.pixel_id_base", 0);
        auto events = readPrimaryHits(primary_path, pixel_id_base);
        appendNsb(events, config, detector, n_pixels);
        if (events.empty()) {
            throw std::runtime_error(
                "no primary p.e. input and standalone NSB is disabled");
        }
        const std::filesystem::path root(output_path);
        std::filesystem::create_directories(root);
        for (const auto& event : events) {
            const auto result = lact::electronics::runDetectorPipeline(
                detector, n_pixels, event.second);
            const std::filesystem::path directory =
                events.size() == 1
                    ? root
                    : root / ("event_" + std::to_string(event.first.first) +
                              "_tel_" + std::to_string(event.first.second));
            std::filesystem::create_directories(directory);
            writeResult(directory, result, detector, pixel_id_base);
        }
        std::cout << "electronics diagnostics written: " << root.string()
                  << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
