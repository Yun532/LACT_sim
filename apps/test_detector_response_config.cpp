#include "app/OpticalSimCommon.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

using namespace lact;

namespace {

bool check(bool cond, const std::string& msg)
{
    if (!cond) {
        std::cerr << msg << "\n";
        return false;
    }
    return cond;
}

bool expectInvalid(const std::map<std::string, std::string>& cfg,
                   const std::string& label)
{
    try {
        (void)buildNsbConfig(cfg);
        (void)buildTriggerConfig(cfg);
    } catch (...) {
        return true;
    }
    std::cerr << "expected invalid config: " << label << "\n";
    return false;
}

} // namespace

int main()
{
    bool ok = true;

    auto default_nsb = buildNsbConfig({});
    ok &= check(!default_nsb.enabled, "NSB should default to disabled");
    ok &= check(default_nsb.model == "constant_rate", "NSB default model");
    ok &= check(default_nsb.spatial_model == "uniform", "NSB default spatial model");
    ok &= check(std::abs(default_nsb.window_ns - 16.0) < 1e-12, "NSB default window");

    std::map<std::string, std::string> nsb_cfg{
        {"nsb.enabled", "true"},
        {"nsb.model", "constant_rate"},
        {"nsb.rate_pe_per_ns_per_pixel", "0.05"},
        {"nsb.window_ns", "20"},
        {"nsb.seed", "7"},
    };
    auto nsb = buildNsbConfig(nsb_cfg);
    ok &= check(nsb.enabled, "NSB enabled parsed");
    ok &= check(std::abs(nsb.rate_pe_per_ns_per_pixel - 0.05) < 1e-12,
                "NSB rate parsed");
    ok &= check(std::abs(nsb.window_ns - 20.0) < 1e-12, "NSB window parsed");
    ok &= check(nsb.seed == 7, "NSB seed parsed");

    const std::string pixel_scale_path = "test_nsb_pixel_scale.csv";
    {
        std::ofstream out(pixel_scale_path);
        out << "telescope_id,pixel_id,relative_scale\n"
            << "-1,0,1.0\n"
            << "3,0,0.75\n";
    }
    auto scaled_nsb_cfg = nsb_cfg;
    scaled_nsb_cfg["nsb.spatial_model"] = "pixel_scale";
    scaled_nsb_cfg["nsb.pixel_scale_csv"] = pixel_scale_path;
    const auto scaled_nsb = buildNsbConfig(scaled_nsb_cfg);
    ok &= check(scaled_nsb.spatial_model == "pixel_scale",
                "NSB spatial model parsed");
    ok &= check(scaled_nsb.pixel_relative_scale.size() == 2,
                "NSB pixel scale rows parsed");
    ok &= check(std::abs(scaled_nsb.pixel_relative_scale.at({3, 0}) - 0.75) < 1e-12,
                "NSB telescope pixel override parsed");
    std::remove(pixel_scale_path.c_str());

    {
        std::ofstream out(pixel_scale_path);
        out << "telescope_id,pixel_id,relative_scale\n"
            << "-1,0,1.0\n"
            << "-1,0,0.5\n";
    }
    try {
        (void)buildNsbConfig(scaled_nsb_cfg);
        std::cerr << "duplicate NSB pixel scale row should fail\n";
        ok = false;
    } catch (...) {
    }
    std::remove(pixel_scale_path.c_str());

    const auto event_ids =
        parseIntList("603, 201,603,402", "source.filter_event_ids");
    ok &= check(event_ids == std::vector<int>({201, 402, 603}),
                "event ID list is parsed, sorted, and deduplicated");
    try {
        (void)parseIntList("201,bad", "source.filter_event_ids");
        std::cerr << "expected invalid event ID list\n";
        ok = false;
    } catch (...) {
    }
    const auto selected_events = buildSourceRuntimeConfig({
        {"source.filter_event_ids", "603,201,603,402"},
    });
    ok &= check(selected_events.selected_event_ids ==
                    std::vector<int>({201, 402, 603}),
                "multiple source event filters are retained");
    try {
        (void)buildSourceRuntimeConfig({
            {"source.filter_event_id", "201"},
            {"source.filter_event_ids", "201,402"},
        });
        std::cerr << "expected mutually exclusive event filters\n";
        ok = false;
    } catch (...) {
    }

    std::map<std::string, std::string> spectral_nsb_cfg{
        {"nsb.enabled", "true"},
        {"nsb.model", "spectral_flux"},
        {"nsb.spectrum_csv", "configs/nsb/nsb_spectrum.csv"},
        {"nsb.spectrum_unit", "ph_s_nm_sr_m2"},
        {"nsb.effective_area_m2", "24.576860"},
        {"nsb.pixel_solid_angle", "auto"},
    };
    auto spectral_nsb = buildNsbConfig(spectral_nsb_cfg);
    ok &= check(spectral_nsb.enabled, "spectral NSB enabled parsed");
    ok &= check(spectral_nsb.model == "spectral_flux", "spectral NSB model parsed");
    ok &= check(spectral_nsb.spectrum_csv == "configs/nsb/nsb_spectrum.csv",
                "spectral NSB spectrum path parsed");
    ok &= check(std::abs(spectral_nsb.effective_area_m2 - 24.576860) < 1e-12,
                "spectral NSB effective area parsed");

    auto trigger_default = buildTriggerConfig({});
    ok &= check(!trigger_default.enabled, "trigger should default to disabled");
    ok &= check(std::abs(trigger_default.pixel_threshold_pe - 5.0) < 1e-12,
                "trigger default threshold");
    const auto split_windows = buildTriggerConfig({
        {"trigger.coincidence_window_ns", "50"},
        {"trigger.camera_coincidence_window_ns", "8"},
        {"trigger.array_coincidence_window_ns", "24"},
    });
    ok &= check(std::abs(split_windows.camera_coincidence_window_ns - 8.0) < 1e-12,
                "camera coincidence window parsed");
    ok &= check(std::abs(split_windows.array_coincidence_window_ns - 24.0) < 1e-12,
                "array coincidence window parsed");

    std::map<std::string, std::string> trigger_cfg{
        {"trigger.enabled", "true"},
        {"trigger.pixel_threshold_pe", "3.5"},
        {"trigger.camera_multiplicity", "4"},
        {"trigger.array_multiplicity", "2"},
        {"trigger.coincidence_window_ns", "30"},
    };
    auto trigger = buildTriggerConfig(trigger_cfg);
    ok &= check(trigger.enabled, "trigger enabled parsed");
    ok &= check(std::abs(trigger.pixel_threshold_pe - 3.5) < 1e-12,
                "trigger threshold parsed");
    ok &= check(trigger.camera_multiplicity == 4, "camera multiplicity parsed");
    ok &= check(trigger.array_multiplicity == 2, "array multiplicity parsed");
    ok &= check(std::abs(trigger.coincidence_window_ns - 30.0) < 1e-12,
                "coincidence window parsed");

    ok &= expectInvalid({{"nsb.rate_pe_per_ns_per_pixel", "-1"}}, "negative NSB rate");
    ok &= expectInvalid({{"nsb.spatial_model", "invalid"}}, "invalid NSB spatial model");
    ok &= expectInvalid({{"nsb.spatial_model", "pixel_scale"}},
                        "pixel scale model missing CSV");
    ok &= expectInvalid({{"nsb.enabled", "true"}, {"nsb.model", "spectral_flux"}},
                        "spectral NSB missing spectrum");
    ok &= expectInvalid({{"trigger.camera_multiplicity", "0"}}, "zero camera multiplicity");

    std::cout << "detector response config checks passed\n";
    return ok ? 0 : 1;
}
