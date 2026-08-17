#include "app/OpticalSimCommon.hpp"

#include <cmath>
#include <filesystem>
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

bool expectInvalidEfficiency(const std::map<std::string, std::string>& cfg,
                             const std::string& label)
{
    try {
        (void)buildEfficiencyConfig(cfg);
    } catch (...) {
        return true;
    }
    std::cerr << "expected invalid efficiency config: " << label << "\n";
    return false;
}

bool expectInvalidDetector(const std::map<std::string, std::string>& cfg,
                           const std::string& label)
{
    try {
        (void)buildDetectorPipelineConfig(cfg);
    } catch (...) {
        return true;
    }
    std::cerr << "expected invalid detector config: " << label << "\n";
    return false;
}

bool expectInvalidCamera(const std::map<std::string, std::string>& cfg,
                         const std::string& label)
{
    try {
        (void)buildCameraConfig(cfg);
    } catch (...) {
        return true;
    }
    std::cerr << "expected invalid camera config: " << label << "\n";
    return false;
}

} // namespace

int main()
{
    bool ok = true;

    const auto explicit_whiteboard = buildCameraConfig({
        {"camera.mode", "whiteboard"},
    });
    ok &= check(explicit_whiteboard.whiteboard &&
                    !explicit_whiteboard.enabled &&
                    !explicit_whiteboard.implicit_whiteboard_legacy,
                "explicit whiteboard mode should be distinct from a camera");
    const auto legacy_whiteboard = buildCameraConfig({});
    ok &= check(!legacy_whiteboard.whiteboard &&
                    !legacy_whiteboard.enabled &&
                    legacy_whiteboard.implicit_whiteboard_legacy,
                "missing camera mode should retain legacy whiteboard behavior");
    ok &= expectInvalidCamera(
        {{"camera.mode", "whiteboard"}, {"camera.enabled", "true"}},
        "whiteboard with camera enabled");
    ok &= expectInvalidCamera(
        {{"camera.mode", "whiteboard"}, {"output.mode", "pixel"}},
        "whiteboard with pixel output");
    ok &= expectInvalidCamera(
        {{"camera.mode", "whiteboard"}, {"electronics.enabled", "true"}},
        "whiteboard with electronics");
    ok &= expectInvalidCamera(
        {{"camera.mode", "whiteboard"}, {"trigger.enabled", "true"}},
        "whiteboard with trigger");

    const auto default_detector = buildDetectorPipelineConfig({});
    ok &= check(!default_detector.enabled,
                "electronics should default to disabled");
    ok &= check(
        default_detector.microcell.pde_includes_inter_channel_gaps,
        "inter-channel-gap PDE convention should default to included");

    const auto enabled_detector = buildDetectorPipelineConfig({
        {"electronics.enabled", "true"},
    });
    ok &= check(enabled_detector.enabled,
                "electronics.enabled should enable detector simulation");

    ok &= expectInvalidDetector(
        {{"electronics.pipeline.enabled", "true"}},
        "removed electronics.pipeline.enabled key");
    ok &= expectInvalidDetector(
        {{"electronics.enabled", "false"},
         {"electronics.pipeline.enabled", "true"}},
        "removed key must fail even when electronics.enabled is present");

    auto default_nsb = buildNsbConfig({});
    ok &= check(!default_nsb.enabled, "NSB should default to disabled");
    ok &= check(default_nsb.model == "constant_rate", "NSB default model");
    ok &= check(std::abs(default_nsb.window_ns - 32.0) < 1e-12, "NSB default window");
    auto enabled_default_nsb = default_nsb;
    enabled_default_nsb.enabled = true;
    ok &= check(nsbTimeInImageWindow(enabled_default_nsb, 0.0),
                "NSB image gate includes its lower edge");
    ok &= check(nsbTimeInImageWindow(enabled_default_nsb, 31.999),
                "NSB image gate includes times below its upper edge");
    ok &= check(!nsbTimeInImageWindow(enabled_default_nsb, -0.001) &&
                    !nsbTimeInImageWindow(enabled_default_nsb, 32.0),
                "NSB image gate excludes times outside [0, window_ns)");

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
        {"nsb.collector_mean_transmission", "0.923436437"},
        {"nsb.pixel_solid_angle", "auto"},
    };
    auto spectral_nsb = buildNsbConfig(spectral_nsb_cfg);
    ok &= check(spectral_nsb.enabled, "spectral NSB enabled parsed");
    ok &= check(spectral_nsb.model == "spectral_flux", "spectral NSB model parsed");
    ok &= check(spectral_nsb.spectrum_csv == "configs/nsb/nsb_spectrum.csv",
                "spectral NSB spectrum path parsed");
    ok &= check(std::abs(spectral_nsb.effective_area_m2 - 24.576860) < 1e-12,
                "spectral NSB effective area parsed");
    ok &= check(std::abs(spectral_nsb.collector_mean_transmission -
                         0.923436437) < 1e-12,
                "spectral NSB collector transmission parsed");

    const auto s17351 = buildDetectorPipelineConfig({
        {"electronics.enabled", "true"},
        {"electronics.microcell.enabled", "true"},
        {"electronics.microcell.saturation_enabled", "false"},
        {"electronics.microcell.layout", "s17351_tiled_2x4"},
        {"electronics.microcell.sensor_size_x_m", "0.0134"},
        {"electronics.microcell.sensor_size_y_m", "0.0134"},
        {"electronics.microcell.channel_columns", "2"},
        {"electronics.microcell.channel_rows", "4"},
        {"electronics.microcell.channel_size_x_m", "0.0066"},
        {"electronics.microcell.channel_size_y_m", "0.0032"},
        {"electronics.microcell.channel_gap_x_m", "0.0002"},
        {"electronics.microcell.channel_gap_y_m", "0.0002"},
        {"electronics.microcell.microcell_columns_per_channel", "264"},
        {"electronics.microcell.microcell_rows_per_channel", "128"},
        {"electronics.microcell.pde_includes_inter_channel_gaps", "true"},
    });
    ok &= check(s17351.microcell.enabled,
                "S17351 microcell geometry enabled parsed");
    ok &= check(!s17351.microcell.saturation_enabled,
                "independent saturation switch parsed");
    ok &= check(s17351.microcell.layout == "s17351_tiled_2x4",
                "physical S17351 layout parsed");
    ok &= check(s17351.microcell.pde_includes_inter_channel_gaps,
                "inter-channel-gap PDE convention parsed");
    CameraConfig camera_without_collector;
    camera_without_collector.enabled = true;
    try {
        validateCameraDetectorCompatibility(camera_without_collector, s17351);
        std::cerr << "explicit S17351 accepted a camera without collector\n";
        ok = false;
    } catch (...) {
    }
    CameraConfig camera_with_collector = camera_without_collector;
    camera_with_collector.collector = "bezier";
    try {
        validateCameraDetectorCompatibility(camera_with_collector, s17351);
    } catch (const std::exception& error) {
        std::cerr << "explicit S17351 rejected configured collector: "
                  << error.what() << '\n';
        ok = false;
    }

    const auto spectrum_path = std::filesystem::temp_directory_path() /
        "lact_nsb_collector_factor_test.csv";
    {
        std::ofstream spectrum(spectrum_path);
        spectrum << "wavelength_nm,flux\n300,1\n400,1\n";
    }
    CameraGeometry one_pixel;
    CameraPixel pixel;
    pixel.id = 0;
    pixel.center = {0.0, 0.0, 0.0};
    pixel.size = 0.0244;
    pixel.shape = PixelShape::Square;
    one_pixel.addPixel(pixel);
    TelescopeConfig nsb_telescope;
    auto make_spectral = [&]() {
        return buildNsbConfig({
            {"nsb.enabled", "true"},
            {"nsb.model", "spectral_flux"},
            {"nsb.spectrum_csv", spectrum_path.string()},
            {"nsb.effective_area_m2", "1"},
            {"nsb.collector_mean_transmission", "0.5"},
            {"nsb.pixel_solid_angle", "manual"},
            {"nsb.pixel_solid_angle_sr", "1"},
        });
    };
    auto package_pde_nsb = make_spectral();
    resolveNsbSpectralRate(package_pde_nsb, {}, one_pixel, nsb_telescope,
                           &s17351);
    ok &= check(std::abs(package_pde_nsb.rate_pe_per_ns_per_pixel - 5.0e-8) <
                    1.0e-18,
                "collector transmission was not applied to spectral NSB");
    ok &= check(std::abs(package_pde_nsb.microcell_geometric_acceptance - 1.0) <
                    1.0e-12,
                "full-package PDE must not lose the channel gap twice");
    auto intrinsic_pde_detector = s17351;
    intrinsic_pde_detector.microcell.pde_includes_inter_channel_gaps = false;
    auto intrinsic_pde_nsb = make_spectral();
    resolveNsbSpectralRate(intrinsic_pde_nsb, {}, one_pixel, nsb_telescope,
                           &intrinsic_pde_detector);
    const double active_fraction =
        electronics::interChannelActiveFraction(s17351.microcell);
    ok &= check(std::abs(intrinsic_pde_nsb.rate_pe_per_ns_per_pixel /
                             package_pde_nsb.rate_pe_per_ns_per_pixel -
                         active_fraction) < 1.0e-12,
                "intrinsic PDE spectral NSB did not apply channel active area");
    std::filesystem::remove(spectrum_path);

    auto trigger_default = buildTriggerConfig({});
    ok &= check(!trigger_default.enabled, "trigger should default to disabled");
    ok &= check(std::abs(trigger_default.pixel_threshold_pe - 5.0) < 1e-12,
                "trigger default threshold");
    ok &= check(std::abs(trigger_default.camera_coincidence_window_ns - 20.0) < 1e-12,
                "camera trigger should default to a 20 ns window");
    ok &= check(std::abs(trigger_default.array_coincidence_window_ns) < 1e-12,
                "raw array timing should default to disabled");
    ok &= check(trigger_default.array_time_correction == "none",
                "array timing correction should default to none");
    ok &= check(std::abs(trigger_default.array_wavefront_speed_m_per_ns) < 1e-12,
                "array wavefront speed should default to observation-altitude auto mode");
    const auto split_windows = buildTriggerConfig({
        {"trigger.coincidence_window_ns", "50"},
        {"trigger.camera_coincidence_window_ns", "8"},
        {"trigger.array_coincidence_window_ns", "24"},
        {"trigger.array_time_correction", "plane_wave"},
        {"trigger.array_wavefront_speed_m_per_ns", "0.2997"},
    });
    ok &= check(std::abs(split_windows.camera_coincidence_window_ns - 8.0) < 1e-12,
                "camera coincidence window parsed");
    ok &= check(std::abs(split_windows.array_coincidence_window_ns - 24.0) < 1e-12,
                "array coincidence window parsed");
    ok &= check(split_windows.array_time_correction == "plane_wave",
                "plane-wave array timing correction parsed");
    ok &= check(std::abs(split_windows.array_wavefront_speed_m_per_ns - 0.2997) < 1e-12,
                "array wavefront speed parsed");

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
    ok &= check(std::abs(trigger.camera_coincidence_window_ns - 30.0) < 1e-12,
                "legacy coincidence window should configure the camera window");
    ok &= check(std::abs(trigger.array_coincidence_window_ns - 30.0) < 1e-12,
                "legacy coincidence window should configure the array window");

    ok &= expectInvalid({{"nsb.rate_pe_per_ns_per_pixel", "-1"}}, "negative NSB rate");
    ok &= expectInvalid({{"nsb.collector_mean_transmission", "0"}},
                        "zero collector transmission");
    ok &= expectInvalid({{"nsb.collector_mean_transmission", "1.01"}},
                        "collector transmission above one");
    ok &= expectInvalid({{"nsb.model", "spectral_flux"},
                         {"nsb.spectrum_csv", "spectrum.csv"},
                         {"nsb.effective_area_m2", "1"}},
                        "spectral NSB missing collector transmission");
    ok &= expectInvalid({{"nsb.enabled", "true"}, {"nsb.model", "spectral_flux"}},
                        "spectral NSB missing spectrum");
    ok &= expectInvalid({{"trigger.camera_multiplicity", "0"}}, "zero camera multiplicity");
    ok &= expectInvalid({{"trigger.array_time_correction", "unknown"}},
                        "unknown array timing correction");
    ok &= expectInvalidEfficiency(
        {{"efficiency.mirror_reflectivity", "1.1"}},
        "reflectivity above one");
    ok &= expectInvalidEfficiency(
        {{"atmosphere.model", "modtran_tau"},
         {"efficiency.atmosphere_transmission", "0.9"}},
        "duplicate atmosphere response");
    ok &= expectInvalidEfficiency(
        {{"camera.collector", "bezier"},
         {"efficiency.funnel_acceptance", "true"}},
        "duplicate collector angular response");
    const auto intentional_atmosphere = buildEfficiencyConfig(
        {{"atmosphere.model", "modtran_tau"},
         {"efficiency.atmosphere_transmission", "0.9"},
         {"atmosphere.allow_additional_spectral_factor", "true"}});
    ok &= check(intentional_atmosphere.atmosphere_transmission.enabled,
                "intentional additional atmosphere factor should be explicit");

    std::cout << "detector response config checks passed\n";
    return ok ? 0 : 1;
}
