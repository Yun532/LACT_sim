#include "app/OpticalSimCommon.hpp"

#include <exception>
#include <iostream>

using namespace lact;

int main(int argc, char** argv)
{
    try {
        const auto command = parseConfigCommandLine(argc, argv);
        if (command.help || command.positional.size() != 1) {
            std::cerr
                << "usage: compute_nsb_rate CONFIG.cfg [-C key=value ...]\n";
            return command.help ? 0 : 2;
        }
        const std::string config_path = command.positional[0];
        ComponentConfigPaths component_paths;
        auto main_cfg = readKeyValueConfig(config_path);
        applyConfigOverrides(main_cfg, command.overrides);
        auto cfg = expandConfig(main_cfg, config_path, component_paths);
        applyConfigOverrides(cfg, command.overrides);

        TelescopeConfig telescope_cfg = buildTelescopeConfig(cfg);
        CameraConfig camera_cfg = buildCameraConfig(cfg);
        CameraGeometry camera = buildCameraGeometry(camera_cfg);
        OpticalEfficiencyConfig efficiency_cfg = buildEfficiencyConfig(cfg);
        const auto detector_cfg = buildDetectorPipelineConfig(cfg);
        NsbConfig nsb_cfg = buildNsbConfig(cfg);
        resolveNsbSpectralRate(nsb_cfg, efficiency_cfg, camera, telescope_cfg,
                               &detector_cfg);

        printSection("NSB spectral rate");
        printField("config", config_path);
        for (const auto& [key, value] : command.overrides) {
            printField("override", key + "=" + value);
        }
        if (!component_paths.nsb.empty()) printField("nsb_config", component_paths.nsb);
        if (!component_paths.camera.empty()) printField("camera_config", component_paths.camera);
        if (!component_paths.sipm.empty()) printField("sipm_config", component_paths.sipm);
        if (!component_paths.efficiency.empty()) {
            printField("efficiency_config", component_paths.efficiency);
        }
        printField("enabled", nsb_cfg.enabled ? "true" : "false");
        printField("model", nsb_cfg.model);
        printField("spectrum_csv", nsb_cfg.spectrum_csv);
        printField("spectrum_unit", nsb_cfg.spectrum_unit);
        printField("effective_area_m2", doubleToString(nsb_cfg.effective_area_m2, 9));
        printField("collector_mean_transmission",
                   doubleToString(nsb_cfg.collector_mean_transmission, 9));
        printField("microcell_geometric_acceptance",
                   doubleToString(nsb_cfg.microcell_geometric_acceptance, 9));
        printField("pixel_solid_angle_sr", doubleToString(nsb_cfg.pixel_solid_angle_sr, 12));
        printField("spectral_integral_pe_s_sr_m2",
                   doubleToString(nsb_cfg.spectral_integral_pe_s_sr_m2, 3));
        printField("rate_pe_per_ns_per_pixel",
                   doubleToString(nsb_cfg.rate_pe_per_ns_per_pixel, 9));
        printField("mean_pe_per_pixel_25ns",
                   doubleToString(nsb_cfg.rate_pe_per_ns_per_pixel * 25.0, 9));
        printField("mean_pe_per_pixel_30ns",
                   doubleToString(nsb_cfg.rate_pe_per_ns_per_pixel * 30.0, 9));
    } catch (const std::exception& e) {
        std::cerr << "compute_nsb_rate error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
