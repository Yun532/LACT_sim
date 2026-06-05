#include "app/OpticalSimCommon.hpp"

#include <exception>
#include <iostream>

using namespace lact;

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "usage: compute_nsb_rate CONFIG.cfg\n";
        return 2;
    }

    try {
        ComponentConfigPaths component_paths;
        auto main_cfg = readKeyValueConfig(argv[1]);
        auto cfg = expandConfig(main_cfg, argv[1], component_paths);

        TelescopeConfig telescope_cfg = buildTelescopeConfig(cfg);
        CameraConfig camera_cfg = buildCameraConfig(cfg);
        CameraGeometry camera = buildCameraGeometry(camera_cfg);
        OpticalEfficiencyConfig efficiency_cfg = buildEfficiencyConfig(cfg);
        NsbConfig nsb_cfg = buildNsbConfig(cfg);
        resolveNsbSpectralRate(nsb_cfg, efficiency_cfg, camera, telescope_cfg);

        printSection("NSB spectral rate");
        printField("config", argv[1]);
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
