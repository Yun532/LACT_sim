#include "io/CorsikaTraceEventIOInput.hpp"

#include <algorithm>
#include <cmath>

namespace lact {

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
        applyEventIORotationCenter(
            out, source_runtime_cfg.eventio_rotation_center_local_m);
    }
    return out;
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

} // namespace lact
