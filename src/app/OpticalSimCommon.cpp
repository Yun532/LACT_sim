#include "app/OpticalSimCommon.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "core/DeterministicSeed.hpp"
#include "geometry/MirrorFacetValidation.hpp"
#include "io/CorsikaTraceOutputTypes.hpp"
#include "io/MirrorFacetCsvReader.hpp"

namespace lact {

extern const double DEG_TO_RAD = 3.14159265358979323846 / 180.0;

namespace {

std::uint64_t mixNsbSeed(std::uint64_t seed, std::uint64_t value)
{
    seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
    return seed;
}

std::mt19937_64 makeNsbEventTelescopeRng(const NsbConfig& nsb,
                                         int event_id,
                                         int telescope_id)
{
    std::uint64_t seed = nsb.seed;
    seed = mixNsbSeed(seed, static_cast<std::uint64_t>(
                                static_cast<std::int64_t>(event_id) + 0x80000000LL));
    seed = mixNsbSeed(seed, static_cast<std::uint64_t>(telescope_id + 1));
    return std::mt19937_64(seed);
}

} // namespace

std::string trim(const std::string& s) {
    auto first = std::find_if_not(s.begin(), s.end(), [](unsigned char c) {
        return std::isspace(c);
    });
    auto last = std::find_if_not(s.rbegin(), s.rend(), [](unsigned char c) {
        return std::isspace(c);
    }).base();
    if (first >= last) {
        return "";
    }
    return std::string(first, last);
}

std::string lowerCopy(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return s;
}

bool startsWith(const std::string& text, const std::string& prefix) {
    return text.rfind(prefix, 0) == 0;
}

std::string expandEnvironmentVariables(const std::string& text)
{
    std::string out;
    for (std::size_t i = 0; i < text.size();) {
        if (text[i] != '$') {
            out.push_back(text[i++]);
            continue;
        }

        if (i + 1 < text.size() && text[i + 1] == '{') {
            const std::size_t close = text.find('}', i + 2);
            if (close == std::string::npos) {
                out.push_back(text[i++]);
                continue;
            }
            const std::string name = text.substr(i + 2, close - i - 2);
            const char* value = std::getenv(name.c_str());
            if (!value) {
                throw std::runtime_error("environment variable not set: " + name);
            }
            out += value;
            i = close + 1;
            continue;
        }

        std::size_t j = i + 1;
        while (j < text.size()) {
            const unsigned char c = static_cast<unsigned char>(text[j]);
            if (!(std::isalnum(c) || text[j] == '_')) {
                break;
            }
            ++j;
        }
        if (j == i + 1) {
            out.push_back(text[i++]);
            continue;
        }
        const std::string name = text.substr(i + 1, j - i - 1);
        const char* value = std::getenv(name.c_str());
        if (!value) {
            throw std::runtime_error("environment variable not set: " + name);
        }
        out += value;
        i = j;
    }
    return out;
}

std::string parentDirectory(const std::string& path) {
    auto pos = path.find_last_of("/\\");
    if (pos == std::string::npos) {
        return ".";
    }
    if (pos == 0) {
        return path.substr(0, 1);
    }
    return path.substr(0, pos);
}

bool isAbsolutePath(const std::string& path) {
    return !path.empty() && path.front() == '/';
}

std::string resolveRelativePath(const std::string& base_config_path,
                                const std::string& path)
{
    if (path.empty() || isAbsolutePath(path)) {
        return path;
    }
    std::filesystem::path base(base_config_path);
    std::filesystem::path resolved = base.parent_path() / path;
    return std::filesystem::absolute(resolved).lexically_normal().string();
}

std::map<std::string, std::string> readKeyValueConfig(const std::string& path) {
    std::ifstream ifs(path);
    if (!ifs) {
        throw std::runtime_error("failed to open config: " + path);
    }

    std::map<std::string, std::string> values;
    std::string line;
    int line_no = 0;
    while (std::getline(ifs, line)) {
        ++line_no;
        auto comment = line.find('#');
        if (comment != std::string::npos) {
            line = line.substr(0, comment);
        }
        line = trim(line);
        if (line.empty()) {
            continue;
        }

        auto eq = line.find('=');
        if (eq == std::string::npos) {
            std::ostringstream oss;
            oss << "config line " << line_no << " missing '='";
            throw std::runtime_error(oss.str());
        }

        std::string key = lowerCopy(trim(line.substr(0, eq)));
        std::string value = trim(line.substr(eq + 1));
        if (key.empty()) {
            std::ostringstream oss;
            oss << "config line " << line_no << " has empty key";
            throw std::runtime_error(oss.str());
        }
        values[key] = value;
    }

    return values;
}

ConfigCommandLine parseConfigCommandLine(int argc, char** argv)
{
    ConfigCommandLine out;
    bool parse_options = true;
    auto append_override = [&](const std::string& text) {
        const auto equal = text.find('=');
        if (equal == std::string::npos) {
            throw std::runtime_error(
                "configuration override must be key=value: " + text);
        }
        const std::string key = lowerCopy(trim(text.substr(0, equal)));
        if (key.empty()) {
            throw std::runtime_error(
                "configuration override has an empty key: " + text);
        }
        out.overrides.emplace_back(key, trim(text.substr(equal + 1)));
    };

    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        if (parse_options && argument == "--") {
            parse_options = false;
        } else if (parse_options &&
                   (argument == "--help" || argument == "-h")) {
            out.help = true;
        } else if (parse_options &&
                   (argument == "-C" || argument == "--set")) {
            if (i + 1 >= argc) {
                throw std::runtime_error(argument + " requires key=value");
            }
            append_override(argv[++i]);
        } else if (parse_options && startsWith(argument, "-C") &&
                   argument.size() > 2) {
            append_override(argument.substr(2));
        } else if (parse_options && startsWith(argument, "--set=")) {
            append_override(argument.substr(6));
        } else if (parse_options && !argument.empty() && argument[0] == '-') {
            throw std::runtime_error("unknown command-line option: " + argument);
        } else {
            out.positional.push_back(argument);
        }
    }
    return out;
}

void applyConfigOverrides(
    std::map<std::string, std::string>& cfg,
    const std::vector<std::pair<std::string, std::string>>& overrides)
{
    for (const auto& [key, value] : overrides) {
        cfg[key] = value;
    }
}

std::string scopedComponentKey(const std::string& key, const std::string& prefix) {
    if (startsWith(key, "telescope.") ||
        startsWith(key, "mirror.") || startsWith(key, "source.") ||
        startsWith(key, "output.") || startsWith(key, "camera.") ||
        startsWith(key, "sipm.") ||
        startsWith(key, "electronics.") ||
        startsWith(key, "efficiency.") ||
        startsWith(key, "atmosphere.") ||
        startsWith(key, "nsb.") ||
        startsWith(key, "trigger.") ||
        startsWith(key, "error.") ||
        startsWith(key, "obstruction.") ||
        startsWith(key, "dish.") ||
        startsWith(key, "facet.")) {
        return key;
    }
    return prefix + key;
}

Vec3 TelescopeFrame::rotateVector(const Vec3& local) const {
    return x_axis * local.x + y_axis * local.y + z_axis * local.z;
}

Vec3 TelescopeFrame::pointToGlobal(const Vec3& local) const {
    return origin + rotateVector(local);
}

Vec3 TelescopeFrame::rotateVectorToLocal(const Vec3& global) const {
    return {global.dot(x_axis), global.dot(y_axis), global.dot(z_axis)};
}

Vec3 TelescopeFrame::pointToLocal(const Vec3& global) const {
    return rotateVectorToLocal(global - origin);
}

bool isIncludeConfigKey(const std::string& key) {
    return key == "mirror.config" || key == "source.config" ||
           key == "output.config" || key == "camera.config" ||
           key == "sipm.config" ||
           key == "electronics.config" ||
           key == "efficiency.config" ||
           key == "atmosphere.config" ||
           key == "nsb.config" ||
           key == "trigger.config" ||
           key == "error.config" ||
           key == "obstruction.config";
}

void mergeTelescopeConfig(std::map<std::string, std::string>& dst,
                          ComponentConfigPaths& paths,
                          const std::map<std::string, std::string>& main_cfg,
                          const std::string& main_config_path)
{
    auto it = main_cfg.find("telescope.config");
    if (it == main_cfg.end() || trim(it->second).empty()) {
        return;
    }

    const std::string path = resolveRelativePath(main_config_path, trim(it->second));
    auto telescope_cfg = readKeyValueConfig(path);
    for (const auto& kv : telescope_cfg) {
        const std::string key = scopedComponentKey(kv.first, "telescope.");
        std::string value = kv.second;
        if (isIncludeConfigKey(key)) {
            value = resolveRelativePath(path, value);
        }
        dst[key] = value;
    }
    paths.telescope = path;
}

void mergeComponentConfig(std::map<std::string, std::string>& dst,
                          ComponentConfigPaths& paths,
                          const std::map<std::string, std::string>& assembly_cfg,
                          const std::string& main_config_path,
                          const std::string& include_key,
                          const std::string& prefix)
{
    auto it = assembly_cfg.find(include_key);
    if (it == assembly_cfg.end() || trim(it->second).empty()) {
        return;
    }

    const std::string path = resolveRelativePath(main_config_path, trim(it->second));
    auto component_cfg = readKeyValueConfig(path);
    for (const auto& kv : component_cfg) {
        const std::string scoped = scopedComponentKey(kv.first, prefix);
        std::string value = kv.second;
        if (scoped == "mirror.series_csv_path" || scoped == "mirror.series_csv_pattern" ||
            scoped == "efficiency.mirror_reflectivity_scale_csv") {
            value = resolveRelativePath(path, value);
        }
        if (scoped == "error.structural_deformation_config") {
            value = resolveRelativePath(path, value);
        }
        if (scoped == "obstruction.mask_csv" ||
            scoped == "obstruction.primitives_csv") {
            value = resolveRelativePath(path, value);
        }
        if (scoped == "atmosphere.tau_table") {
            value = resolveRelativePath(path, value);
        }
        if (scoped == "nsb.spectrum_csv") {
            value = resolveRelativePath(path, value);
        }
        if (scoped == "electronics.microcell.device" ||
            scoped == "electronics.single_pe.csv" ||
            scoped == "electronics.single_pe.csv_path" ||
            scoped == "electronics.single_pe.template" ||
            scoped == "electronics.single_pe.charge_fluctuation.csv_path" ||
            scoped == "electronics.single_pe.charge_fluctuation.samples") {
            value = resolveRelativePath(path, value);
        }
        dst[scoped] = value;
    }

    if (include_key == "mirror.config") {
        paths.mirror = path;
    } else if (include_key == "source.config") {
        paths.source = path;
    } else if (include_key == "output.config") {
        paths.output = path;
    } else if (include_key == "camera.config") {
        paths.camera = path;
    } else if (include_key == "sipm.config") {
        paths.sipm = path;
    } else if (include_key == "electronics.config") {
        paths.electronics = path;
    } else if (include_key == "efficiency.config") {
        paths.efficiency = path;
    } else if (include_key == "atmosphere.config") {
        paths.atmosphere = path;
    } else if (include_key == "nsb.config") {
        paths.nsb = path;
    } else if (include_key == "trigger.config") {
        paths.trigger = path;
    } else if (include_key == "error.config") {
        paths.error = path;
    } else if (include_key == "obstruction.config") {
        paths.obstruction = path;
    }
}

std::map<std::string, std::string> expandConfig(const std::map<std::string, std::string>& main_cfg,
                                                const std::string& main_config_path,
                                                ComponentConfigPaths& paths)
{
    std::map<std::string, std::string> expanded;
    std::map<std::string, std::string> assembly_cfg;
    mergeTelescopeConfig(assembly_cfg, paths, main_cfg, main_config_path);

    // Values in the main file intentionally override telescope defaults.
    for (const auto& kv : main_cfg) {
        assembly_cfg[kv.first] = kv.second;
    }

    mergeComponentConfig(expanded, paths, assembly_cfg, main_config_path,
                         "mirror.config", "mirror.");
    mergeComponentConfig(expanded, paths, assembly_cfg, main_config_path,
                         "source.config", "source.");
    mergeComponentConfig(expanded, paths, assembly_cfg, main_config_path,
                         "output.config", "output.");
    mergeComponentConfig(expanded, paths, assembly_cfg, main_config_path,
                         "camera.config", "camera.");
    mergeComponentConfig(expanded, paths, assembly_cfg, main_config_path,
                         "sipm.config", "sipm.");
    mergeComponentConfig(expanded, paths, assembly_cfg, main_config_path,
                         "electronics.config", "electronics.");
    mergeComponentConfig(expanded, paths, assembly_cfg, main_config_path,
                         "efficiency.config", "efficiency.");
    mergeComponentConfig(expanded, paths, assembly_cfg, main_config_path,
                         "atmosphere.config", "atmosphere.");
    mergeComponentConfig(expanded, paths, assembly_cfg, main_config_path,
                         "nsb.config", "nsb.");
    mergeComponentConfig(expanded, paths, assembly_cfg, main_config_path,
                         "trigger.config", "trigger.");
    mergeComponentConfig(expanded, paths, assembly_cfg, main_config_path,
                         "error.config", "error.");
    mergeComponentConfig(expanded, paths, assembly_cfg, main_config_path,
                         "obstruction.config", "obstruction.");

    // Assembly values intentionally override component defaults.
    for (const auto& kv : assembly_cfg) {
        expanded[kv.first] = kv.second;
    }
    return expanded;
}

std::string getString(const std::map<std::string, std::string>& cfg,
                      const std::string& key,
                      const std::string& fallback)
{
    auto it = cfg.find(key);
    return expandEnvironmentVariables(it == cfg.end() ? fallback : it->second);
}

double getDouble(const std::map<std::string, std::string>& cfg,
                 const std::string& key,
                 double fallback)
{
    auto it = cfg.find(key);
    if (it == cfg.end()) {
        return fallback;
    }
    std::size_t pos = 0;
    double value = std::stod(it->second, &pos);
    if (pos != it->second.size()) {
        throw std::runtime_error("invalid numeric config value for " + key + ": " + it->second);
    }
    return value;
}

int getInt(const std::map<std::string, std::string>& cfg,
           const std::string& key,
           int fallback)
{
    auto it = cfg.find(key);
    if (it == cfg.end()) {
        return fallback;
    }
    std::size_t pos = 0;
    int value = std::stoi(it->second, &pos);
    if (pos != it->second.size()) {
        throw std::runtime_error("invalid integer config value for " + key + ": " + it->second);
    }
    return value;
}

std::uint64_t getUInt64(const std::map<std::string, std::string>& cfg,
                        const std::string& key,
                        std::uint64_t fallback)
{
    auto it = cfg.find(key);
    if (it == cfg.end()) {
        return fallback;
    }
    std::size_t pos = 0;
    auto value = static_cast<std::uint64_t>(std::stoull(it->second, &pos));
    if (pos != it->second.size()) {
        throw std::runtime_error("invalid integer config value for " + key + ": " + it->second);
    }
    return value;
}

bool getBool(const std::map<std::string, std::string>& cfg,
             const std::string& key,
             bool fallback)
{
    auto it = cfg.find(key);
    if (it == cfg.end()) {
        return fallback;
    }
    std::string value = lowerCopy(trim(it->second));
    if (value == "1" || value == "true" || value == "yes" || value == "on") {
        return true;
    }
    if (value == "0" || value == "false" || value == "no" || value == "off") {
        return false;
    }
    throw std::runtime_error("invalid boolean config value for " + key + ": " + it->second);
}

Vec3 parseVec3(const std::string& text, const std::string& key) {
    std::stringstream ss(text);
    std::string cell;
    std::vector<double> values;
    while (std::getline(ss, cell, ',')) {
        cell = trim(cell);
        if (cell.empty()) {
            throw std::runtime_error("empty component in vector config value for " + key);
        }
        std::size_t pos = 0;
        double value = std::stod(cell, &pos);
        if (pos != cell.size()) {
            throw std::runtime_error("invalid vector component in " + key + ": " + cell);
        }
        values.push_back(value);
    }
    if (values.size() != 3) {
        throw std::runtime_error("vector config value for " + key + " must have 3 components");
    }
    return {values[0], values[1], values[2]};
}

Vec3 getVec3(const std::map<std::string, std::string>& cfg,
             const std::string& key,
             const Vec3& fallback)
{
    auto it = cfg.find(key);
    return it == cfg.end() ? fallback : parseVec3(it->second, key);
}

DishType parseDishType(const std::string& text) {
    std::string s = lowerCopy(trim(text));
    if (s == "daviescotton" || s == "davies-cotton" || s == "dc") {
        return DishType::DaviesCotton;
    }
    if (s == "parabolic" || s == "paraboloid") {
        return DishType::Parabolic;
    }
    throw std::runtime_error("unsupported dish.type: " + text);
}

SurfaceType parseSurfaceType(const std::string& text) {
    std::string s = lowerCopy(trim(text));
    if (s == "spherical") return SurfaceType::Spherical;
    if (s == "parabolic") return SurfaceType::Parabolic;
    if (s == "planar") return SurfaceType::Planar;
    if (s == "polynomial") return SurfaceType::Polynomial;
    throw std::runtime_error("unsupported facet.surface_type: " + text);
}

ApertureShape parseApertureShape(const std::string& text) {
    std::string s = lowerCopy(trim(text));
    if (s == "circular") return ApertureShape::Circular;
    if (s == "hexagon") return ApertureShape::Hexagon;
    if (s == "square") return ApertureShape::Square;
    throw std::runtime_error("unsupported facet.aperture_shape: " + text);
}

PixelShape parsePixelShape(const std::string& text) {
    std::string s = lowerCopy(trim(text));
    if (s == "circular" || s == "circle") return PixelShape::Circular;
    if (s == "hexagonal" || s == "hexagon" || s == "hex") return PixelShape::Hexagonal;
    if (s == "square") return PixelShape::Square;
    throw std::runtime_error("unsupported camera.pixel_shape: " + text);
}

std::string pixelShapeName(PixelShape shape) {
    switch (shape) {
        case PixelShape::Circular: return "Circular";
        case PixelShape::Hexagonal: return "Hexagonal";
        case PixelShape::Square: return "Square";
    }
    return "Unknown";
}

SyntheticMode parseSyntheticMode(const std::string& text) {
    std::string s = lowerCopy(trim(text));
    if (s == "parallelbeam" || s == "parallel") {
        return SyntheticMode::ParallelBeam;
    }
    if (s == "pointsource" || s == "point") {
        return SyntheticMode::PointSource;
    }
    throw std::runtime_error("unsupported source.mode: " + text);
}

bool isPhotonCsvMode(const std::string& text) {
    std::string s = lowerCopy(trim(text));
    return s == "photoncsv" || s == "photon_csv" || s == "csv" || s == "file";
}

bool isEventIOMode(const std::string& text) {
    std::string s = lowerCopy(trim(text));
    return s == "eventio" || s == "corsika" || s == "corsikaeventio" ||
           s == "corsika_eventio" || s == "iact";
}

std::string vec3ToString(const Vec3& v) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6)
        << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    return oss.str();
}

std::string sourceModeName(SyntheticMode mode) {
    switch (mode) {
    case SyntheticMode::ParallelBeam:
        return "ParallelBeam";
    case SyntheticMode::PointSource:
        return "PointSource";
    }
    return "Unknown";
}

TelescopeFrame buildTelescopeFrame(const TelescopeConfig& telescope)
{
    // Canonical LACT telescope structure frame in global NWU coordinates:
    //   global +X = North, +Y = West, +Z = Up
    //   local  +z = boresight toward the sky
    //   local  +x = mirror/obstruction CSV +x as seen camera -> mirror
    //               (West for a telescope pointing North)
    //   local  +y = increasing elevation / sky-up
    // Pointing azimuth is zero at North and increases clockwise toward East.
    // Keep this identical to the CORSIKA adapter so synthetic and EventIO
    // photons enter one optical geometry with the same physical axes.
    const double az = telescope.pointing_az_deg * DEG_TO_RAD;
    const double el = telescope.pointing_el_deg * DEG_TO_RAD;
    const double sin_el = std::sin(el);
    const double cos_el = std::cos(el);
    const double sin_az = std::sin(az);
    const double cos_az = std::cos(az);

    TelescopeFrame frame;
    frame.origin = telescope.position_m;
    frame.x_axis = Vec3{sin_az, cos_az, 0.0}.normalized();
    frame.y_axis = Vec3{-sin_el * cos_az, sin_el * sin_az, cos_el}.normalized();
    frame.z_axis = Vec3{cos_el * cos_az, -cos_el * sin_az, sin_el}.normalized();
    return frame;
}

TelescopeFrame buildCorsikaNwuTelescopeFrame(const TelescopeConfig& telescope)
{
    // CORSIKA/sim_telarray input adapter in magnetic North-West-Up
    // coordinates.  The returned axes are the canonical telescope-local
    // optical axes used by the mirror, obstruction, output plane, and camera:
    //   local +z: boresight toward the sky
    //   local +y: increasing elevation (toward zenith / sky-up)
    //   local +x: completes the right-handed frame; opposite increasing
    //             azimuth (West at az=0)
    // This is right-handed: x cross y = z.
    const double az = telescope.pointing_az_deg * DEG_TO_RAD;
    const double el = telescope.pointing_el_deg * DEG_TO_RAD;
    const double sin_el = std::sin(el);
    const double cos_el = std::cos(el);
    const double sin_az = std::sin(az);
    const double cos_az = std::cos(az);

    TelescopeFrame frame;
    frame.origin = telescope.position_m;
    frame.x_axis = Vec3{sin_az, cos_az, 0.0}.normalized();
    frame.y_axis = Vec3{-sin_el * cos_az, sin_el * sin_az, cos_el}.normalized();
    frame.z_axis = Vec3{cos_el * cos_az, -cos_el * sin_az, sin_el}.normalized();
    return frame;
}

TelescopeFrame buildEnuEastTelescopeFrame(const TelescopeConfig& telescope)
{
    // Explicit East-North-Up input frame:
    //   global +x: East, +y: North, +z: Up
    //   pointing az=0: East, az=90: North
    //
    // Keep the same canonical right-handed telescope-local optical axes as
    // the CORSIKA adapter: local +y is sky-up, local +x is the transverse
    // horizontal axis, and local +z is the boresight.
    TelescopeConfig nwu_telescope = telescope;
    nwu_telescope.position_m = {
        telescope.position_m.y, -telescope.position_m.x, telescope.position_m.z};
    nwu_telescope.pointing_az_deg = 90.0 - telescope.pointing_az_deg;
    const TelescopeFrame nwu = buildCorsikaNwuTelescopeFrame(nwu_telescope);
    TelescopeFrame frame;
    frame.origin = telescope.position_m;
    frame.x_axis = {-nwu.x_axis.y, nwu.x_axis.x, nwu.x_axis.z};
    frame.y_axis = {-nwu.y_axis.y, nwu.y_axis.x, nwu.y_axis.z};
    frame.z_axis = {-nwu.z_axis.y, nwu.z_axis.x, nwu.z_axis.z};
    return frame;
}

std::string normalizeSourceCoordinateFrame(const std::string& frame_name)
{
    const std::string frame = lowerCopy(trim(frame_name));
    if (frame.empty() || frame == "telescope_local" || frame == "local" ||
        frame == "optical_local") {
        return "telescope_local";
    }
    if (frame == "corsika_nwu_relative" || frame == "corsika_iact" ||
        frame == "corsika" || frame == "simtelarray") {
        return "corsika_nwu_relative";
    }
    if (frame == "corsika_nwu_global" || frame == "corsika_global") {
        return "corsika_nwu_global";
    }
    if (frame == "enu_east_relative" || frame == "enu_relative" ||
        frame == "east_north_up_relative" || frame == "east_start_relative") {
        return "enu_east_relative";
    }
    if (frame == "enu_east_global" || frame == "enu_global" ||
        frame == "east_north_up_global" || frame == "east_start_global") {
        return "enu_east_global";
    }
    if (frame == "lact_generic_global" || frame == "generic_global" ||
        frame == "array_global" || frame == "global") {
        return "lact_generic_global";
    }
    throw std::runtime_error(
        "unsupported source.coordinate_frame: " + frame_name +
        "; expected telescope_local, corsika_nwu_relative, "
        "corsika_nwu_global, enu_east_relative, enu_east_global, "
        "or lact_generic_global");
}

std::string sourceCoordinateFrameDescription(const std::string& frame_name)
{
    const std::string frame = normalizeSourceCoordinateFrame(frame_name);
    if (frame == "telescope_local") {
        return "existing telescope-local optical coordinates; axis definitions are unchanged";
    }
    if (frame == "corsika_nwu_relative") {
        return "CORSIKA NWU; positions are relative to the selected telescope";
    }
    if (frame == "corsika_nwu_global") {
        return "CORSIKA NWU absolute array positions; telescope.position_m is subtracted";
    }
    if (frame == "enu_east_relative") {
        return "ENU (+x East, +y North, +z Up), relative to the telescope; "
               "pointing azimuth starts at East and increases toward North";
    }
    if (frame == "enu_east_global") {
        return "ENU (+x East, +y North, +z Up) absolute array positions; "
               "telescope.position_m is subtracted; pointing azimuth starts at "
               "East and increases toward North";
    }
    return "LACT global NWU (+x North, +y West, +z Up); azimuth is 0 at North and increases toward East";
}

PhotonBunch transformBunchToTelescopeLocal(const PhotonBunch& input,
                                           const TelescopeConfig& telescope,
                                           const std::string& frame_name)
{
    PhotonBunch out = input;
    const std::string frame_name_normalized = normalizeSourceCoordinateFrame(frame_name);
    if (frame_name_normalized == "telescope_local") {
        out.photon.normalizeDirection();
        return out;
    }

    if (frame_name_normalized == "corsika_nwu_relative") {
        const TelescopeFrame frame = buildCorsikaNwuTelescopeFrame(telescope);
        out.photon.pos = frame.rotateVectorToLocal(input.photon.pos);
        out.photon.dir = frame.rotateVectorToLocal(input.photon.dir).normalized();
        return out;
    }

    if (frame_name_normalized == "corsika_nwu_global") {
        const TelescopeFrame frame = buildCorsikaNwuTelescopeFrame(telescope);
        out.photon.pos = frame.pointToLocal(input.photon.pos);
        out.photon.dir = frame.rotateVectorToLocal(input.photon.dir).normalized();
        return out;
    }

    if (frame_name_normalized == "enu_east_relative") {
        TelescopeConfig nwu_telescope = telescope;
        nwu_telescope.position_m = {
            telescope.position_m.y, -telescope.position_m.x, telescope.position_m.z};
        nwu_telescope.pointing_az_deg = 90.0 - telescope.pointing_az_deg;
        const TelescopeFrame frame = buildCorsikaNwuTelescopeFrame(nwu_telescope);
        const Vec3 nwu_pos{input.photon.pos.y, -input.photon.pos.x,
                           input.photon.pos.z};
        const Vec3 nwu_dir{input.photon.dir.y, -input.photon.dir.x,
                           input.photon.dir.z};
        out.photon.pos = frame.rotateVectorToLocal(nwu_pos);
        out.photon.dir = frame.rotateVectorToLocal(nwu_dir).normalized();
        return out;
    }

    if (frame_name_normalized == "enu_east_global") {
        TelescopeConfig nwu_telescope = telescope;
        nwu_telescope.position_m = {
            telescope.position_m.y, -telescope.position_m.x, telescope.position_m.z};
        nwu_telescope.pointing_az_deg = 90.0 - telescope.pointing_az_deg;
        const TelescopeFrame frame = buildCorsikaNwuTelescopeFrame(nwu_telescope);
        const Vec3 nwu_pos{input.photon.pos.y, -input.photon.pos.x,
                           input.photon.pos.z};
        const Vec3 nwu_dir{input.photon.dir.y, -input.photon.dir.x,
                           input.photon.dir.z};
        out.photon.pos = frame.pointToLocal(nwu_pos);
        out.photon.dir = frame.rotateVectorToLocal(nwu_dir).normalized();
        return out;
    }

    const TelescopeFrame frame = buildTelescopeFrame(telescope);
    out.photon.pos = frame.pointToLocal(input.photon.pos);
    out.photon.dir = frame.rotateVectorToLocal(input.photon.dir).normalized();
    return out;
}

void applyEventIOReferenceZOffset(PhotonBunch& bunch,
                                  double eventio_reference_z_m)
{
    if (!std::isfinite(eventio_reference_z_m)) {
        throw std::runtime_error(
            "source.eventio_reference_z_m must be finite");
    }
    // The user config remains one scalar. Internally it is the
    // telescope-local translation (0,0,z); the physical x/y carried by the
    // bunch are deliberately left unchanged.
    bunch.photon.pos.z += eventio_reference_z_m;
}

Vec3 sourceDirectionInWorld(const PhotonBunch& input,
                            const TelescopeConfig& telescope,
                            const std::string& frame_name)
{
    const std::string frame = normalizeSourceCoordinateFrame(frame_name);
    if (frame == "telescope_local") {
        return buildTelescopeFrame(telescope)
            .rotateVector(input.photon.dir)
            .normalized();
    }
    // CORSIKA and LACT global frames already carry physical NWU directions.
    // ENU input has the same Up component, which is the component used by the
    // horizontally symmetric atmosphere model.
    return input.photon.dir.normalized();
}

void applyTelescopeFrame(std::vector<MirrorFacet>& facets,
                         OutputPlane& plane,
                         const TelescopeFrame& frame)
{
    // 镜片和输出面都先在本地坐标里定义，再整体转到全局坐标。
    // 这样 source / mirror / output 的相对几何关系在不同 pointing 下保持一致。
    for (auto& facet : facets) {
        Vec3 local_u;
        Vec3 local_v;
        facet.apertureFrame(local_u, local_v);
        facet.center = frame.pointToGlobal(facet.center);
        facet.normal = frame.rotateVector(facet.normal).normalized();
        facet.aperture_u_axis = frame.rotateVector(local_u).normalized();
        facet.aperture_v_axis = frame.rotateVector(local_v).normalized();
    }

    plane.point = frame.pointToGlobal(plane.point);
    plane.normal = frame.rotateVector(plane.normal).normalized();
    plane.u_axis = frame.rotateVector(plane.u_axis).normalized();
    plane.v_axis = frame.rotateVector(plane.v_axis).normalized();
}

void applyTelescopeFrame(Photon& photon, const TelescopeFrame& frame)
{
    // 光子源同样先在本地望远镜坐标里采样，再整体转到全局。
    // 因此平行光撒点圆盘和点源 target plane 会随 pointing 一起转动。
    photon.pos = frame.pointToGlobal(photon.pos);
    photon.dir = frame.rotateVector(photon.dir).normalized();
}

double elapsedSeconds(std::chrono::steady_clock::time_point start,
                      std::chrono::steady_clock::time_point stop)
{
    return std::chrono::duration<double>(stop - start).count();
}

void printSection(const std::string& title) {
    std::cout << "\n[" << title << "]\n";
}

void printField(const std::string& label, const std::string& value) {
    std::cout << "  " << std::left << std::setw(24) << label
              << ": " << value << "\n" << std::right;
}

std::string doubleToString(double value, int precision) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << value;
    return oss.str();
}

std::string intToString(std::uint64_t value) {
    std::ostringstream oss;
    oss << value;
    return oss.str();
}

bool parseDoubleStrict(const std::string& text, double& out) {
    std::size_t pos = 0;
    try {
        out = std::stod(text, &pos);
    } catch (...) {
        return false;
    }
    return pos == text.size();
}

bool isDisabledText(const std::string& text) {
    std::string value = lowerCopy(trim(text));
    return value.empty() || value == "none" || value == "off" ||
           value == "false" || value == "no";
}

std::string factorDescription(const EfficiencyFactorConfig& factor) {
    if (!factor.enabled) {
        return "not set -> 1";
    }
    if (factor.use_curve) {
        return "curve: " + factor.csv_path;
    }
    return "constant: " + doubleToString(factor.constant);
}

EfficiencyFactorConfig parseEfficiencyFactor(const std::map<std::string, std::string>& cfg,
                                             const std::string& base_key)
{
    EfficiencyFactorConfig factor;
    const std::string combined = getString(cfg, base_key, "");
    const std::string legacy_csv = getString(cfg, base_key + "_csv", "");
    const std::string constant_text = getString(cfg, base_key + "_constant", "");

    std::string value;
    if (!combined.empty()) {
        value = combined;
    } else if (!constant_text.empty()) {
        value = constant_text;
    } else if (!legacy_csv.empty()) {
        value = legacy_csv;
    }

    if (isDisabledText(value)) {
        return factor;
    }

    double constant = 1.0;
    if (parseDoubleStrict(value, constant)) {
        if (!std::isfinite(constant) || constant < 0.0 || constant > 1.0) {
            throw std::runtime_error(
                base_key + " constant must be finite and in [0,1]");
        }
        factor.enabled = true;
        factor.use_curve = false;
        factor.constant = constant;
        return factor;
    }

    factor.enabled = true;
    factor.use_curve = true;
    factor.csv_path = value;
    return factor;
}

bool hasEfficiencyFactorConfig(const std::map<std::string, std::string>& cfg,
                               const std::string& base_key)
{
    return cfg.find(base_key) != cfg.end() ||
           cfg.find(base_key + "_csv") != cfg.end() ||
           cfg.find(base_key + "_constant") != cfg.end();
}

EfficiencyFactorConfig parseFirstEfficiencyFactor(
    const std::map<std::string, std::string>& cfg,
    const std::vector<std::string>& base_keys)
{
    for (const auto& base_key : base_keys) {
        if (hasEfficiencyFactorConfig(cfg, base_key)) {
            return parseEfficiencyFactor(cfg, base_key);
        }
    }
    return {};
}

std::vector<double> parseDoubleList(const std::string& text, const std::string& key) {
    std::stringstream ss(text);
    std::string cell;
    std::vector<double> values;
    while (std::getline(ss, cell, ',')) {
        cell = trim(cell);
        if (cell.empty()) {
            throw std::runtime_error("empty component in list config value for " + key);
        }
        std::size_t pos = 0;
        double value = std::stod(cell, &pos);
        if (pos != cell.size()) {
            throw std::runtime_error("invalid list component in " + key + ": " + cell);
        }
        values.push_back(value);
    }
    if (values.empty()) {
        throw std::runtime_error("list config value for " + key + " is empty");
    }
    return values;
}

std::vector<int> parseIntList(const std::string& text, const std::string& key) {
    std::stringstream ss(text);
    std::string cell;
    std::vector<int> values;
    while (std::getline(ss, cell, ',')) {
        cell = trim(cell);
        if (cell.empty()) {
            throw std::runtime_error("empty component in list config value for " + key);
        }
        std::size_t pos = 0;
        const int value = std::stoi(cell, &pos);
        if (pos != cell.size()) {
            throw std::runtime_error("invalid list component in " + key + ": " + cell);
        }
        values.push_back(value);
    }
    if (values.empty()) {
        throw std::runtime_error("list config value for " + key + " is empty");
    }
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    return values;
}

std::string formatAngleToken(double angle_deg) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6) << angle_deg;
    std::string token = oss.str();
    while (!token.empty() && token.back() == '0') {
        token.pop_back();
    }
    if (!token.empty() && token.back() == '.') {
        token.pop_back();
    }
    if (token.empty()) {
        token = "0";
    }
    return token;
}

std::string expandAnglePattern(const std::string& pattern, double angle_deg) {
    const std::string token = "{angle}";
    std::string path = pattern;
    auto pos = path.find(token);
    if (pos == std::string::npos) {
        throw std::runtime_error("mirror.series_csv_pattern must contain {angle}");
    }
    path.replace(pos, token.size(), formatAngleToken(angle_deg));
    return path;
}

struct Raw1229FacetState {
    int cell_index = 0;
    int ring_index = 0;
    Vec3 center;
    Vec3 normal;
};

struct SeriesFacetState {
    int id = -1;
    int cell_index = -1;
    int ring_index = -1;
    Vec3 center;
    Vec3 normal;
    bool has_surface_type = false;
    SurfaceType surface_type = SurfaceType::Spherical;
    bool has_radius_of_curvature = false;
    double radius_of_curvature = 16.0;
    bool has_aperture_shape = false;
    ApertureShape aperture_shape = ApertureShape::Hexagon;
    bool has_size1 = false;
    double size1 = 0.8;
    bool has_size2 = false;
    double size2 = 0.0;
    bool has_aperture_rotation = false;
    double aperture_rotation_rad = 0.0;
    bool has_reflectivity_scale = false;
    double reflectivity_scale = 1.0;
    bool has_roughness_sigma_rad = false;
    double roughness_sigma_rad = 0.0;
    bool has_misalign_sigma_rad = false;
    double misalign_sigma_rad = 0.0;
};

std::vector<Raw1229FacetState> readRaw1229FacetStates(const std::string& path, bool swap_xy) {
    std::ifstream ifs(path);
    if (!ifs) {
        throw std::runtime_error("failed to open raw 1229 mirror file: " + path);
    }

    std::vector<Raw1229FacetState> states;
    std::string line;
    int line_no = 0;
    while (std::getline(ifs, line)) {
        ++line_no;
        line = trim(line);
        if (line.empty()) {
            continue;
        }

        std::stringstream ss(line);
        std::string cell;
        std::vector<double> values;
        while (std::getline(ss, cell, ',')) {
            cell = trim(cell);
            if (cell.empty()) {
                throw std::runtime_error("raw 1229 file has empty cell at line " +
                                         std::to_string(line_no));
            }
            values.push_back(std::stod(cell));
        }

        if (values.size() != 8) {
            throw std::runtime_error("raw 1229 file line " + std::to_string(line_no) +
                                     " expected 8 columns, got " +
                                     std::to_string(values.size()));
        }

        Raw1229FacetState state;
        state.cell_index = static_cast<int>(std::lround(values[0]));
        state.ring_index = static_cast<int>(std::lround(values[1]));
        state.center = {values[2] * 0.001, values[3] * 0.001, values[4] * 0.001};
        Vec3 curvature_center = {values[5] * 0.001, values[6] * 0.001, values[7] * 0.001};

        if (swap_xy) {
            std::swap(state.center.x, state.center.y);
            std::swap(curvature_center.x, curvature_center.y);
        }

        state.normal = (curvature_center - state.center).normalized();
        states.push_back(state);
    }

    if (states.empty()) {
        throw std::runtime_error("raw 1229 mirror file is empty: " + path);
    }
    return states;
}

std::vector<std::string> splitCsvCells(const std::string& line) {
    std::stringstream ss(line);
    std::string cell;
    std::vector<std::string> cells;
    while (std::getline(ss, cell, ',')) {
        cells.push_back(trim(cell));
    }
    return cells;
}

int headerIndex(const std::map<std::string, int>& header, const std::string& key)
{
    auto it = header.find(lowerCopy(key));
    return it == header.end() ? -1 : it->second;
}

std::string getOptionalCell(const std::vector<std::string>& cells,
                            const std::map<std::string, int>& header,
                            const std::string& key)
{
    int idx = headerIndex(header, key);
    if (idx < 0 || static_cast<std::size_t>(idx) >= cells.size()) {
        return "";
    }
    return trim(cells[idx]);
}

std::string getRequiredCell(const std::vector<std::string>& cells,
                            const std::map<std::string, int>& header,
                            const std::string& key)
{
    std::string value = getOptionalCell(cells, header, key);
    if (value.empty()) {
        throw std::runtime_error("missing required CSV column: " + key);
    }
    return value;
}

namespace {
void validateCameraGeometry(const CameraGeometry& camera,
                            const std::string& description);
}

CameraGeometry readCameraCsv(const std::string& path)
{
    std::ifstream ifs(path);
    if (!ifs) {
        throw std::runtime_error("failed to open camera CSV: " + path);
    }

    std::string line;
    if (!std::getline(ifs, line)) {
        throw std::runtime_error("camera CSV is empty: " + path);
    }

    auto header_cells = splitCsvCells(line);
    std::map<std::string, int> header;
    for (int i = 0; i < static_cast<int>(header_cells.size()); ++i) {
        header[lowerCopy(header_cells[i])] = i;
    }

    CameraGeometry camera;
    int line_no = 1;
    while (std::getline(ifs, line)) {
        ++line_no;
        if (trim(line).empty()) {
            continue;
        }
        auto cells = splitCsvCells(line);

        CameraPixel pixel;
        const std::string id_text = getRequiredCell(cells, header, "id");
        std::size_t id_consumed = 0;
        pixel.id = std::stoi(id_text, &id_consumed);
        if (id_consumed != id_text.size() || pixel.id < 0) {
            throw std::runtime_error(
                "camera CSV line " + std::to_string(line_no) +
                " has an invalid non-negative integer id");
        }
        const std::string x_text = getOptionalCell(cells, header, "x_m").empty()
            ? getRequiredCell(cells, header, "center_x")
            : getOptionalCell(cells, header, "x_m");
        const std::string y_text = getOptionalCell(cells, header, "y_m").empty()
            ? getRequiredCell(cells, header, "center_y")
            : getOptionalCell(cells, header, "y_m");
        double center_x = 0.0;
        double center_y = 0.0;
        if (!parseDoubleStrict(x_text, center_x) ||
            !parseDoubleStrict(y_text, center_y) ||
            !std::isfinite(center_x) || !std::isfinite(center_y)) {
            throw std::runtime_error(
                "camera CSV line " + std::to_string(line_no) +
                " has a non-finite or invalid pixel center");
        }
        pixel.center = {center_x, center_y, 0.0};

        const std::string shape = getOptionalCell(cells, header, "shape");
        pixel.shape = shape.empty() ? PixelShape::Hexagonal : parsePixelShape(shape);

        const std::string size = getOptionalCell(cells, header, "size_m");
        if (size.empty()) {
            throw std::runtime_error("camera CSV line " + std::to_string(line_no) +
                                     " missing size_m");
        }
        if (!parseDoubleStrict(size, pixel.size) ||
            !std::isfinite(pixel.size) || pixel.size <= 0.0) {
            throw std::runtime_error(
                "camera CSV line " + std::to_string(line_no) +
                " requires finite size_m > 0");
        }
        camera.addPixel(pixel);
    }

    if (camera.empty()) {
        throw std::runtime_error("camera CSV has no data rows: " + path);
    }
    validateCameraGeometry(camera, "camera CSV " + path);
    return camera;
}

namespace {

std::vector<Vec3> cameraPixelPolygon(const CameraPixel& pixel)
{
    if (pixel.shape == PixelShape::Square) {
        const double h = 0.5 * pixel.size;
        return {
            {pixel.center.x - h, pixel.center.y - h, 0.0},
            {pixel.center.x + h, pixel.center.y - h, 0.0},
            {pixel.center.x + h, pixel.center.y + h, 0.0},
            {pixel.center.x - h, pixel.center.y + h, 0.0},
        };
    }
    if (pixel.shape == PixelShape::Hexagonal) {
        const double a = 0.5 * pixel.size;
        const double y1 = a / std::sqrt(3.0);
        const double y2 = 2.0 * y1;
        return {
            {pixel.center.x,     pixel.center.y - y2, 0.0},
            {pixel.center.x + a, pixel.center.y - y1, 0.0},
            {pixel.center.x + a, pixel.center.y + y1, 0.0},
            {pixel.center.x,     pixel.center.y + y2, 0.0},
            {pixel.center.x - a, pixel.center.y + y1, 0.0},
            {pixel.center.x - a, pixel.center.y - y1, 0.0},
        };
    }
    return {};
}

std::pair<double, double> polygonProjection(
    const std::vector<Vec3>& polygon,
    const Vec3& axis)
{
    double minimum = std::numeric_limits<double>::max();
    double maximum = std::numeric_limits<double>::lowest();
    for (const auto& point : polygon) {
        const double value = point.dot(axis);
        minimum = std::min(minimum, value);
        maximum = std::max(maximum, value);
    }
    return {minimum, maximum};
}

bool projectionsHaveInteriorOverlap(double min_a,
                                    double max_a,
                                    double min_b,
                                    double max_b,
                                    double tolerance)
{
    return max_a > min_b + tolerance &&
           max_b > min_a + tolerance;
}

std::vector<Vec3> polygonAxes(const std::vector<Vec3>& polygon)
{
    std::vector<Vec3> axes;
    axes.reserve(polygon.size());
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const Vec3 edge = polygon[(i + 1) % polygon.size()] - polygon[i];
        const Vec3 axis{-edge.y, edge.x, 0.0};
        if (axis.norm2() > 0.0) {
            axes.push_back(axis.normalized());
        }
    }
    return axes;
}

bool polygonInteriorsOverlap(const std::vector<Vec3>& first,
                             const std::vector<Vec3>& second,
                             double tolerance)
{
    auto axes = polygonAxes(first);
    auto second_axes = polygonAxes(second);
    axes.insert(axes.end(), second_axes.begin(), second_axes.end());
    for (const auto& axis : axes) {
        const auto [min_a, max_a] = polygonProjection(first, axis);
        const auto [min_b, max_b] = polygonProjection(second, axis);
        if (!projectionsHaveInteriorOverlap(
                min_a, max_a, min_b, max_b, tolerance)) {
            return false;
        }
    }
    return true;
}

bool circlePolygonInteriorsOverlap(const CameraPixel& circle,
                                   const std::vector<Vec3>& polygon,
                                   double tolerance)
{
    std::vector<Vec3> axes = polygonAxes(polygon);
    const Vec3 center = circle.center;
    const auto closest = std::min_element(
        polygon.begin(), polygon.end(),
        [&center](const Vec3& lhs, const Vec3& rhs) {
            return (lhs - center).norm2() < (rhs - center).norm2();
        });
    if (closest != polygon.end() && (*closest - center).norm2() > 0.0) {
        axes.push_back((*closest - center).normalized());
    }
    const double radius = 0.5 * circle.size;
    for (const auto& axis : axes) {
        const auto [min_polygon, max_polygon] =
            polygonProjection(polygon, axis);
        const double projected_center = center.dot(axis);
        if (!projectionsHaveInteriorOverlap(
                projected_center - radius,
                projected_center + radius,
                min_polygon,
                max_polygon,
                tolerance)) {
            return false;
        }
    }
    return true;
}

bool cameraPixelInteriorsOverlap(const CameraPixel& first,
                                 const CameraPixel& second)
{
    const double tolerance = std::max(
        1.0e-12, std::max(first.size, second.size) * 1.0e-10);
    const auto boundingRadius = [](const CameraPixel& pixel) {
        if (pixel.shape == PixelShape::Square) {
            return pixel.size / std::sqrt(2.0);
        }
        if (pixel.shape == PixelShape::Hexagonal) {
            return pixel.size / std::sqrt(3.0);
        }
        return 0.5 * pixel.size;
    };
    const double center_distance = (first.center - second.center).norm();
    if (center_distance >=
        boundingRadius(first) + boundingRadius(second) - tolerance) {
        return false;
    }
    const bool first_circle = first.shape == PixelShape::Circular;
    const bool second_circle = second.shape == PixelShape::Circular;
    if (first_circle && second_circle) {
        const double radius_sum = 0.5 * (first.size + second.size);
        return (first.center - second.center).norm() <
               radius_sum - tolerance;
    }
    if (first_circle) {
        return circlePolygonInteriorsOverlap(
            first, cameraPixelPolygon(second), tolerance);
    }
    if (second_circle) {
        return circlePolygonInteriorsOverlap(
            second, cameraPixelPolygon(first), tolerance);
    }
    return polygonInteriorsOverlap(
        cameraPixelPolygon(first), cameraPixelPolygon(second), tolerance);
}

void validateCameraGeometry(const CameraGeometry& camera,
                            const std::string& description)
{
    if (camera.empty()) {
        throw std::runtime_error(description + " has no pixels");
    }
    std::set<int> ids;
    const auto& pixels = camera.pixels();
    for (const auto& pixel : pixels) {
        if (pixel.id < 0 || !ids.insert(pixel.id).second) {
            throw std::runtime_error(
                description + " has a negative or duplicate pixel id " +
                std::to_string(pixel.id));
        }
        if (!std::isfinite(pixel.center.x) ||
            !std::isfinite(pixel.center.y) ||
            !std::isfinite(pixel.center.z) ||
            !std::isfinite(pixel.size) || pixel.size <= 0.0) {
            throw std::runtime_error(
                description + " pixel " + std::to_string(pixel.id) +
                " has invalid geometry");
        }
    }

    // CameraGeometry returns one pixel for a point. Overlapping interiors
    // therefore make the selected pixel depend on tie-breaking. Boundary
    // contact is allowed; only positive-area overlap is rejected.
    for (std::size_t i = 0; i < pixels.size(); ++i) {
        for (std::size_t j = i + 1; j < pixels.size(); ++j) {
            if (cameraPixelInteriorsOverlap(pixels[i], pixels[j])) {
                throw std::runtime_error(
                    description + " pixels " +
                    std::to_string(pixels[i].id) + " and " +
                    std::to_string(pixels[j].id) +
                    " have overlapping interiors");
            }
        }
    }
}

} // namespace

std::vector<std::pair<double, double>> readCollectorReflectivityCsv(const std::string& path)
{
    std::ifstream ifs(path);
    if (!ifs) {
        throw std::runtime_error("failed to open collector reflectivity CSV: " + path);
    }

    std::string line;
    if (!std::getline(ifs, line)) {
        throw std::runtime_error("collector reflectivity CSV is empty: " + path);
    }

    auto header_cells = splitCsvCells(line);
    std::map<std::string, int> header;
    for (int i = 0; i < static_cast<int>(header_cells.size()); ++i) {
        header[lowerCopy(header_cells[i])] = i;
    }

    std::vector<std::pair<double, double>> points;
    int line_no = 1;
    while (std::getline(ifs, line)) {
        ++line_no;
        if (trim(line).empty()) {
            continue;
        }
        auto cells = splitCsvCells(line);
        std::string theta = getOptionalCell(cells, header, "theta_deg");
        if (theta.empty()) theta = getOptionalCell(cells, header, "angle_deg");
        if (theta.empty()) theta = getOptionalCell(cells, header, "incidence_angle_deg");

        std::string reflectivity = getOptionalCell(cells, header, "reflectivity");
        if (reflectivity.empty()) reflectivity = getOptionalCell(cells, header, "efficiency");

        if (theta.empty() || reflectivity.empty()) {
            throw std::runtime_error(
                "collector reflectivity CSV line " + std::to_string(line_no) +
                " must provide theta_deg/angle_deg and reflectivity/efficiency");
        }
        double theta_deg = 0.0;
        double value = 0.0;
        if (!parseDoubleStrict(theta, theta_deg) ||
            !parseDoubleStrict(reflectivity, value) ||
            !std::isfinite(theta_deg) || theta_deg < 0.0 ||
            theta_deg > 90.0 ||
            !std::isfinite(value) || value < 0.0 || value > 1.0) {
            throw std::runtime_error(
                "collector reflectivity CSV line " +
                std::to_string(line_no) +
                " requires theta_deg in [0,90] and reflectivity in [0,1]");
        }
        points.emplace_back(theta_deg, value);
    }

    if (points.empty()) {
        throw std::runtime_error("collector reflectivity CSV has no data rows: " + path);
    }
    std::sort(points.begin(), points.end());
    for (std::size_t i = 1; i < points.size(); ++i) {
        if (points[i - 1].first == points[i].first) {
            throw std::runtime_error(
                "collector reflectivity CSV has duplicate angle " +
                doubleToString(points[i].first));
        }
    }
    return points;
}

CameraGeometry makeGeneratedCamera(const CameraConfig& cfg)
{
    if (cfg.radius_m <= 0.0 || cfg.pixel_size_m <= 0.0 || cfg.pixel_pitch_m <= 0.0) {
        throw std::runtime_error("camera radius, pixel_size_m, and pixel_pitch_m must be > 0");
    }

    CameraGeometry camera;
    const PixelShape shape = parsePixelShape(cfg.pixel_shape);
    const std::string mode = lowerCopy(trim(cfg.mode));
    int id = 0;

    if (mode == "square_grid" || mode == "grid") {
        int n = static_cast<int>(std::ceil(cfg.radius_m / cfg.pixel_pitch_m));
        for (int iy = -n; iy <= n; ++iy) {
            for (int ix = -n; ix <= n; ++ix) {
                double x = ix * cfg.pixel_pitch_m;
                double y = iy * cfg.pixel_pitch_m;
                if (x * x + y * y > cfg.radius_m * cfg.radius_m + 1e-14) {
                    continue;
                }
                CameraPixel p;
                p.id = id++;
                p.center = {x, y, 0.0};
                p.shape = shape;
                p.size = cfg.pixel_size_m;
                camera.addPixel(p);
            }
        }
        return camera;
    }

    if (mode == "hex_grid" || mode == "generated") {
        int n = static_cast<int>(std::ceil(cfg.radius_m / cfg.pixel_pitch_m)) + 2;
        const double y_step = cfg.pixel_pitch_m * std::sqrt(3.0) * 0.5;
        for (int r = -n; r <= n; ++r) {
            for (int q = -n; q <= n; ++q) {
                double x = cfg.pixel_pitch_m * (static_cast<double>(q) + 0.5 * r);
                double y = y_step * r;
                if (x * x + y * y > cfg.radius_m * cfg.radius_m + 1e-14) {
                    continue;
                }
                CameraPixel p;
                p.id = id++;
                p.center = {x, y, 0.0};
                p.shape = shape;
                p.size = cfg.pixel_size_m;
                camera.addPixel(p);
            }
        }
        return camera;
    }

    throw std::runtime_error("unsupported camera.mode: " + cfg.mode);
}

std::map<double, std::vector<SeriesFacetState>>
readElevationSeriesCsv(const std::string& path)
{
    std::ifstream ifs(path);
    if (!ifs) {
        throw std::runtime_error("failed to open elevation-series CSV: " + path);
    }

    std::string line;
    if (!std::getline(ifs, line)) {
        throw std::runtime_error("elevation-series CSV is empty: " + path);
    }

    auto header_cells = splitCsvCells(line);
    std::map<std::string, int> header;
    for (int i = 0; i < static_cast<int>(header_cells.size()); ++i) {
        header[lowerCopy(header_cells[i])] = i;
    }

    std::map<double, std::vector<SeriesFacetState>> by_angle;
    while (std::getline(ifs, line)) {
        line = trim(line);
        if (line.empty()) {
            continue;
        }
        auto cells = splitCsvCells(line);

        SeriesFacetState state;
        state.id = std::stoi(getRequiredCell(cells, header, "id"));
        const double elevation_deg = std::stod(getRequiredCell(cells, header, "elevation_deg"));
        state.center = parseVec3(
            getRequiredCell(cells, header, "center_x") + "," +
            getRequiredCell(cells, header, "center_y") + "," +
            getRequiredCell(cells, header, "center_z"),
            "series center");
        state.normal = parseVec3(
            getRequiredCell(cells, header, "normal_x") + "," +
            getRequiredCell(cells, header, "normal_y") + "," +
            getRequiredCell(cells, header, "normal_z"),
            "series normal").normalized();

        const std::string cell_index = getOptionalCell(cells, header, "cell_index");
        const std::string ring_index = getOptionalCell(cells, header, "ring_index");
        if (!cell_index.empty()) state.cell_index = std::stoi(cell_index);
        if (!ring_index.empty()) state.ring_index = std::stoi(ring_index);

        const std::string surface_type = getOptionalCell(cells, header, "surface_type");
        if (!surface_type.empty()) {
            state.has_surface_type = true;
            state.surface_type = parseSurfaceType(surface_type);
        }
        const std::string radius = getOptionalCell(cells, header, "radius_of_curvature");
        if (!radius.empty()) {
            double value = 0.0;
            if (!parseDoubleStrict(radius, value) || !std::isfinite(value) || value < 0.0) {
                throw std::runtime_error(
                    "radius_of_curvature must be finite and >= 0");
            }
            if (value > 0.0) {
                state.has_radius_of_curvature = true;
                state.radius_of_curvature = value;
            }
        }
        const std::string aperture_shape = getOptionalCell(cells, header, "aperture_shape");
        if (!aperture_shape.empty()) {
            state.has_aperture_shape = true;
            state.aperture_shape = parseApertureShape(aperture_shape);
        }
        const std::string size1 = getOptionalCell(cells, header, "size1");
        if (!size1.empty()) {
            state.has_size1 = true;
            state.size1 = std::stod(size1);
        }
        const std::string size2 = getOptionalCell(cells, header, "size2");
        if (!size2.empty()) {
            state.has_size2 = true;
            state.size2 = std::stod(size2);
        }
        const std::string rotation = getOptionalCell(cells, header, "aperture_rotation_rad");
        if (!rotation.empty()) {
            state.has_aperture_rotation = true;
            state.aperture_rotation_rad = std::stod(rotation);
        }

        auto read_non_negative = [&](const std::string& key,
                                     bool& has_value,
                                     double& value) {
            const std::string text = getOptionalCell(cells, header, key);
            if (text.empty()) {
                return;
            }
            if (!parseDoubleStrict(text, value) || !std::isfinite(value) || value < 0.0) {
                throw std::runtime_error(key + " must be finite and >= 0");
            }
            has_value = true;
        };
        read_non_negative("reflectivity_scale",
                          state.has_reflectivity_scale,
                          state.reflectivity_scale);
        read_non_negative("roughness_sigma_rad",
                          state.has_roughness_sigma_rad,
                          state.roughness_sigma_rad);
        read_non_negative("misalign_sigma_rad",
                          state.has_misalign_sigma_rad,
                          state.misalign_sigma_rad);

        by_angle[elevation_deg].push_back(state);
    }

    if (by_angle.empty()) {
        throw std::runtime_error("elevation-series CSV has no data rows: " + path);
    }

    for (auto& kv : by_angle) {
        auto& states = kv.second;
        std::sort(states.begin(), states.end(), [](const SeriesFacetState& a, const SeriesFacetState& b) {
            return a.id < b.id;
        });
    }
    return by_angle;
}

Vec3 slerpUnitVectors(const Vec3& a_in, const Vec3& b_in, double t) {
    Vec3 a = a_in.normalized();
    Vec3 b = b_in.normalized();
    double dot = std::clamp(a.dot(b), -1.0, 1.0);
    if (dot > 0.999999) {
        return (a * (1.0 - t) + b * t).normalized();
    }
    if (dot < -0.999999) {
        return (a * (1.0 - t) - b * t).normalized();
    }
    double theta = std::acos(dot);
    double s = std::sin(theta);
    return (a * (std::sin((1.0 - t) * theta) / s) +
            b * (std::sin(t * theta) / s)).normalized();
}

std::vector<MirrorFacet> buildElevationSeriesFacets(const std::map<std::string, std::string>& cfg) {
    const std::string format = lowerCopy(getString(cfg, "mirror.series_format", "raw1229"));
    if (format != "raw1229" && format != "facet_csv" && format != "csv") {
        throw std::runtime_error("unsupported mirror.series_format: " + format);
    }
    const double elevation_deg =
        getDouble(cfg, "mirror.series_elevation_deg",
                  getDouble(cfg, "telescope.pointing_el_deg", 90.0));

    SurfaceType surface_type =
        parseSurfaceType(getString(cfg, "mirror.surface_type", "Spherical"));
    ApertureShape aperture_shape =
        parseApertureShape(getString(cfg, "mirror.aperture_shape", "Hexagon"));
    const double radius_of_curvature =
        getDouble(cfg, "mirror.radius_of_curvature", 16.0);
    const double size1 = getDouble(cfg, "mirror.size1", 0.8);
    const double size2 = getDouble(cfg, "mirror.size2", 0.0);
    const double aperture_rotation_rad =
        getDouble(cfg, "mirror.aperture_rotation_rad", 0.0);

    if (format == "raw1229") {
        const std::string pattern = getString(cfg, "mirror.series_csv_pattern", "");
        if (pattern.empty()) {
            throw std::runtime_error("mirror.series_csv_pattern is required when mirror.series_format=raw1229");
        }

        const std::vector<double> angles_deg =
            parseDoubleList(getString(cfg, "mirror.series_angles_deg", ""), "mirror.series_angles_deg");
        const bool swap_xy = getBool(cfg, "mirror.series_swap_xy", true);

        std::vector<double> sorted_angles = angles_deg;
        std::sort(sorted_angles.begin(), sorted_angles.end());

        double lower_angle = sorted_angles.front();
        double upper_angle = sorted_angles.back();
        for (double angle : sorted_angles) {
            if (angle <= elevation_deg) lower_angle = angle;
            if (angle >= elevation_deg) {
                upper_angle = angle;
                break;
            }
        }
        if (elevation_deg <= sorted_angles.front()) {
            lower_angle = upper_angle = sorted_angles.front();
        } else if (elevation_deg >= sorted_angles.back()) {
            lower_angle = upper_angle = sorted_angles.back();
        }

        const auto lower_states = readRaw1229FacetStates(expandAnglePattern(pattern, lower_angle), swap_xy);
        const auto upper_states = (upper_angle == lower_angle)
            ? lower_states
            : readRaw1229FacetStates(expandAnglePattern(pattern, upper_angle), swap_xy);

        if (lower_states.size() != upper_states.size()) {
            throw std::runtime_error("mirror elevation series files have mismatched facet counts");
        }

        const double t = (upper_angle == lower_angle)
            ? 0.0
            : (elevation_deg - lower_angle) / (upper_angle - lower_angle);

        std::vector<MirrorFacet> facets;
        facets.reserve(lower_states.size());
        for (std::size_t i = 0; i < lower_states.size(); ++i) {
            const auto& lo = lower_states[i];
            const auto& hi = upper_states[i];
            if (lo.cell_index != hi.cell_index || lo.ring_index != hi.ring_index) {
                throw std::runtime_error("mirror elevation series files have inconsistent facet ordering");
            }

            MirrorFacet facet;
            facet.id = static_cast<int>(i);
            facet.center = lo.center * (1.0 - t) + hi.center * t;
            facet.normal = slerpUnitVectors(lo.normal, hi.normal, t);
            facet.surface_type = surface_type;
            facet.radius_of_curvature = radius_of_curvature;
            facet.aperture_shape = aperture_shape;
            facet.size1 = size1;
            facet.size2 = size2;
            facet.aperture_rotation_rad = aperture_rotation_rad;
            facet.reflectivity_scale = 1.0;
            facets.push_back(facet);
        }
        return facets;
    }

    const std::string csv_path = getString(cfg, "mirror.series_csv_path", "");
    if (csv_path.empty()) {
        throw std::runtime_error("mirror.series_csv_path is required when mirror.series_format=facet_csv");
    }
    auto by_angle = readElevationSeriesCsv(csv_path);
    std::vector<double> sorted_angles;
    for (const auto& kv : by_angle) sorted_angles.push_back(kv.first);

    double lower_angle = sorted_angles.front();
    double upper_angle = sorted_angles.back();
    for (double angle : sorted_angles) {
        if (angle <= elevation_deg) lower_angle = angle;
        if (angle >= elevation_deg) {
            upper_angle = angle;
            break;
        }
    }
    if (elevation_deg <= sorted_angles.front()) {
        lower_angle = upper_angle = sorted_angles.front();
    } else if (elevation_deg >= sorted_angles.back()) {
        lower_angle = upper_angle = sorted_angles.back();
    }

    const auto& lower_states = by_angle.at(lower_angle);
    const auto& upper_states = by_angle.at(upper_angle);
    if (lower_states.size() != upper_states.size()) {
        throw std::runtime_error("mirror series CSV has mismatched facet counts between elevations");
    }
    const double t = (upper_angle == lower_angle)
        ? 0.0
        : (elevation_deg - lower_angle) / (upper_angle - lower_angle);

    std::vector<MirrorFacet> facets;
    facets.reserve(lower_states.size());
    for (std::size_t i = 0; i < lower_states.size(); ++i) {
        const auto& lo = lower_states[i];
        const auto& hi = upper_states[i];
        if (lo.id != hi.id) {
            throw std::runtime_error("mirror series CSV has inconsistent facet ordering");
        }

        MirrorFacet facet;
        facet.id = lo.id;
        facet.center = lo.center * (1.0 - t) + hi.center * t;
        facet.normal = slerpUnitVectors(lo.normal, hi.normal, t);
        facet.surface_type = lo.has_surface_type ? lo.surface_type : surface_type;
        if (lo.has_radius_of_curvature != hi.has_radius_of_curvature) {
            throw std::runtime_error(
                "mirror series CSV has inconsistent radius_of_curvature for facet " +
                std::to_string(lo.id));
        }
        facet.radius_of_curvature = lo.has_radius_of_curvature
            ? lo.radius_of_curvature * (1.0 - t) + hi.radius_of_curvature * t
            : radius_of_curvature;
        facet.aperture_shape = lo.has_aperture_shape ? lo.aperture_shape : aperture_shape;
        facet.size1 = lo.has_size1 ? lo.size1 : size1;
        facet.size2 = lo.has_size2 ? lo.size2 : size2;
        facet.aperture_rotation_rad =
            lo.has_aperture_rotation ? lo.aperture_rotation_rad : aperture_rotation_rad;
        if (lo.has_reflectivity_scale != hi.has_reflectivity_scale ||
            lo.has_roughness_sigma_rad != hi.has_roughness_sigma_rad ||
            lo.has_misalign_sigma_rad != hi.has_misalign_sigma_rad) {
            throw std::runtime_error(
                "mirror series CSV has inconsistent optional per-facet fields for facet " +
                std::to_string(lo.id));
        }
        if (lo.has_reflectivity_scale) {
            facet.reflectivity_scale =
                lo.reflectivity_scale * (1.0 - t) + hi.reflectivity_scale * t;
            facet.has_reflectivity_scale = true;
        }
        if (lo.has_roughness_sigma_rad) {
            facet.roughness_sigma_rad =
                lo.roughness_sigma_rad * (1.0 - t) + hi.roughness_sigma_rad * t;
            facet.has_roughness_sigma_rad = true;
        }
        if (lo.has_misalign_sigma_rad) {
            facet.misalign_sigma_rad =
                lo.misalign_sigma_rad * (1.0 - t) + hi.misalign_sigma_rad * t;
            facet.has_misalign_sigma_rad = true;
        }
        facets.push_back(facet);
    }
    return facets;
}

std::map<std::string, std::string> loadScopedMirrorConfig(const std::string& path) {
    auto cfg = readKeyValueConfig(path);
    std::map<std::string, std::string> scoped;
    for (const auto& kv : cfg) {
        scoped[scopedComponentKey(kv.first, "mirror.")] = kv.second;
    }
    return scoped;
}

Vec3 perturbVectorOnSphere(const Vec3& direction,
                           double sigma_rad,
                           std::mt19937_64& rng)
{
    if (sigma_rad <= 0.0) {
        return direction.normalized();
    }
    Vec3 w = direction.normalized();
    Vec3 ref = (std::abs(w.z) < 0.9) ? Vec3{0.0, 0.0, 1.0}
                                     : Vec3{0.0, 1.0, 0.0};
    Vec3 u = ref.cross(w).normalized();
    Vec3 v = w.cross(u).normalized();
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double theta = sigma_rad * std::sqrt(-2.0 * std::log(std::max(1e-16, 1.0 - uniform(rng))));
    theta = std::min(theta, 3.14159265358979323846);
    double phi = 2.0 * 3.14159265358979323846 * uniform(rng);
    return (w * std::cos(theta) +
            u * (std::sin(theta) * std::cos(phi)) +
            v * (std::sin(theta) * std::sin(phi))).normalized();
}

std::vector<MirrorFacet> buildFacetsFromConfig(const std::map<std::string, std::string>& cfg) {
    std::string mode = lowerCopy(getString(cfg, "mirror.mode", "generated"));
    std::vector<MirrorFacet> facets;

    if (mode == "generated" || mode == "auto") {
        DishPrescription dish;
        dish.type = parseDishType(getString(cfg, "dish.type", "DaviesCotton"));
        dish.telescope_focal_length = getDouble(cfg, "dish.telescope_focal_length", 5.0);
        dish.dish_shape_length = getDouble(cfg, "dish.dish_shape_length",
                                           dish.telescope_focal_length);
        dish.dish_radius = getDouble(cfg, "dish.dish_radius", 2.0);
        dish.vertex = getVec3(cfg, "dish.vertex", {0.0, 0.0, 0.0});
        dish.optical_axis = getVec3(cfg, "dish.optical_axis", {0.0, 0.0, 1.0});

        FacetGridConfig grid;
        grid.facet_spacing = getDouble(cfg, "facet.spacing", 0.55);
        grid.facet_radius = getDouble(cfg, "facet.radius", 0.22);
        grid.aperture_shape = parseApertureShape(getString(cfg, "facet.aperture_shape", "Circular"));
        grid.surface_type = parseSurfaceType(getString(cfg, "facet.surface_type", "Spherical"));

        facets = FacetFactory::buildFacets(dish, grid);
    } else if (mode == "csv" || mode == "imported") {
        std::string path = getString(cfg, "mirror.csv_path", "");
        if (path.empty()) {
            throw std::runtime_error("mirror.csv_path is required when mirror.mode=csv");
        }
        MirrorFacetCsvDefaults defaults;
        defaults.surface_type =
            parseSurfaceType(getString(cfg, "mirror.surface_type", "Spherical"));
        defaults.radius_of_curvature =
            getDouble(cfg, "mirror.radius_of_curvature", 16.0);
        defaults.aperture_shape =
            parseApertureShape(getString(cfg, "mirror.aperture_shape", "Circular"));
        defaults.size1 = getDouble(cfg, "mirror.size1", 0.25);
        defaults.size2 = getDouble(cfg, "mirror.size2", 0.0);
        defaults.aperture_rotation_rad =
            getDouble(cfg, "mirror.aperture_rotation_rad", 0.0);
        std::string error;
        if (!readMirrorFacetsCSV(path, facets, &error, &defaults)) {
            throw std::runtime_error("failed to read mirror CSV: " + error);
        }
    } else if (mode == "elevation_series" || mode == "series") {
        facets = buildElevationSeriesFacets(cfg);
    } else {
        throw std::runtime_error("unsupported mirror.mode: " + mode);
    }

    std::string error;
    if (!validateMirrorFacets(facets, &error)) {
        throw std::runtime_error("invalid mirror facets: " + error);
    }
    return facets;
}

SyntheticPhotonConfig buildSourceConfig(const std::map<std::string, std::string>& cfg) {
    SyntheticPhotonConfig source;
    const std::string mode_text = getString(cfg, "source.mode", "ParallelBeam");
    if (!isPhotonCsvMode(mode_text)) {
        if (!isEventIOMode(mode_text)) {
            source.mode = parseSyntheticMode(mode_text);
        }
    }
    source.n_bunches = getInt(cfg, "source.n_bunches", source.n_bunches);
    source.multiplicity = getDouble(cfg, "source.multiplicity", source.multiplicity);
    source.wavelength_nm = getDouble(cfg, "source.wavelength_nm", source.wavelength_nm);
    source.time_ns = getDouble(cfg, "source.time_ns", source.time_ns);
    source.photon_weight = getDouble(cfg, "source.photon_weight", source.photon_weight);

    source.source_plane_z = getDouble(cfg, "source.source_plane_z", source.source_plane_z);
    source.beam_radius_m = getDouble(cfg, "source.beam_radius_m", source.beam_radius_m);
    source.beam_direction = getVec3(cfg, "source.beam_direction", source.beam_direction);

    source.source_position = getVec3(cfg, "source.source_position", source.source_position);
    source.aperture_z = getDouble(cfg, "source.aperture_z", source.aperture_z);
    source.aperture_radius_m = getDouble(cfg, "source.aperture_radius_m",
                                         source.aperture_radius_m);

    auto theta_it = cfg.find("source.beam_theta_deg");
    auto phi_it = cfg.find("source.beam_phi_deg");
    if (theta_it != cfg.end() || phi_it != cfg.end()) {
        double theta_deg = getDouble(cfg, "source.beam_theta_deg", 0.0);
        double phi_deg = getDouble(cfg, "source.beam_phi_deg", 0.0);
        double theta = theta_deg * DEG_TO_RAD;
        double phi = phi_deg * DEG_TO_RAD;
        source.beam_direction = {
            std::sin(theta) * std::cos(phi),
            std::sin(theta) * std::sin(phi),
            -std::cos(theta)
        };
    }

    source.event_id = getInt(cfg, "source.event_id", source.event_id);
    source.telescope_id = getInt(cfg, "source.telescope_id", source.telescope_id);
    source.random_seed = getUInt64(cfg, "source.random_seed", source.random_seed);
    return source;
}

SourceRuntimeConfig buildSourceRuntimeConfig(const std::map<std::string, std::string>& cfg) {
    SourceRuntimeConfig runtime;
    const std::string mode_text = getString(cfg, "source.mode", "ParallelBeam");
    runtime.use_photon_csv = isPhotonCsvMode(mode_text);
    runtime.use_eventio = isEventIOMode(mode_text);
    runtime.csv_path = getString(cfg, "source.csv_path", "");
    runtime.eventio_path = getString(cfg, "source.eventio_path", "");
    runtime.event_id_mode = getString(cfg, "source.event_id_mode", "event");
    const auto coordinate_frame_it = cfg.find("source.coordinate_frame");
    if (coordinate_frame_it != cfg.end()) {
        runtime.coordinate_frame = normalizeSourceCoordinateFrame(coordinate_frame_it->second);
    } else if (runtime.use_eventio) {
        runtime.coordinate_frame = normalizeSourceCoordinateFrame(
            getString(cfg, "source.eventio_coordinate_frame", "corsika_iact"));
    } else if (runtime.use_photon_csv &&
               cfg.find("source.local_telescope_frame") != cfg.end()) {
        runtime.coordinate_frame = getBool(cfg, "source.local_telescope_frame", true)
            ? "telescope_local"
            : "lact_generic_global";
    } else {
        runtime.coordinate_frame = "telescope_local";
    }
    runtime.eventio_coordinate_frame = runtime.coordinate_frame;
    const auto reference_z = cfg.find("source.eventio_reference_z_m");
    const auto legacy_2d_z = cfg.find("source.eventio_2d_input_plane_z_m");
    if (reference_z != cfg.end() && legacy_2d_z != cfg.end()) {
        const double new_value = getDouble(
            cfg, "source.eventio_reference_z_m", -16.0);
        const double legacy_value = getDouble(
            cfg, "source.eventio_2d_input_plane_z_m", -16.0);
        if (std::abs(new_value - legacy_value) > 1e-12) {
            throw std::runtime_error(
                "source.eventio_reference_z_m conflicts with deprecated "
                "source.eventio_2d_input_plane_z_m");
        }
    }
    runtime.eventio_reference_z_m = getDouble(
        cfg,
        "source.eventio_reference_z_m",
        getDouble(cfg, "source.eventio_2d_input_plane_z_m", -16.0));
    if (!std::isfinite(runtime.eventio_reference_z_m)) {
        throw std::runtime_error(
            "source.eventio_reference_z_m must be finite");
    }
    runtime.eventio_2d_plane_mode =
        lowerCopy(getString(cfg, "source.eventio_2d_plane_mode", "auto"));
    if (runtime.eventio_2d_plane_mode != "auto" &&
        runtime.eventio_2d_plane_mode != "forward" &&
        runtime.eventio_2d_plane_mode != "backproject") {
        throw std::runtime_error(
            "source.eventio_2d_plane_mode must be auto, forward, or backproject");
    }
    runtime.use_eventio_telescope_position =
        getBool(cfg, "source.use_eventio_telescope_position", true);
    runtime.csv_local_telescope_frame = runtime.coordinate_frame == "telescope_local";
    const std::string filter_value = getString(cfg, "source.filter_telescope_id", "");
    if (!trim(filter_value).empty() && !isDisabledText(filter_value)) {
        runtime.filter_telescope_id = true;
        runtime.selected_telescope_id = std::stoi(filter_value);
    }
    const std::string event_filter_value = getString(cfg, "source.filter_event_id", "");
    if (!trim(event_filter_value).empty() && !isDisabledText(event_filter_value)) {
        runtime.filter_event_id = true;
        runtime.selected_event_id = std::stoi(event_filter_value);
    }
    const std::string event_filters_value =
        getString(cfg, "source.filter_event_ids", "");
    if (!trim(event_filters_value).empty() &&
        !isDisabledText(event_filters_value)) {
        runtime.selected_event_ids =
            parseIntList(event_filters_value, "source.filter_event_ids");
    }
    if (runtime.filter_event_id && !runtime.selected_event_ids.empty()) {
        throw std::runtime_error(
            "source.filter_event_id and source.filter_event_ids are mutually exclusive");
    }
    const std::string shower_filter_value = getString(cfg, "source.filter_shower_event_id", "");
    if (!trim(shower_filter_value).empty() && !isDisabledText(shower_filter_value)) {
        runtime.filter_shower_event_id = true;
        runtime.selected_shower_event_id = std::stoi(shower_filter_value);
    }
    runtime.max_shower_events =
        getInt(cfg, "source.max_shower_events",
               getInt(cfg, "source.max_events", runtime.max_shower_events));
    return runtime;
}

#ifdef LACT_HAS_HESSIO
EventIOPhotonConfig buildEventIOPhotonConfig(const std::map<std::string, std::string>& cfg,
                                             const SyntheticPhotonConfig& source_cfg,
                                             const SourceRuntimeConfig& runtime_cfg) {
    EventIOPhotonConfig eventio;
    eventio.path = runtime_cfg.eventio_path;
    eventio.local_telescope_frame = runtime_cfg.csv_local_telescope_frame;
    eventio.event_id_mode = runtime_cfg.event_id_mode;
    eventio.default_wavelength_nm = source_cfg.wavelength_nm;
    eventio.missing_wavelength_model =
        getString(cfg, "source.missing_wavelength_model",
                  getString(cfg, "source.eventio_missing_wavelength_model",
                            eventio.missing_wavelength_model));
    eventio.missing_wavelength_min_nm =
        getDouble(cfg, "source.missing_wavelength_min_nm",
                  getDouble(cfg, "source.wavelength_min_nm",
                            eventio.missing_wavelength_min_nm));
    eventio.missing_wavelength_max_nm =
        getDouble(cfg, "source.missing_wavelength_max_nm",
                  getDouble(cfg, "source.wavelength_max_nm",
                            eventio.missing_wavelength_max_nm));
    eventio.missing_wavelength_seed =
        getUInt64(cfg, "source.missing_wavelength_seed",
                  eventio.missing_wavelength_seed);
    eventio.default_time_ns = source_cfg.time_ns;
    eventio.default_weight = source_cfg.photon_weight;
    eventio.default_multiplicity = source_cfg.multiplicity;
    eventio.default_event_id = source_cfg.event_id;
    eventio.default_telescope_id = source_cfg.telescope_id;
    eventio.filter_telescope_id = runtime_cfg.filter_telescope_id;
    eventio.selected_telescope_id = runtime_cfg.selected_telescope_id;
    eventio.filter_event_id = runtime_cfg.filter_event_id;
    eventio.selected_event_id = runtime_cfg.selected_event_id;
    eventio.selected_event_ids = runtime_cfg.selected_event_ids;
    eventio.filter_shower_event_id = runtime_cfg.filter_shower_event_id;
    eventio.selected_shower_event_id = runtime_cfg.selected_shower_event_id;
    eventio.max_shower_events = runtime_cfg.max_shower_events;
    eventio.read_emitter_info =
        getBool(cfg, "source.read_emitter_info",
                getBool(cfg, "source.eventio_read_emitter_info",
                        eventio.read_emitter_info));
    return eventio;
}

void printEventIOMetadata(const EventIOMetadata& metadata,
                          const SourceRuntimeConfig& runtime_cfg) {
    printSection("EventIO metadata");
    printField("input_lines", intToString(metadata.input_lines.size()));
    for (std::size_t i = 0; i < metadata.input_lines.size() && i < 12; ++i) {
        printField("input[" + intToString(i) + "]", metadata.input_lines[i]);
    }
    if (metadata.input_lines.size() > 12) {
        printField("input_more", intToString(metadata.input_lines.size() - 12));
    }

    printField("telescopes", intToString(metadata.telescopes.size()));
    for (const auto& tel : metadata.telescopes) {
        std::ostringstream label;
        label << "tel[" << tel.telescope_id << "]";
        std::ostringstream value;
        value << "pos_m=(" << doubleToString(tel.x_m) << ", "
              << doubleToString(tel.y_m) << ", "
              << doubleToString(tel.z_m) << "), radius_m="
              << doubleToString(tel.radius_m);
        printField(label.str(), value.str());
    }

    if (metadata.selected_event) {
        const auto& event = *metadata.selected_event;
        printField("selected_shower_event",
                   intToString(static_cast<std::uint64_t>(event.shower_event_id)));
        printField("primary_type",
                   intToString(static_cast<std::uint64_t>(event.primary_type)));
        printField("energy_gev", doubleToString(event.energy_gev));
        printField("theta_deg", doubleToString(event.theta_deg));
        printField("phi_deg", doubleToString(event.phi_deg));
        printField("az_north_to_east_deg",
                   doubleToString(event.azimuth_north_to_east_deg));
        printField("core_position_m",
                   "(" + doubleToString(event.core_x_m) + ", " +
                       doubleToString(event.core_y_m) + ")");
        printField("array_rotation_deg", doubleToString(event.array_rotation_deg));
    } else {
        printField("selected_event", "not found");
    }

    if (metadata.selected_event_offsets) {
        const auto& offsets = *metadata.selected_event_offsets;
        printField("selected_array_id",
                   intToString(static_cast<std::uint64_t>(metadata.selected_array_id)));
        printField("array_offsets",
                   intToString(static_cast<std::uint64_t>(offsets.x_m.size())));
        printField("array_time_offset_ns", doubleToString(offsets.time_offset_ns));
        for (std::size_t i = 0; i < offsets.x_m.size() && i < 20; ++i) {
            std::ostringstream label;
            label << "array[" << i << "]";
            std::ostringstream value;
            value << "offset_m=(" << doubleToString(offsets.x_m[i]) << ", "
                  << doubleToString(offsets.y_m[i]) << ")";
            if (i < offsets.weight.size()) {
                value << ", weight=" << doubleToString(offsets.weight[i]);
            }
            printField(label.str(), value.str());
        }
    } else {
        printField("array_offsets", "not found");
    }

    if (runtime_cfg.use_eventio_telescope_position) {
        printField("position_source", "EventIO telescope table");
    } else {
        printField("position_source", "telescope config");
    }
}
#endif


PhotonCsvConfig buildPhotonCsvConfig(const std::map<std::string, std::string>& cfg,
                                     const SyntheticPhotonConfig& source_cfg,
                                     const SourceRuntimeConfig& runtime_cfg) {
    PhotonCsvConfig csv;
    csv.csv_path = runtime_cfg.csv_path;
    csv.local_telescope_frame = runtime_cfg.csv_local_telescope_frame;
    csv.default_wavelength_nm = source_cfg.wavelength_nm;
    csv.default_time_ns = source_cfg.time_ns;
    csv.default_weight = source_cfg.photon_weight;
    csv.default_multiplicity = source_cfg.multiplicity;
    csv.default_event_id = source_cfg.event_id;
    csv.default_telescope_id = source_cfg.telescope_id;
    csv.default_eventio_2d =
        getBool(cfg, "source.eventio_2d", csv.default_eventio_2d);
    csv.filter_telescope_id = runtime_cfg.filter_telescope_id;
    csv.selected_telescope_id = runtime_cfg.selected_telescope_id;
    csv.filter_event_id = runtime_cfg.filter_event_id;
    csv.selected_event_id = runtime_cfg.selected_event_id;
    return csv;
}

OutputPlane buildOutputPlane(const std::map<std::string, std::string>& cfg) {
    OutputPlane plane;
    plane.point = getVec3(cfg, "output.plane_point", {0.0, 0.0, 5.0});
    plane.normal = getVec3(cfg, "output.plane_normal", {0.0, 0.0, 1.0});
    plane.buildLocalFrame();
    const auto u_it = cfg.find("output.plane_u_axis");
    const auto v_it = cfg.find("output.plane_v_axis");
    if (u_it != cfg.end() || v_it != cfg.end()) {
        if (u_it == cfg.end() || v_it == cfg.end()) {
            throw std::runtime_error("output.plane_u_axis and output.plane_v_axis must be set together");
        }
        Vec3 u = parseVec3(u_it->second, "output.plane_u_axis").normalized();
        Vec3 v = parseVec3(v_it->second, "output.plane_v_axis").normalized();
        const double uv_dot = std::abs(u.dot(v));
        const double un_dot = std::abs(u.dot(plane.normal));
        const double vn_dot = std::abs(v.dot(plane.normal));
        if (uv_dot > 1e-6 || un_dot > 1e-6 || vn_dot > 1e-6) {
            throw std::runtime_error(
                "output.plane_u_axis and output.plane_v_axis must be orthogonal to each other and to output.plane_normal");
        }
        const double handedness = u.cross(v).dot(plane.normal.normalized());
        if (handedness < 1.0 - 1e-6) {
            throw std::runtime_error(
                "output plane axes must satisfy u cross v = plane_normal; "
                "for the LACT focal plane normal (0,0,-1), use u=(-1,0,0), v=(0,1,0)");
        }
        plane.u_axis = u;
        plane.v_axis = v;
    }
    return plane;
}

CameraConfig buildCameraConfig(const std::map<std::string, std::string>& cfg) {
    CameraConfig camera;
    camera.enabled = getBool(cfg, "camera.enabled", false);
    camera.mode = getString(cfg, "camera.mode", camera.enabled ? "hex_grid" : "none");
    camera.csv_path = getString(cfg, "camera.csv_path", "");
    camera.pixel_shape = getString(cfg, "camera.pixel_shape", camera.pixel_shape);
    camera.pixel_size_m = getDouble(cfg, "camera.pixel_size_m", camera.pixel_size_m);
    camera.pixel_pitch_m = getDouble(cfg, "camera.pixel_pitch_m", camera.pixel_pitch_m);
    camera.radius_m = getDouble(cfg, "camera.radius_m", camera.radius_m);
    camera.collector = getString(cfg, "camera.collector", camera.collector);
    camera.collector_material =
        getString(cfg, "camera.collector_material", camera.collector_material);
    camera.collector_reflectivity_csv =
        getString(cfg, "camera.collector_reflectivity_csv", "");
    camera.collector_entrance_size_m =
        getDouble(cfg, "camera.collector_entrance_size_m",
                  camera.collector_entrance_size_m);
    camera.collector_exit_size_m =
        getDouble(cfg, "camera.collector_exit_size_m",
                  camera.collector_exit_size_m);
    camera.collector_height_m =
        getDouble(cfg, "camera.collector_height_m", camera.collector_height_m);
    if (isDisabledText(camera.mode)) {
        camera.enabled = false;
    }
    if (camera.enabled) {
        if (!std::isfinite(camera.pixel_size_m) || camera.pixel_size_m <= 0.0) {
            throw std::runtime_error("camera.pixel_size_m must be finite and > 0");
        }
        if (!std::isfinite(camera.pixel_pitch_m) || camera.pixel_pitch_m <= 0.0) {
            throw std::runtime_error("camera.pixel_pitch_m must be finite and > 0");
        }
        if (!std::isfinite(camera.radius_m) || camera.radius_m <= 0.0) {
            throw std::runtime_error("camera.radius_m must be finite and > 0");
        }
        if (!std::isfinite(camera.collector_entrance_size_m) ||
            camera.collector_entrance_size_m < 0.0) {
            throw std::runtime_error("camera.collector_entrance_size_m must be finite and >= 0");
        }
        if (!std::isfinite(camera.collector_exit_size_m) ||
            camera.collector_exit_size_m <= 0.0) {
            throw std::runtime_error("camera.collector_exit_size_m must be finite and > 0");
        }
        if (!std::isfinite(camera.collector_height_m) ||
            camera.collector_height_m <= 0.0) {
            throw std::runtime_error("camera.collector_height_m must be finite and > 0");
        }
    }
    return camera;
}

SipmConfig buildSipmConfig(const std::map<std::string, std::string>& cfg) {
    SipmConfig sipm;
    sipm.size_m = getDouble(cfg, "sipm.size_m", sipm.size_m);
    if (!std::isfinite(sipm.size_m) || sipm.size_m <= 0.0) {
        throw std::runtime_error("sipm.size_m must be finite and > 0");
    }
    return sipm;
}

ElectronicsResponse::ElectronicsResponse(const ElectronicsConfig& cfg)
    : cfg_(cfg)
{
    if (cfg_.channels_per_pixel <= 0) {
        throw std::runtime_error(
            "electronics.sipm.channels_per_pixel must be > 0");
    }
    if (cfg_.microcells_per_channel <= 0) {
        throw std::runtime_error(
            "electronics.sipm.microcells_per_channel must be > 0");
    }
    if (cfg_.saturation_model != "hard_no_recovery") {
        throw std::runtime_error(
            "electronics.sipm.saturation_model must be hard_no_recovery");
    }
}

double ElectronicsResponse::peConversion(double wavelength_nm) const
{
    (void)wavelength_nm;
    return 1.0;
}

double ElectronicsResponse::totalMicrocellsPerPixel() const
{
    return static_cast<double>(cfg_.channels_per_pixel) *
           static_cast<double>(cfg_.microcells_per_channel);
}

double ElectronicsResponse::saturatedPe(double primary_pe) const
{
    if (!std::isfinite(primary_pe) || primary_pe < 0.0) {
        throw std::runtime_error("integrated primary p.e. must be finite and >= 0");
    }
    if (!cfg_.saturation_enabled || primary_pe == 0.0) {
        return primary_pe;
    }
    const double cells = totalMicrocellsPerPixel();
    return -cells * std::expm1(-primary_pe / cells);
}

std::vector<double> ElectronicsResponse::saturatedWaveform(
    const std::vector<double>& primary_pe_by_time_bin) const
{
    for (const double primary_pe : primary_pe_by_time_bin) {
        if (!std::isfinite(primary_pe) || primary_pe < 0.0) {
            throw std::runtime_error(
                "primary p.e. waveform samples must be finite and >= 0");
        }
    }
    if (!cfg_.saturation_enabled) {
        return primary_pe_by_time_bin;
    }
    std::vector<double> fired_pe_by_time_bin;
    fired_pe_by_time_bin.reserve(primary_pe_by_time_bin.size());
    double cumulative_primary_pe = 0.0;
    double cumulative_fired_pe = 0.0;
    for (const double primary_pe : primary_pe_by_time_bin) {
        cumulative_primary_pe += primary_pe;
        const double next_cumulative_fired_pe =
            saturatedPe(cumulative_primary_pe);
        fired_pe_by_time_bin.push_back(std::max(
            0.0, next_cumulative_fired_pe - cumulative_fired_pe));
        cumulative_fired_pe = next_cumulative_fired_pe;
    }
    return fired_pe_by_time_bin;
}

ElectronicsConfig buildElectronicsConfig(const std::map<std::string, std::string>& cfg) {
    ElectronicsConfig electronics;
    electronics.channels_per_pixel = getInt(
        cfg, "electronics.sipm.channels_per_pixel",
        electronics.channels_per_pixel);
    electronics.microcells_per_channel = getInt(
        cfg, "electronics.sipm.microcells_per_channel",
        electronics.microcells_per_channel);
    electronics.saturation_enabled = getBool(
        cfg, "electronics.sipm.saturation_enabled",
        electronics.saturation_enabled);
    electronics.saturation_model = lowerCopy(trim(getString(
        cfg, "electronics.sipm.saturation_model",
        electronics.saturation_model)));
    // Validate the parsed values even when only the configuration object is
    // used later by an output writer.
    (void)ElectronicsResponse(electronics);
    return electronics;
}

std::vector<double> readNormalizedPositiveCsvColumn(
    const std::string& path,
    const std::string& requested_column)
{
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error(
            "failed to open empirical single-p.e. charge CSV: " + path);
    }
    std::string line;
    if (!std::getline(input, line)) {
        throw std::runtime_error(
            "empty empirical single-p.e. charge CSV: " + path);
    }
    const auto header = splitCsvCells(line);
    const auto column_it = std::find(
        header.begin(), header.end(), requested_column);
    if (column_it == header.end()) {
        throw std::runtime_error(
            "missing column '" + requested_column +
            "' in empirical single-p.e. charge CSV: " + path);
    }
    const std::size_t column = static_cast<std::size_t>(
        std::distance(header.begin(), column_it));
    std::vector<double> values;
    while (std::getline(input, line)) {
        if (trim(line).empty() || trim(line).front() == '#') continue;
        const auto cells = splitCsvCells(line);
        if (column >= cells.size()) continue;
        try {
            const double value = std::stod(cells[column]);
            if (std::isfinite(value) && value > 0.0) {
                values.push_back(value);
            }
        } catch (const std::exception&) {
            // Ignore non-numeric rows; the final empty check remains strict.
        }
    }
    if (values.empty()) {
        throw std::runtime_error(
            "no positive numeric values in empirical single-p.e. charge CSV: " +
            path);
    }
    double mean = 0.0;
    for (const double value : values) mean += value;
    mean /= static_cast<double>(values.size());
    for (double& value : values) value /= mean;
    return values;
}

void applyMicrocellDeviceConfig(
    electronics::MicrocellConfig& microcell,
    const std::string& path)
{
    if (trim(path).empty()) return;
    const auto device = readKeyValueConfig(path);
    const auto int_value = [&](const std::string& key, int fallback) {
        return getInt(
            device, key,
            getInt(device, "microcell." + key, fallback));
    };
    const auto double_value = [&](const std::string& key, double fallback) {
        return getDouble(
            device, key,
            getDouble(device, "microcell." + key, fallback));
    };
    const auto string_value = [&](const std::string& key,
                                  const std::string& fallback) {
        return getString(
            device, key,
            getString(device, "microcell." + key, fallback));
    };

    microcell.layout = lowerCopy(trim(string_value("layout", microcell.layout)));
    microcell.channel_columns = int_value(
        "channel_columns", microcell.channel_columns);
    microcell.channel_rows = int_value("channel_rows", microcell.channel_rows);
    microcell.channel_size_x_m = double_value(
        "channel_size_x_m", microcell.channel_size_x_m);
    microcell.channel_size_y_m = double_value(
        "channel_size_y_m", microcell.channel_size_y_m);
    microcell.channel_gap_x_m = double_value(
        "channel_gap_x_m", microcell.channel_gap_x_m);
    microcell.channel_gap_y_m = double_value(
        "channel_gap_y_m", microcell.channel_gap_y_m);
    microcell.microcell_columns_per_channel = int_value(
        "microcell_columns_per_channel",
        microcell.microcell_columns_per_channel);
    microcell.microcell_rows_per_channel = int_value(
        "microcell_rows_per_channel",
        microcell.microcell_rows_per_channel);

    // All redundant totals are derived from the independent device inputs.
    microcell.channels_per_pixel =
        microcell.channel_columns * microcell.channel_rows;
    microcell.microcells_per_channel =
        microcell.microcell_columns_per_channel *
        microcell.microcell_rows_per_channel;
    microcell.grid_columns =
        microcell.channel_columns *
        microcell.microcell_columns_per_channel;
    microcell.grid_rows =
        microcell.channel_rows *
        microcell.microcell_rows_per_channel;
    microcell.sensor_size_x_m =
        microcell.channel_columns * microcell.channel_size_x_m +
        (microcell.channel_columns - 1) * microcell.channel_gap_x_m;
    microcell.sensor_size_y_m =
        microcell.channel_rows * microcell.channel_size_y_m +
        (microcell.channel_rows - 1) * microcell.channel_gap_y_m;
}

electronics::DetectorPipelineConfig buildDetectorPipelineConfig(
    const std::map<std::string, std::string>& cfg)
{
    electronics::DetectorPipelineConfig out;
    if (cfg.find("electronics.pipeline.enabled") != cfg.end()) {
        throw std::runtime_error(
            "electronics.pipeline.enabled has been removed; "
            "use electronics.enabled instead");
    }
    out.enabled = getBool(cfg, "electronics.enabled", out.enabled);
    out.random_seed = getUInt64(
        cfg, "electronics.random_seed", out.random_seed);

    const std::string microcell_device = getString(
        cfg, "electronics.microcell.device", "");
    applyMicrocellDeviceConfig(out.microcell, microcell_device);
    out.microcell.enabled = getBool(
        cfg, "electronics.microcell.enabled", out.microcell.enabled);
    out.microcell.saturation_enabled = getBool(
        cfg, "electronics.microcell.saturation_enabled",
        out.microcell.saturation_enabled);
    out.microcell.model = lowerCopy(trim(getString(
        cfg, "electronics.microcell.model", out.microcell.model)));
    out.microcell.layout = lowerCopy(trim(getString(
        cfg, "electronics.microcell.layout", out.microcell.layout)));
    out.microcell.sensor_size_x_m = getDouble(
        cfg, "electronics.microcell.sensor_size_x_m",
        getDouble(cfg, "sipm.size_m", out.microcell.sensor_size_x_m));
    out.microcell.sensor_size_y_m = getDouble(
        cfg, "electronics.microcell.sensor_size_y_m",
        getDouble(cfg, "sipm.size_m", out.microcell.sensor_size_y_m));
    out.microcell.grid_columns = getInt(
        cfg, "electronics.microcell.grid_columns",
        out.microcell.grid_columns);
    out.microcell.grid_rows = getInt(
        cfg, "electronics.microcell.grid_rows", out.microcell.grid_rows);
    out.microcell.channels_per_pixel = getInt(
        cfg, "electronics.microcell.channels_per_pixel",
        out.microcell.channels_per_pixel);
    out.microcell.microcells_per_channel = getInt(
        cfg, "electronics.microcell.microcells_per_channel",
        out.microcell.microcells_per_channel);
    out.microcell.channel_columns = getInt(
        cfg, "electronics.microcell.channel_columns",
        out.microcell.channel_columns);
    out.microcell.channel_rows = getInt(
        cfg, "electronics.microcell.channel_rows",
        out.microcell.channel_rows);
    out.microcell.channel_size_x_m = getDouble(
        cfg, "electronics.microcell.channel_size_x_m",
        out.microcell.channel_size_x_m);
    out.microcell.channel_size_y_m = getDouble(
        cfg, "electronics.microcell.channel_size_y_m",
        out.microcell.channel_size_y_m);
    out.microcell.channel_gap_x_m = getDouble(
        cfg, "electronics.microcell.channel_gap_x_m",
        out.microcell.channel_gap_x_m);
    out.microcell.channel_gap_y_m = getDouble(
        cfg, "electronics.microcell.channel_gap_y_m",
        out.microcell.channel_gap_y_m);
    out.microcell.microcell_columns_per_channel = getInt(
        cfg, "electronics.microcell.microcell_columns_per_channel",
        out.microcell.microcell_columns_per_channel);
    out.microcell.microcell_rows_per_channel = getInt(
        cfg, "electronics.microcell.microcell_rows_per_channel",
        out.microcell.microcell_rows_per_channel);
    out.microcell.pde_includes_inter_channel_gaps = getBool(
        cfg, "electronics.microcell.pde_includes_inter_channel_gaps",
        out.microcell.pde_includes_inter_channel_gaps);
    if (getBool(cfg, "electronics.microcell.recovery_enabled", false)) {
        throw std::runtime_error(
            "microcell recovery is not available; use "
            "electronics.microcell.recovery_enabled=false");
    }
    out.single_pe.enabled = getBool(
        cfg, "electronics.single_pe.enabled", out.single_pe.enabled);
    out.single_pe.model = lowerCopy(trim(getString(
        cfg, "electronics.single_pe.model", out.single_pe.model)));
    out.single_pe.csv_path = getString(
        cfg, "electronics.single_pe.template",
        getString(cfg, "electronics.single_pe.csv",
            getString(cfg, "electronics.single_pe.csv_path",
                      out.single_pe.csv_path)));
    out.single_pe.template_time_reference = lowerCopy(trim(getString(
        cfg, "electronics.single_pe.template_time_reference",
        out.single_pe.template_time_reference)));
    out.single_pe.unit = lowerCopy(trim(getString(
        cfg, "electronics.single_pe.unit", out.single_pe.unit)));
    out.single_pe.rise_ns = getDouble(
        cfg, "electronics.single_pe.rise_ns", out.single_pe.rise_ns);
    out.single_pe.fall_ns = getDouble(
        cfg, "electronics.single_pe.fall_ns", out.single_pe.fall_ns);
    out.single_pe.support_ns = getDouble(
        cfg, "electronics.single_pe.support_ns", out.single_pe.support_ns);
    out.single_pe.amplitude_scale = getDouble(
        cfg, "electronics.single_pe.amplitude_scale",
        out.single_pe.amplitude_scale);

    auto& charge = out.single_pe.charge_fluctuation;
    charge.enabled = getBool(
        cfg, "electronics.single_pe.charge_fluctuation.enabled",
        charge.enabled);
    charge.model = lowerCopy(trim(getString(
        cfg, "electronics.single_pe.charge_fluctuation.model",
        charge.model)));
    charge.csv_path = getString(
        cfg, "electronics.single_pe.charge_fluctuation.samples",
        getString(cfg,
            "electronics.single_pe.charge_fluctuation.csv_path",
            charge.csv_path));
    charge.csv_column = getString(
        cfg, "electronics.single_pe.charge_fluctuation.csv_column",
        charge.csv_column);
    if (charge.enabled && !charge.csv_path.empty()) {
        charge.empirical_samples = readNormalizedPositiveCsvColumn(
            charge.csv_path, charge.csv_column);
    }

    auto& jitter = out.single_pe.time_jitter;
    jitter.enabled = getBool(
        cfg, "electronics.single_pe.time_jitter.enabled", jitter.enabled);
    jitter.model = lowerCopy(trim(getString(
        cfg, "electronics.single_pe.time_jitter.model", jitter.model)));
    jitter.sigma_ns = getDouble(
        cfg, "electronics.single_pe.time_jitter.sigma_ns",
        jitter.sigma_ns);

    out.sampling.width_ns = getDouble(
        cfg, "electronics.sampling.width_ns",
        getDouble(cfg, "waveform.time_bin_width_ns",
                  out.sampling.width_ns));
    out.sampling.start_ns = getDouble(
        cfg, "electronics.sampling.start_ns",
        getDouble(cfg, "waveform.time_window_start_ns",
                  out.sampling.start_ns));
    out.sampling.end_ns = getDouble(
        cfg, "electronics.sampling.end_ns",
        getDouble(cfg, "waveform.time_window_end_ns",
                  out.sampling.end_ns));

    out.camera_trigger.enabled = getBool(
        cfg, "trigger.camera.enabled",
        getBool(cfg, "trigger.enabled", out.camera_trigger.enabled));
    out.camera_trigger.mode = lowerCopy(trim(getString(
        cfg, "trigger.camera.mode", out.camera_trigger.mode)));
    out.camera_trigger.pixel_threshold_pe = getDouble(
        cfg, "trigger.camera.pixel_threshold_pe",
        getDouble(cfg, "trigger.pixel_threshold_pe",
                  out.camera_trigger.pixel_threshold_pe));
    out.camera_trigger.pixel_threshold_mv = getDouble(
        cfg, "trigger.camera.pixel_threshold_mv",
        out.camera_trigger.pixel_threshold_mv);
    out.camera_trigger.multiplicity = getInt(
        cfg, "trigger.camera.multiplicity",
        getInt(cfg, "trigger.camera_multiplicity",
               out.camera_trigger.multiplicity));
    out.camera_trigger.coincidence_window_ns = getDouble(
        cfg, "trigger.camera.coincidence_window_ns",
        getDouble(cfg, "trigger.camera_coincidence_window_ns",
                  out.camera_trigger.coincidence_window_ns));
    out.save_primary_sequence = getBool(
        cfg, "electronics.output.save_primary_sequence",
        out.save_primary_sequence);
    out.save_fired_sequence = getBool(
        cfg, "electronics.output.save_fired_sequence",
        out.save_fired_sequence);
    out.save_microcell_decisions = getBool(
        cfg, "electronics.output.save_microcell_decisions",
        out.save_microcell_decisions);
    out.save_channel_waveforms = getBool(
        cfg, "electronics.output.save_channel_waveforms",
        out.save_channel_waveforms);
    electronics::validateDetectorPipelineConfig(out);
    return out;
}

void validateCameraDetectorCompatibility(
    const CameraConfig& camera,
    const electronics::DetectorPipelineConfig& detector)
{
    if (!detector.enabled || !detector.microcell.enabled ||
        detector.microcell.layout != "s17351_tiled_2x4") {
        return;
    }
    if (!camera.enabled) {
        throw std::runtime_error(
            "explicit S17351 microcell geometry requires camera.enabled=true");
    }
    if (isDisabledText(camera.collector)) {
        throw std::runtime_error(
            "explicit S17351 microcell geometry requires an enabled light "
            "collector because microcell positions use collector-exit "
            "coordinates; refusing the previous silent (0,0) gap rejection");
    }
}

NsbConfig buildNsbConfig(const std::map<std::string, std::string>& cfg) {
    NsbConfig nsb;
    nsb.enabled = getBool(cfg, "nsb.enabled", nsb.enabled);
    nsb.model = lowerCopy(trim(getString(cfg, "nsb.model", nsb.model)));
    nsb.rate_pe_per_ns_per_pixel =
        getDouble(cfg, "nsb.rate_pe_per_ns_per_pixel", nsb.rate_pe_per_ns_per_pixel);
    nsb.window_ns = getDouble(cfg, "nsb.window_ns", nsb.window_ns);
    nsb.seed = getUInt64(cfg, "nsb.seed", nsb.seed);
    nsb.spectrum_csv = getString(cfg, "nsb.spectrum_csv", nsb.spectrum_csv);
    nsb.spectrum_unit =
        lowerCopy(trim(getString(cfg, "nsb.spectrum_unit", nsb.spectrum_unit)));
    nsb.effective_area_m2 =
        getDouble(cfg, "nsb.effective_area_m2", nsb.effective_area_m2);
    nsb.collector_mean_transmission = getDouble(
        cfg, "nsb.collector_mean_transmission",
        nsb.collector_mean_transmission);
    nsb.collector_mean_transmission_configured =
        cfg.find("nsb.collector_mean_transmission") != cfg.end();
    nsb.pixel_solid_angle =
        lowerCopy(trim(getString(cfg, "nsb.pixel_solid_angle", nsb.pixel_solid_angle)));
    nsb.pixel_solid_angle_sr =
        getDouble(cfg, "nsb.pixel_solid_angle_sr", nsb.pixel_solid_angle_sr);

    if (!(nsb.model == "constant_rate" || nsb.model == "spectral_flux" ||
          nsb.model == "none" || nsb.model == "off")) {
        throw std::runtime_error("nsb.model must be constant_rate, spectral_flux, or none");
    }
    if (!std::isfinite(nsb.rate_pe_per_ns_per_pixel) ||
        nsb.rate_pe_per_ns_per_pixel < 0.0) {
        throw std::runtime_error("nsb.rate_pe_per_ns_per_pixel must be finite and >= 0");
    }
    if (!std::isfinite(nsb.window_ns) || nsb.window_ns < 0.0) {
        throw std::runtime_error("nsb.window_ns must be finite and >= 0");
    }
    if (!std::isfinite(nsb.effective_area_m2) || nsb.effective_area_m2 < 0.0) {
        throw std::runtime_error("nsb.effective_area_m2 must be finite and >= 0");
    }
    if (!std::isfinite(nsb.collector_mean_transmission) ||
        nsb.collector_mean_transmission <= 0.0 ||
        nsb.collector_mean_transmission > 1.0) {
        throw std::runtime_error(
            "nsb.collector_mean_transmission must be finite and in (0, 1]");
    }
    if (!std::isfinite(nsb.pixel_solid_angle_sr) || nsb.pixel_solid_angle_sr < 0.0) {
        throw std::runtime_error("nsb.pixel_solid_angle_sr must be finite and >= 0");
    }
    if (!(nsb.pixel_solid_angle == "auto" || nsb.pixel_solid_angle == "manual" ||
          nsb.pixel_solid_angle == "fixed")) {
        throw std::runtime_error("nsb.pixel_solid_angle must be auto, manual, or fixed");
    }
    if (!(nsb.spectrum_unit == "ph_s_nm_sr_m2" ||
          nsb.spectrum_unit == "photons_m2_s_sr_nm")) {
        throw std::runtime_error(
            "nsb.spectrum_unit must be ph_s_nm_sr_m2 or photons_m2_s_sr_nm");
    }
    if (nsb.model == "spectral_flux") {
        if (trim(nsb.spectrum_csv).empty()) {
            throw std::runtime_error("nsb.spectrum_csv is required for nsb.model=spectral_flux");
        }
        if (nsb.effective_area_m2 <= 0.0) {
            throw std::runtime_error(
                "nsb.effective_area_m2 must be > 0 for nsb.model=spectral_flux");
        }
        if (!nsb.collector_mean_transmission_configured) {
            throw std::runtime_error(
                "nsb.collector_mean_transmission must be explicitly set for "
                "nsb.model=spectral_flux; use a measured/simulated collector "
                "average, or set 1 only when effective_area_m2 already "
                "includes the collector");
        }
    }
    if (isDisabledText(nsb.model)) {
        nsb.enabled = false;
    }
    return nsb;
}

std::vector<std::pair<double, double>> readNsbSpectrumCsv(const std::string& path)
{
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("failed to open NSB spectrum CSV: " + path);
    }
    std::vector<std::pair<double, double>> rows;
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') {
            continue;
        }
        const auto cells = splitCsvCells(line);
        if (cells.size() < 2) {
            continue;
        }
        double wavelength = 0.0;
        double flux = 0.0;
        if (!parseDoubleStrict(cells[0], wavelength) ||
            !parseDoubleStrict(cells[1], flux)) {
            continue;
        }
        if (!std::isfinite(wavelength) || !std::isfinite(flux) ||
            wavelength <= 0.0 || flux < 0.0) {
            throw std::runtime_error("invalid NSB spectrum row in " + path + ": " + line);
        }
        rows.push_back({wavelength, flux});
    }
    if (rows.size() < 2) {
        throw std::runtime_error("NSB spectrum needs at least two numeric rows: " + path);
    }
    std::sort(rows.begin(), rows.end());
    return rows;
}

void resolveNsbSpectralRate(NsbConfig& nsb,
                            const OpticalEfficiencyConfig& efficiency_cfg,
                            const CameraGeometry& camera,
                            const TelescopeConfig& telescope,
                            const electronics::DetectorPipelineConfig* detector)
{
    nsb.computed_from_spectrum = false;
    nsb.spectral_integral_pe_s_sr_m2 = 0.0;
    if (!nsb.enabled || nsb.model != "spectral_flux") {
        return;
    }

    if (nsb.pixel_solid_angle == "auto") {
        if (camera.size() == 0) {
            throw std::runtime_error(
                "nsb.pixel_solid_angle=auto requires a configured camera");
        }
        const double pixel_size_m = camera.pixels().front().size;
        if (!std::isfinite(pixel_size_m) || pixel_size_m <= 0.0 ||
            !std::isfinite(telescope.focal_length_m) || telescope.focal_length_m <= 0.0) {
            throw std::runtime_error(
                "cannot compute automatic NSB pixel solid angle from camera/focal length");
        }
        nsb.pixel_solid_angle_sr =
            (pixel_size_m / telescope.focal_length_m) *
            (pixel_size_m / telescope.focal_length_m);
    } else if (nsb.pixel_solid_angle_sr <= 0.0) {
        throw std::runtime_error(
            "nsb.pixel_solid_angle_sr must be > 0 when nsb.pixel_solid_angle is manual/fixed");
    }

    OpticalEfficiency efficiency(efficiency_cfg);
    const auto spectrum = readNsbSpectrumCsv(nsb.spectrum_csv);
    double integral = 0.0;
    for (std::size_t i = 1; i < spectrum.size(); ++i) {
        const double w0 = spectrum[i - 1].first;
        const double w1 = spectrum[i].first;
        const double f0 = spectrum[i - 1].second * efficiency.total(w0, 0.0);
        const double f1 = spectrum[i].second * efficiency.total(w1, 0.0);
        integral += 0.5 * (f0 + f1) * (w1 - w0);
    }
    nsb.spectral_integral_pe_s_sr_m2 = integral;
    nsb.microcell_geometric_acceptance = 1.0;
    if (detector && detector->enabled && detector->microcell.enabled &&
        detector->microcell.layout == "s17351_tiled_2x4" &&
        !detector->microcell.pde_includes_inter_channel_gaps) {
        nsb.microcell_geometric_acceptance =
            electronics::interChannelActiveFraction(detector->microcell);
    }
    nsb.rate_pe_per_ns_per_pixel =
        1.0e-9 * integral * nsb.effective_area_m2 *
        nsb.pixel_solid_angle_sr * nsb.collector_mean_transmission *
        nsb.microcell_geometric_acceptance;
    nsb.computed_from_spectrum = true;
}

void generateIntegratedNsbPe(const NsbConfig& nsb,
                             int event_id,
                             int telescope_id,
                             std::size_t n_pixels,
                             double window_ns,
                             const std::function<void(std::size_t, float)>& add_sample)
{
    if (!nsb.enabled || nsb.rate_pe_per_ns_per_pixel <= 0.0 ||
        window_ns <= 0.0 || n_pixels == 0) {
        return;
    }

    auto rng = makeNsbEventTelescopeRng(nsb, event_id, telescope_id);
    std::poisson_distribution<int> poisson(nsb.rate_pe_per_ns_per_pixel * window_ns);
    for (std::size_t col = 0; col < n_pixels; ++col) {
        const auto pe = static_cast<float>(poisson(rng));
        if (pe > 0.0f) {
            add_sample(col, pe);
        }
    }
}

void generateTimeBinnedNsbPe(const NsbConfig& nsb,
                             const WaveformOutputConfig& waveform_cfg,
                             int event_id,
                             int telescope_id,
                             std::size_t n_pixels,
                             std::size_t n_bins,
                             const std::function<void(std::size_t, std::size_t, float)>& add_sample)
{
    if (!nsb.enabled || nsb.rate_pe_per_ns_per_pixel <= 0.0 ||
        waveform_cfg.time_bin_width_ns <= 0.0 || n_pixels == 0 || n_bins == 0) {
        return;
    }

    auto rng = makeNsbEventTelescopeRng(nsb, event_id, telescope_id);
    const double total_mean =
        nsb.rate_pe_per_ns_per_pixel *
        waveform_cfg.time_bin_width_ns *
        static_cast<double>(n_pixels) *
        static_cast<double>(n_bins);
    std::poisson_distribution<unsigned long long> total_poisson(total_mean);
    const auto total_pe = static_cast<std::uint64_t>(total_poisson(rng));
    if (total_pe == 0) {
        return;
    }

    std::uniform_int_distribution<unsigned long long> pixel_dist(
        0ULL, static_cast<unsigned long long>(n_pixels - 1));
    std::uniform_int_distribution<unsigned long long> bin_dist(
        0ULL, static_cast<unsigned long long>(n_bins - 1));
    for (std::uint64_t i = 0; i < total_pe; ++i) {
        add_sample(static_cast<std::size_t>(pixel_dist(rng)),
                   static_cast<std::size_t>(bin_dist(rng)),
                   1.0f);
    }
}

TriggerConfig buildTriggerConfig(const std::map<std::string, std::string>& cfg) {
    TriggerConfig trigger;
    trigger.enabled = getBool(
        cfg, "trigger.camera.enabled",
        getBool(cfg, "trigger.enabled", trigger.enabled));
    trigger.array_enabled = getBool(
        cfg, "trigger.array.enabled", trigger.array_enabled);
    trigger.pixel_threshold_pe = getDouble(
        cfg, "trigger.camera.pixel_threshold_pe",
        getDouble(cfg, "trigger.pixel_threshold_pe",
                  trigger.pixel_threshold_pe));
    trigger.camera_multiplicity = getInt(
        cfg, "trigger.camera.multiplicity",
        getInt(cfg, "trigger.camera_multiplicity",
               trigger.camera_multiplicity));
    trigger.array_multiplicity = getInt(
        cfg, "trigger.array.multiplicity",
        getInt(cfg, "trigger.array_multiplicity",
               trigger.array_multiplicity));
    const bool has_legacy_coincidence_window =
        cfg.find("trigger.coincidence_window_ns") != cfg.end();
    trigger.coincidence_window_ns = getDouble(
        cfg, "trigger.coincidence_window_ns", trigger.coincidence_window_ns);
    trigger.camera_coincidence_window_ns = getDouble(
        cfg, "trigger.camera.coincidence_window_ns",
        getDouble(
            cfg,
            "trigger.camera_coincidence_window_ns",
            has_legacy_coincidence_window
                ? trigger.coincidence_window_ns
                : trigger.camera_coincidence_window_ns));
    trigger.array_coincidence_window_ns = getDouble(
        cfg,
        "trigger.array_coincidence_window_ns",
        has_legacy_coincidence_window
            ? trigger.coincidence_window_ns
            : trigger.array_coincidence_window_ns);
    trigger.array_time_correction = lowerCopy(trim(getString(
        cfg, "trigger.array_time_correction",
        trigger.array_time_correction)));
    trigger.array_wavefront_speed_m_per_ns = getDouble(
        cfg, "trigger.array_wavefront_speed_m_per_ns",
        trigger.array_wavefront_speed_m_per_ns);

    if (!std::isfinite(trigger.pixel_threshold_pe) || trigger.pixel_threshold_pe < 0.0) {
        throw std::runtime_error("trigger.pixel_threshold_pe must be finite and >= 0");
    }
    if (trigger.camera_multiplicity <= 0) {
        throw std::runtime_error("trigger.camera_multiplicity must be > 0");
    }
    if (trigger.array_multiplicity <= 0) {
        throw std::runtime_error("trigger.array_multiplicity must be > 0");
    }
    if (!std::isfinite(trigger.coincidence_window_ns) ||
        trigger.coincidence_window_ns < 0.0) {
        throw std::runtime_error("trigger.coincidence_window_ns must be finite and >= 0");
    }
    if (!std::isfinite(trigger.camera_coincidence_window_ns) ||
        trigger.camera_coincidence_window_ns < 0.0) {
        throw std::runtime_error(
            "trigger.camera_coincidence_window_ns must be finite and >= 0");
    }
    if (!std::isfinite(trigger.array_coincidence_window_ns) ||
        trigger.array_coincidence_window_ns < 0.0) {
        throw std::runtime_error(
            "trigger.array_coincidence_window_ns must be finite and >= 0");
    }
    if (trigger.array_time_correction != "none" &&
        trigger.array_time_correction != "plane_wave") {
        throw std::runtime_error(
            "trigger.array_time_correction must be none or plane_wave");
    }
    if (!std::isfinite(trigger.array_wavefront_speed_m_per_ns) ||
        trigger.array_wavefront_speed_m_per_ns < 0.0) {
        throw std::runtime_error(
            "trigger.array_wavefront_speed_m_per_ns must be finite and >= 0");
    }
    return trigger;
}

CameraGeometry buildCameraGeometry(const CameraConfig& cfg) {
    if (!cfg.enabled) {
        return CameraGeometry{};
    }
    const std::string mode = lowerCopy(trim(cfg.mode));
    if (mode == "csv" || mode == "imported") {
        if (cfg.csv_path.empty()) {
            throw std::runtime_error("camera.csv_path is required when camera.mode=csv");
        }
        return readCameraCsv(cfg.csv_path);
    }
    auto camera = makeGeneratedCamera(cfg);
    validateCameraGeometry(camera, "generated camera");
    return camera;
}

double cameraPixelSizeForCollector(const CameraConfig& cfg, const CameraGeometry& camera)
{
    if (cfg.collector_entrance_size_m > 0.0) {
        return cfg.collector_entrance_size_m;
    }
    if (!camera.empty()) {
        return camera.pixels().front().size;
    }
    return cfg.pixel_size_m;
}

std::unique_ptr<Cone::SquareCone> buildLightCollector(const CameraConfig& cfg,
                                                      const CameraGeometry& camera)
{
    const std::string collector = lowerCopy(trim(cfg.collector));
    if (isDisabledText(collector)) {
        return nullptr;
    }
    if (camera.empty()) {
        throw std::runtime_error(
            "camera.collector requires an enabled, non-empty camera");
    }
    const double reference_size = camera.pixels().front().size;
    for (const auto& pixel : camera.pixels()) {
        if (pixel.shape != PixelShape::Square) {
            throw std::runtime_error(
                "square-cone camera.collector currently requires square pixels; "
                "pixel " + std::to_string(pixel.id) + " is not square");
        }
        if (std::abs(pixel.size - reference_size) >
            std::max(1.0e-12, reference_size * 1.0e-9)) {
            throw std::runtime_error(
                "camera.collector requires all pixels to have one common size");
        }
    }

    // 这里仅把 output-plane 局部坐标转换为单个 square cone 的局部坐标。
    std::unique_ptr<Cone::SquareCone> cone;
    const double entrance_size_m = cameraPixelSizeForCollector(cfg, camera);
    if (!std::isfinite(entrance_size_m) || entrance_size_m <= 0.0 ||
        !std::isfinite(cfg.collector_exit_size_m) ||
        cfg.collector_exit_size_m <= 0.0 ||
        !std::isfinite(cfg.collector_height_m) ||
        cfg.collector_height_m <= 0.0) {
        throw std::runtime_error(
            "camera collector entrance, exit, and height must be finite and > 0");
    }
    if (entrance_size_m <= cfg.collector_exit_size_m) {
        throw std::runtime_error(
            "camera collector entrance size must exceed its exit size");
    }
    if (collector == "bezier" || collector == "square_cone_bezier") {
        cone.reset(Cone::create_cone_bezier_surface(
            entrance_size_m * 1000.0,
            cfg.collector_exit_size_m * 1000.0,
            cfg.collector_height_m * 1000.0,
            0.0));
    } else if (collector == "parabolic" || collector == "square_cone_parabolic") {
        if (std::abs(entrance_size_m - 0.0243992) > 1.0e-9 ||
            std::abs(cfg.collector_exit_size_m - 0.0150014) > 1.0e-9 ||
            std::abs(cfg.collector_height_m - 0.0252718) > 1.0e-9) {
            throw std::runtime_error(
                "parabolic collector has fixed geometry: entrance=0.0243992 m, "
                "exit=0.0150014 m, height=0.0252718 m");
        }
        cone.reset(Cone::create_cone_parabolic_surface());
    } else {
        throw std::runtime_error("unsupported camera.collector: " + cfg.collector);
    }

    if (!isDisabledText(cfg.collector_reflectivity_csv)) {
        cone->set_material(std::make_unique<Cone::TableReflectMaterial>(
            readCollectorReflectivityCsv(cfg.collector_reflectivity_csv)));
        return cone;
    }

    const std::string material = lowerCopy(trim(cfg.collector_material));
    if (material == "true_reflect") {
        cone->set_material(std::make_unique<Cone::TrueReflectMaterial>());
    } else if (material == "ideal" || material == "mirror_reflect" || material == "mirror") {
        cone->set_material(std::make_unique<Cone::MirrorReflectMaterial>());
    } else {
        throw std::runtime_error("unsupported camera.collector_material: " +
                                 cfg.collector_material);
    }
    return cone;
}

const CameraPixel* findContainingPixel(const CameraGeometry& camera, double x, double y)
{
    return camera.findContainingPixelPtr(x, y);
}

CollectorTraceResult traceLightCollector(const Cone::SquareCone& cone,
                                         const OutputPlane& plane,
                                         const OpticalSurfaceHit& hit,
                                         const CameraPixel& pixel,
                                         const SipmConfig& sipm)
{
    Cone::Position pos{
        (hit.u_m - pixel.center.x) * 1000.0,
        (hit.v_m - pixel.center.y) * 1000.0,
        cone.height_ + cone.height_start_
    };
    Cone::DirectionVecter dir{
        hit.out_dir.dot(plane.u_axis),
        hit.out_dir.dot(plane.v_axis),
        hit.out_dir.dot(plane.normal.normalized())
    };
    double intensity = 1.0;
    const auto trace_result = cone.ray_trace_impl(pos, dir, intensity);
    const auto& exit_pos = trace_result.position;
    const auto& exit_dir = trace_result.direction;
    const double sipm_half_mm = 0.5 * sipm.size_m * 1000.0;
    const bool hits_sipm =
        std::abs(exit_pos.x_) <= sipm_half_mm + 1e-9 &&
        std::abs(exit_pos.y_) <= sipm_half_mm + 1e-9;

    CollectorTraceResult result;
    result.hit_sipm = trace_result.exits_collector && hits_sipm &&
                      trace_result.intensity > 0.0;
    result.intensity = trace_result.intensity;
    result.reflections = trace_result.reflections;
    result.path_length_m = trace_result.path_length * 0.001;
    result.reflection_limit_reached =
        trace_result.reflection_limit_reached;
    result.exit_position = exit_pos;
    result.exit_direction = exit_dir;
    return result;
}

void applyCameraResponse(const CameraGeometry& camera,
                         const Cone::SquareCone* light_collector,
                         const OutputPlane& plane,
                         const SipmConfig& sipm,
                         const ElectronicsResponse& electronics,
                         OpticalSurfaceHit& hit,
                         double speed_of_light_m_per_ns,
                         const ::lact::electronics::MicrocellConfig* microcell,
                         double post_geometry_pde_scale)
{
    if (!std::isfinite(speed_of_light_m_per_ns) ||
        speed_of_light_m_per_ns <= 0.0) {
        throw std::runtime_error(
            "collector propagation speed must be finite and > 0");
    }
    if (!std::isfinite(post_geometry_pde_scale) ||
        post_geometry_pde_scale <= 0.0) {
        throw std::runtime_error(
            "post-geometry PDE scale must be finite and > 0");
    }
    hit.camera_enabled = true;
    const double front_face_cosine =
        hit.out_dir.dot(plane.normal.normalized());
    if (!std::isfinite(front_face_cosine) || front_face_cosine >= -1.0e-12) {
        hit.hit_camera = false;
        hit.accepted = false;
        hit.pixel_id = -1;
        return;
    }
    hit.camera_x_m = hit.u_m;
    hit.camera_y_m = hit.v_m;
    const CameraPixel* pixel = findContainingPixel(camera, hit.camera_x_m, hit.camera_y_m);
    hit.pixel_id = pixel ? pixel->id : -1;
    hit.hit_camera = pixel != nullptr;

    if (pixel && light_collector) {
        hit.collector_enabled = true;
        auto collector_result = traceLightCollector(*light_collector, plane, hit, *pixel, sipm);
        hit.hit_collector = collector_result.hit_sipm;
        hit.hit_camera = collector_result.hit_sipm;
        hit.collector_reflections = collector_result.reflections;
        hit.collector_intensity = collector_result.intensity;
        hit.collector_path_length_m = collector_result.path_length_m;
        hit.collector_time_delay_ns =
            collector_result.path_length_m / speed_of_light_m_per_ns;
        hit.collector_reflection_limit_reached =
            collector_result.reflection_limit_reached;
        hit.collector_exit_x_m = collector_result.exit_position.x_ * 0.001;
        hit.collector_exit_y_m = collector_result.exit_position.y_ * 0.001;
        hit.collector_exit_z_m = collector_result.exit_position.z_ * 0.001;
        hit.collector_dir_u = collector_result.exit_direction.x_;
        hit.collector_dir_v = collector_result.exit_direction.y_;
        hit.collector_dir_w = collector_result.exit_direction.z_;
        hit.relative_efficiency *= collector_result.intensity;
        if (collector_result.hit_sipm) {
            hit.time_ns += hit.collector_time_delay_ns;
        }
    }

    double pde = electronics.peConversion(hit.wavelength_nm);
    if (hit.hit_camera && microcell && microcell->enabled &&
        microcell->layout == "s17351_tiled_2x4") {
        if (!light_collector && hit.collector_exit_x_m == 0.0 &&
            hit.collector_exit_y_m == 0.0 &&
            hit.collector_exit_z_m == 0.0) {
            throw std::runtime_error(
                "explicit S17351 microcell geometry has no light-collector "
                "exit coordinates; configure camera.collector or disable the "
                "explicit microcell geometry");
        }
        hit.sipm_geometry_enabled = true;
        const auto address = ::lact::electronics::mapMicrocellPosition(
            *microcell, hit.collector_exit_x_m, hit.collector_exit_y_m);
        hit.sipm_grid_column = address.grid_column;
        hit.sipm_grid_row = address.grid_row;
        hit.sipm_channel_id = address.channel_id;
        hit.sipm_microcell_id = address.microcell_id;
        if (!address.inside_channel) {
            hit.sipm_gap_rejected = true;
            hit.sipm_channel_gap_rejected = address.channel_gap;
            hit.hit_camera = false;
            hit.accepted = false;
            hit.relative_efficiency = 0.0;
            return;
        }
        // In expectation mode the datasheet PDE is still carried by the
        // optical hit weight.  If that PDE was averaged over the full tiled
        // package, condition it on being inside a real channel before the
        // explicit channel-gap geometry is applied.
        hit.relative_efficiency *= post_geometry_pde_scale;
    }
    if (hit.hit_camera) {
        hit.relative_efficiency *= pde;
    }
    hit.accepted = hit.hit_camera && hit.relative_efficiency > 0.0;
}

void accumulatePixelHit(std::map<PixelKey, PixelAccumulator>& pixels,
                        int event_id,
                        int telescope_id,
                        const OpticalSurfaceHit& hit)
{
    if (!hit.hit_camera || hit.pixel_id < 0) {
        return;
    }
    const double signal = hit.weight * hit.relative_efficiency;
    auto& acc = pixels[{event_id, telescope_id, hit.pixel_id}];
    acc.event_id = event_id;
    acc.telescope_id = telescope_id;
    acc.pixel_id = hit.pixel_id;
    acc.photon_count += 1;
    acc.pe += signal;
    acc.signal += signal;
    acc.time_sum += signal * hit.time_ns;
    acc.time2_sum += signal * hit.time_ns * hit.time_ns;
}

void writePixelCsv(const std::string& path,
                   const std::map<PixelKey, PixelAccumulator>& pixels)
{
    const std::filesystem::path out_path(path);
    if (out_path.has_parent_path()) {
        std::filesystem::create_directories(out_path.parent_path());
    }
    std::ofstream ofs(path);
    if (!ofs) {
        throw std::runtime_error("failed to write pixel CSV: " + path);
    }
    ofs << std::setprecision(10);
    ofs << "event_id,telescope_id,pixel_id,photon_count,pe,signal,"
        << "time_mean_ns,time_rms_ns\n";
    for (const auto& kv : pixels) {
        const auto& p = kv.second;
        const double mean = p.signal > 0.0 ? p.time_sum / p.signal : 0.0;
        const double var = p.signal > 0.0
            ? std::max(0.0, p.time2_sum / p.signal - mean * mean)
            : 0.0;
        ofs << p.event_id << ","
            << p.telescope_id << ","
            << p.pixel_id << ","
            << p.photon_count << ","
            << p.pe << ","
            << p.signal << ","
            << mean << ","
            << std::sqrt(var) << "\n";
    }
}

TelescopeConfig buildTelescopeConfig(const std::map<std::string, std::string>& cfg) {
    TelescopeConfig telescope;
    telescope.id = getInt(cfg, "telescope.id", telescope.id);
    telescope.name = getString(cfg, "telescope.name", telescope.name);
    telescope.position_m = getVec3(cfg, "telescope.position_m", telescope.position_m);
    telescope.pointing_az_deg =
        getDouble(cfg, "telescope.pointing_az_deg", telescope.pointing_az_deg);
    telescope.pointing_el_deg =
        getDouble(cfg, "telescope.pointing_el_deg", telescope.pointing_el_deg);
    telescope.focal_length_m =
        getDouble(cfg, "telescope.focal_length_m", telescope.focal_length_m);
    telescope.coordinate_system =
        getString(cfg, "telescope.coordinate_system", telescope.coordinate_system);

    if (!std::isfinite(telescope.position_m.x) ||
        !std::isfinite(telescope.position_m.y) ||
        !std::isfinite(telescope.position_m.z)) {
        throw std::runtime_error("telescope.position_m must contain finite coordinates");
    }
    if (!std::isfinite(telescope.pointing_az_deg) ||
        !std::isfinite(telescope.pointing_el_deg)) {
        throw std::runtime_error("telescope pointing angles must be finite");
    }
    if (!std::isfinite(telescope.focal_length_m) || telescope.focal_length_m <= 0.0) {
        throw std::runtime_error("telescope.focal_length_m must be finite and > 0");
    }
    return telescope;
}

ErrorConfig buildErrorConfig(const std::map<std::string, std::string>& cfg) {
    ErrorConfig error;
    error.random_seed = getUInt64(cfg, "error.random_seed", error.random_seed);
    error.facet_radial_position_sigma_m =
        getDouble(cfg, "error.facet_radial_position_sigma_m", 0.0);
    error.facet_normal_sigma_deg =
        getDouble(cfg, "error.facet_normal_sigma_deg", 0.0);
    error.reflect_direction_sigma_deg =
        getDouble(cfg, "error.reflect_direction_sigma_deg", 0.0);
    error.radius_of_curvature_sigma_m =
        getDouble(cfg, "error.radius_of_curvature_sigma_m", 0.0);
    error.reflectivity_scale_sigma =
        getDouble(cfg, "error.reflectivity_scale_sigma", 0.0);
    error.structural_deformation_config =
        getString(cfg, "error.structural_deformation_config", "");

    auto checkNonNegative = [](double value, const std::string& key) {
        if (!std::isfinite(value) || value < 0.0) {
            throw std::runtime_error(key + " must be finite and >= 0");
        }
    };
    checkNonNegative(error.facet_radial_position_sigma_m,
                     "error.facet_radial_position_sigma_m");
    checkNonNegative(error.facet_normal_sigma_deg,
                     "error.facet_normal_sigma_deg");
    checkNonNegative(error.reflect_direction_sigma_deg,
                     "error.reflect_direction_sigma_deg");
    checkNonNegative(error.radius_of_curvature_sigma_m,
                     "error.radius_of_curvature_sigma_m");
    checkNonNegative(error.reflectivity_scale_sigma,
                     "error.reflectivity_scale_sigma");
    return error;
}

bool ObstructionMask::contains(double x_m, double y_m) const
{
    if (!enabled || nx <= 0 || ny <= 0 || cell_size_m <= 0.0 ||
        blocked.size() != static_cast<std::size_t>(nx * ny)) {
        return false;
    }
    const int ix = static_cast<int>(std::floor((x_m - x_min_m) / cell_size_m));
    const int iy = static_cast<int>(std::floor((y_m - y_min_m) / cell_size_m));
    if (ix < 0 || iy < 0 || ix >= nx || iy >= ny) {
        return false;
    }
    return blocked[static_cast<std::size_t>(iy * nx + ix)] != 0;
}

namespace {

bool segmentIntersectsCylinder(const Vec3& a,
                               const Vec3& b,
                               const Vec3& c0,
                               const Vec3& c1,
                               double radius)
{
    if (radius <= 0.0) {
        return false;
    }
    const Vec3 u = b - a;
    const Vec3 v = c1 - c0;
    const Vec3 w0 = a - c0;
    const double uu = u.dot(u);
    const double vv = v.dot(v);
    if (uu <= 0.0 || vv <= 0.0) {
        return false;
    }
    const double uv = u.dot(v);
    const double uw = u.dot(w0);
    const double vw = v.dot(w0);
    const double ww = w0.dot(w0);
    const double denom = uu * vv - uv * uv;

    auto distance2 = [&](double s, double t) {
        const Vec3 d = w0 + u * s - v * t;
        return d.dot(d);
    };

    double best = std::numeric_limits<double>::max();
    if (std::abs(denom) > 1e-14) {
        const double s = std::clamp((uv * vw - vv * uw) / denom, 0.0, 1.0);
        const double t = std::clamp((uu * vw - uv * uw) / denom, 0.0, 1.0);
        best = std::min(best, distance2(s, t));
    }

    const double candidates_s[] = {0.0, 1.0};
    for (double s : candidates_s) {
        const double t = std::clamp((vw + uv * s) / vv, 0.0, 1.0);
        best = std::min(best, distance2(s, t));
    }
    const double candidates_t[] = {0.0, 1.0};
    for (double t : candidates_t) {
        const double s = std::clamp((uv * t - uw) / uu, 0.0, 1.0);
        best = std::min(best, distance2(s, t));
    }

    // Include spherical end caps. This intentionally models a solid tube with
    // rounded ends, which is conservative for simplified telescope struts.
    return best <= radius * radius;
}

bool segmentIntersectsAabb(const Vec3& a,
                           const Vec3& b,
                           const Vec3& center,
                           const Vec3& half,
                           double* out_t0 = nullptr,
                           double* out_t1 = nullptr)
{
    const Vec3 d = b - a;
    double tmin = 0.0;
    double tmax = 1.0;
    const double amin[3] = {center.x - half.x, center.y - half.y, center.z - half.z};
    const double amax[3] = {center.x + half.x, center.y + half.y, center.z + half.z};
    const double p[3] = {a.x, a.y, a.z};
    const double dir[3] = {d.x, d.y, d.z};
    for (int i = 0; i < 3; ++i) {
        if (std::abs(dir[i]) < 1e-14) {
            if (p[i] < amin[i] || p[i] > amax[i]) {
                return false;
            }
            continue;
        }
        double t1 = (amin[i] - p[i]) / dir[i];
        double t2 = (amax[i] - p[i]) / dir[i];
        if (t1 > t2) {
            std::swap(t1, t2);
        }
        tmin = std::max(tmin, t1);
        tmax = std::min(tmax, t2);
        if (tmin > tmax) {
            return false;
        }
    }
    if (out_t0) {
        *out_t0 = tmin;
    }
    if (out_t1) {
        *out_t1 = tmax;
    }
    return true;
}

std::vector<Vec3> regularPolygonVertices(double radius,
                                          double rotation,
                                          int sides,
                                          double z = 0.0)
{
    std::vector<Vec3> vertices;
    if (radius <= 0.0 || sides < 3) {
        return vertices;
    }
    vertices.reserve(static_cast<std::size_t>(sides));
    for (int i = 0; i < sides; ++i) {
        constexpr double kPi = 3.141592653589793238462643383279502884;
        const double theta = rotation + 2.0 * kPi * static_cast<double>(i) /
                                            static_cast<double>(sides);
        vertices.push_back({radius * std::cos(theta), radius * std::sin(theta), z});
    }
    return vertices;
}

bool pointInsideRegularPolygon(const Vec3& p,
                               const Vec3& center,
                               double radius,
                               double rotation,
                               int sides)
{
    const auto vertices = regularPolygonVertices(radius, rotation, sides);
    if (vertices.empty()) {
        return false;
    }
    bool has_pos = false;
    bool has_neg = false;
    const double x = p.x - center.x;
    const double y = p.y - center.y;
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        const Vec3& a = vertices[i];
        const Vec3& b = vertices[(i + 1) % vertices.size()];
        const double cross = (b.x - a.x) * (y - a.y) -
                             (b.y - a.y) * (x - a.x);
        if (cross > 1e-12) {
            has_pos = true;
        } else if (cross < -1e-12) {
            has_neg = true;
        }
        if (has_pos && has_neg) {
            return false;
        }
    }
    return true;
}

bool segmentIntersectsRegularPrism(const Vec3& a,
                                   const Vec3& b,
                                   const Vec3& center,
                                   double radius,
                                   double height,
                                   double rotation,
                                   int sides)
{
    if (radius <= 0.0 || height <= 0.0 || sides < 3) {
        return false;
    }
    const Vec3 half{radius, radius, 0.5 * height};
    double t0 = 0.0;
    double t1 = 0.0;
    if (!segmentIntersectsAabb(a, b, center, half, &t0, &t1)) {
        return false;
    }
    const Vec3 d = b - a;
    const Vec3 p0 = a + d * t0;
    const Vec3 p1 = a + d * t1;
    const Vec3 pm = a + d * (0.5 * (t0 + t1));
    return pointInsideRegularPolygon(p0, center, radius, rotation, sides) ||
           pointInsideRegularPolygon(p1, center, radius, rotation, sides) ||
           pointInsideRegularPolygon(pm, center, radius, rotation, sides);
}

bool segmentIntersectsBoxWithHole(const Vec3& a,
                                  const Vec3& b,
                                  const ObstructionMask::Primitive& primitive)
{
    double t0 = 0.0;
    double t1 = 0.0;
    if (!segmentIntersectsAabb(a, b, primitive.center, primitive.half_size,
                               &t0, &t1)) {
        return false;
    }
    if (!primitive.has_hole) {
        return true;
    }
    const Vec3 d = b - a;
    const Vec3 entry = a + d * t0;
    const Vec3 exit = a + d * t1;
    const bool entry_inside = pointInsideRegularPolygon(entry,
                                                        primitive.center,
                                                        primitive.hole_radius_m,
                                                        primitive.hole_rotation_rad,
                                                        primitive.hole_sides);
    const bool exit_inside = pointInsideRegularPolygon(exit,
                                                       primitive.center,
                                                       primitive.hole_radius_m,
                                                       primitive.hole_rotation_rad,
                                                       primitive.hole_sides);
    // The hole is a convex vertical prism. If both endpoints of the segment
    // inside the AABB are inside the aperture, the whole segment stays in the
    // aperture and the box should not block it.
    return !(entry_inside && exit_inside);
}

enum class ObstructionLeg {
    Incoming,
    Reflected,
};

bool primitiveAppliesToLeg(const ObstructionMask::Primitive& primitive,
                           ObstructionLeg leg)
{
    const std::string role = lowerCopy(primitive.role);
    const std::string material = lowerCopy(primitive.material_id);
    if (material == "void") {
        return false;
    }
    return leg == ObstructionLeg::Incoming || role == "support_strut";
}

std::map<std::string, std::string> csvRowMap(const std::vector<std::string>& header,
                                             const std::vector<std::string>& cells)
{
    std::map<std::string, std::string> row;
    const std::size_t n = std::min(header.size(), cells.size());
    for (std::size_t i = 0; i < n; ++i) {
        row[lowerCopy(trim(header[i]))] = trim(cells[i]);
    }
    return row;
}

std::string csvGetString(const std::map<std::string, std::string>& row,
                         const std::string& key,
                         const std::string& fallback = "")
{
    const auto it = row.find(lowerCopy(key));
    if (it == row.end() || trim(it->second).empty()) {
        return fallback;
    }
    return it->second;
}

double csvGetDouble(const std::map<std::string, std::string>& row,
                    const std::string& key,
                    double fallback = 0.0)
{
    const std::string value = csvGetString(row, key, "");
    if (value.empty()) {
        return fallback;
    }
    return std::stod(value);
}

int csvGetInt(const std::map<std::string, std::string>& row,
              const std::string& key,
              int fallback = 0)
{
    const std::string value = csvGetString(row, key, "");
    if (value.empty()) {
        return fallback;
    }
    return std::stoi(value);
}

bool segmentBlockedLocal(const Vec3& a,
                         const Vec3& b,
                         const ObstructionMask& obstruction,
                         ObstructionLeg leg)
{
    if (!obstruction.enabled) {
        return false;
    }
    if (lowerCopy(obstruction.mode) == "primitives") {
        for (const auto& primitive : obstruction.primitives) {
            if (!primitiveAppliesToLeg(primitive, leg)) {
                continue;
            }
            const std::string type = lowerCopy(primitive.type);
            if (type == "cylinder") {
                if (segmentIntersectsCylinder(a, b, primitive.p0, primitive.p1,
                                              primitive.radius_m)) {
                    return true;
                }
            } else if (type == "box" || type == "aabb") {
                if (segmentIntersectsBoxWithHole(a, b, primitive)) {
                    return true;
                }
            } else if (type == "polygon_prism") {
                if (segmentIntersectsRegularPrism(a, b, primitive.center,
                                                  primitive.radius_m,
                                                  primitive.height_m,
                                                  primitive.rotation_rad,
                                                  primitive.sides)) {
                    return true;
                }
            }
        }
        return false;
    }
    const Vec3 d = b - a;
    if (std::abs(d.z) < 1e-12) {
        return false;
    }
    const double s = (obstruction.plane_z_m - a.z) / d.z;
    if (s < 0.0 || s > 1.0) {
        return false;
    }
    const Vec3 q = a + d * s;
    return obstruction.contains(q.x, q.y);
}

} // namespace

ObstructionMask buildObstructionMask(const std::map<std::string, std::string>& cfg)
{
    ObstructionMask obstruction;
    obstruction.enabled = getBool(cfg, "obstruction.enabled", false);
    obstruction.mode = lowerCopy(getString(cfg, "obstruction.mode", "mask"));
    obstruction.mask_csv = getString(cfg, "obstruction.mask_csv", "");
    obstruction.primitives_csv = getString(cfg, "obstruction.primitives_csv", "");
    obstruction.check_incoming = getBool(cfg, "obstruction.check_incoming", true);
    obstruction.check_reflected = getBool(cfg, "obstruction.check_reflected",
                                          obstruction.mode == "primitives");
    obstruction.mark_only = getBool(cfg, "obstruction.mark_only", false);
    obstruction.plane_z_m = getDouble(cfg, "obstruction.plane_z_m", obstruction.plane_z_m);

    if (!obstruction.enabled) {
        return obstruction;
    }
    if (obstruction.mode == "primitives") {
        if (isDisabledText(obstruction.primitives_csv)) {
            throw std::runtime_error(
                "obstruction.mode=primitives requires obstruction.primitives_csv");
        }
        std::ifstream in(obstruction.primitives_csv);
        if (!in) {
            throw std::runtime_error("failed to read obstruction primitives CSV: " +
                                     obstruction.primitives_csv);
        }
        std::string line;
        std::vector<std::string> header;
        while (std::getline(in, line)) {
            line = trim(line);
            if (line.empty() || line[0] == '#') {
                continue;
            }
            const auto cells = splitCsvCells(line);
            if (header.empty()) {
                header = cells;
                continue;
            }
            const auto row = csvRowMap(header, cells);
            if (row.empty()) {
                continue;
            }
            ObstructionMask::Primitive p;
            p.type = lowerCopy(csvGetString(row, "type"));
            p.name = csvGetString(row, "name", csvGetString(row, "role", p.type));
            p.role = csvGetString(row, "role", "default");
            p.material_id = csvGetString(row, "material_id", "default");
            p.p0 = {csvGetDouble(row, "x0_m"),
                    csvGetDouble(row, "y0_m"),
                    csvGetDouble(row, "z0_m")};
            p.p1 = {csvGetDouble(row, "x1_m", p.p0.x),
                    csvGetDouble(row, "y1_m", p.p0.y),
                    csvGetDouble(row, "z1_m", p.p0.z)};
            p.center = {csvGetDouble(row, "center_x_m", p.p0.x),
                        csvGetDouble(row, "center_y_m", p.p0.y),
                        csvGetDouble(row, "center_z_m", p.p0.z)};
            p.radius_m = csvGetDouble(row, "radius_m");
            p.height_m = csvGetDouble(row, "height_m");
            p.rotation_rad = csvGetDouble(row, "rotation_rad");
            p.sides = csvGetInt(row, "sides");
            p.half_size = {csvGetDouble(row, "half_x_m"),
                           csvGetDouble(row, "half_y_m"),
                           csvGetDouble(row, "half_z_m")};
            p.bbox_min = {csvGetDouble(row, "bbox_min_x_m", p.center.x - p.half_size.x),
                          csvGetDouble(row, "bbox_min_y_m", p.center.y - p.half_size.y),
                          csvGetDouble(row, "bbox_min_z_m", p.center.z - p.half_size.z)};
            p.bbox_max = {csvGetDouble(row, "bbox_max_x_m", p.center.x + p.half_size.x),
                          csvGetDouble(row, "bbox_max_y_m", p.center.y + p.half_size.y),
                          csvGetDouble(row, "bbox_max_z_m", p.center.z + p.half_size.z)};
            if (p.half_size.norm2() <= 0.0 &&
                (p.bbox_max - p.bbox_min).norm2() > 0.0) {
                p.center = (p.bbox_min + p.bbox_max) * 0.5;
                p.half_size = (p.bbox_max - p.bbox_min) * 0.5;
            }
            p.hole_radius_m = csvGetDouble(row, "hole_radius_m");
            p.hole_rotation_rad = csvGetDouble(row, "hole_rotation_rad");
            p.hole_sides = csvGetInt(row, "hole_sides");
            p.has_hole = p.hole_radius_m > 0.0 && p.hole_sides >= 3;
            obstruction.primitives.push_back(p);
        }
        if (obstruction.primitives.empty()) {
            throw std::runtime_error("obstruction primitives CSV has no primitives: " +
                                     obstruction.primitives_csv);
        }
        return obstruction;
    }
    if (isDisabledText(obstruction.mask_csv)) {
        throw std::runtime_error("obstruction.enabled=true requires obstruction.mask_csv");
    }

    std::ifstream in(obstruction.mask_csv);
    if (!in) {
        throw std::runtime_error("failed to read obstruction mask CSV: " +
                                 obstruction.mask_csv);
    }

    std::string line;
    bool saw_header = false;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) {
            continue;
        }
        if (line[0] == '#') {
            const std::string meta = trim(line.substr(1));
            const auto eq = meta.find('=');
            if (eq == std::string::npos) {
                continue;
            }
            const std::string key = lowerCopy(trim(meta.substr(0, eq)));
            const std::string value = trim(meta.substr(eq + 1));
            if (key == "x_min_m") {
                obstruction.x_min_m = std::stod(value);
            } else if (key == "y_min_m") {
                obstruction.y_min_m = std::stod(value);
            } else if (key == "cell_size_m") {
                obstruction.cell_size_m = std::stod(value);
            } else if (key == "nx") {
                obstruction.nx = std::stoi(value);
            } else if (key == "ny") {
                obstruction.ny = std::stoi(value);
            } else if (key == "plane_z_m") {
                obstruction.plane_z_m = std::stod(value);
            }
            continue;
        }

        const auto cells = splitCsvCells(line);
        if (!saw_header) {
            saw_header = true;
            continue;
        }
        if (obstruction.nx <= 0 || obstruction.ny <= 0) {
            throw std::runtime_error("obstruction mask CSV must define # nx and # ny");
        }
        if (obstruction.blocked.empty()) {
            obstruction.blocked.assign(static_cast<std::size_t>(obstruction.nx *
                                                                obstruction.ny),
                                       0);
        }
        if (cells.size() < 2) {
            continue;
        }
        const int ix = std::stoi(cells[0]);
        const int iy = std::stoi(cells[1]);
        if (ix < 0 || iy < 0 || ix >= obstruction.nx || iy >= obstruction.ny) {
            continue;
        }
        obstruction.blocked[static_cast<std::size_t>(iy * obstruction.nx + ix)] = 1;
    }

    if (!std::isfinite(obstruction.x_min_m) ||
        !std::isfinite(obstruction.y_min_m) ||
        !std::isfinite(obstruction.cell_size_m) ||
        obstruction.cell_size_m <= 0.0 ||
        obstruction.nx <= 0 || obstruction.ny <= 0) {
        throw std::runtime_error("invalid obstruction mask metadata in " +
                                 obstruction.mask_csv);
    }
    if (obstruction.blocked.empty()) {
        obstruction.blocked.assign(static_cast<std::size_t>(obstruction.nx *
                                                            obstruction.ny),
                                   0);
    }
    return obstruction;
}

bool photonBlockedByObstruction(const Photon& photon,
                                const ObstructionMask& obstruction,
                                const TelescopeFrame* trace_to_local_frame)
{
    if (!obstruction.enabled || !obstruction.check_incoming) {
        return false;
    }
    Vec3 p = photon.pos;
    Vec3 d = photon.dir;
    if (trace_to_local_frame) {
        p = trace_to_local_frame->pointToLocal(photon.pos);
        d = trace_to_local_frame->rotateVectorToLocal(photon.dir).normalized();
    }
    return segmentBlockedLocal(p, p + d * 1.0e4, obstruction,
                               ObstructionLeg::Incoming);
}

bool incomingSegmentBlockedByObstruction(const Vec3& a,
                                         const Vec3& b,
                                         const ObstructionMask& obstruction,
                                         const TelescopeFrame* trace_to_local_frame)
{
    if (!obstruction.enabled || !obstruction.check_incoming) {
        return false;
    }
    Vec3 p0 = a;
    Vec3 p1 = b;
    if (trace_to_local_frame) {
        p0 = trace_to_local_frame->pointToLocal(a);
        p1 = trace_to_local_frame->pointToLocal(b);
    }
    return segmentBlockedLocal(p0, p1, obstruction, ObstructionLeg::Incoming);
}

bool incomingRayBlockedByObstruction(const Vec3& mirror_point,
                                     const Vec3& incoming_direction,
                                     const ObstructionMask& obstruction,
                                     const TelescopeFrame* trace_to_local_frame)
{
    if (!obstruction.enabled || !obstruction.check_incoming) {
        return false;
    }

    Vec3 p = mirror_point;
    Vec3 d = incoming_direction;
    if (trace_to_local_frame) {
        p = trace_to_local_frame->pointToLocal(mirror_point);
        d = trace_to_local_frame->rotateVectorToLocal(incoming_direction);
    }
    if (d.norm2() <= 0.0) {
        return false;
    }
    d = d.normalized();

    double upstream_length = 0.0;
    if (lowerCopy(obstruction.mode) == "primitives") {
        double obstruction_radius = 0.0;
        for (const auto& primitive : obstruction.primitives) {
            const std::string type = lowerCopy(primitive.type);
            if (type == "cylinder") {
                obstruction_radius = std::max(
                    obstruction_radius,
                    std::max(primitive.p0.norm(), primitive.p1.norm()) +
                        std::max(0.0, primitive.radius_m));
            } else if (type == "box" || type == "aabb") {
                obstruction_radius = std::max(
                    obstruction_radius,
                    primitive.center.norm() + primitive.half_size.norm());
            } else if (type == "polygon_prism") {
                const double half_height = 0.5 * std::max(0.0, primitive.height_m);
                const double local_radius =
                    std::hypot(std::max(0.0, primitive.radius_m), half_height);
                obstruction_radius = std::max(
                    obstruction_radius,
                    primitive.center.norm() + local_radius);
            }
        }
        upstream_length = p.norm() + obstruction_radius + 1.0;
    } else {
        if (std::abs(d.z) < 1e-12) {
            return false;
        }
        const double distance_to_plane = (p.z - obstruction.plane_z_m) / d.z;
        if (distance_to_plane <= 0.0) {
            return false;
        }
        upstream_length = distance_to_plane + 1.0;
    }

    return segmentBlockedLocal(p - d * upstream_length, p, obstruction,
                               ObstructionLeg::Incoming);
}

bool segmentBlockedByObstruction(const Vec3& a,
                                 const Vec3& b,
                                 const ObstructionMask& obstruction,
                                 const TelescopeFrame* trace_to_local_frame)
{
    if (!obstruction.enabled || !obstruction.check_reflected) {
        return false;
    }
    Vec3 p0 = a;
    Vec3 p1 = b;
    if (trace_to_local_frame) {
        p0 = trace_to_local_frame->pointToLocal(a);
        p1 = trace_to_local_frame->pointToLocal(b);
    }
    return segmentBlockedLocal(p0, p1, obstruction, ObstructionLeg::Reflected);
}

void applyStructuralDeformation(std::vector<MirrorFacet>& facets,
                                const ErrorConfig& error,
                                const TelescopeConfig& telescope)
{
    if (isDisabledText(error.structural_deformation_config)) {
        return;
    }

    auto deformation_cfg = loadScopedMirrorConfig(error.structural_deformation_config);
    deformation_cfg["telescope.pointing_el_deg"] = doubleToString(telescope.pointing_el_deg);
    deformation_cfg["mirror.series_elevation_deg"] =
        getString(deformation_cfg, "mirror.series_elevation_deg",
                  doubleToString(telescope.pointing_el_deg));

    const std::string mode = lowerCopy(getString(deformation_cfg, "mirror.mode", ""));
    if (!(mode == "elevation_series" || mode == "series")) {
        throw std::runtime_error(
            "error.structural_deformation_config must point to a mirror elevation_series config");
    }

    const auto deformation_facets = buildElevationSeriesFacets(deformation_cfg);
    if (deformation_facets.size() != facets.size()) {
        throw std::runtime_error(
            "structural deformation facet count does not match base mirror layout");
    }

    for (std::size_t i = 0; i < facets.size(); ++i) {
        facets[i].center = deformation_facets[i].center;
        facets[i].normal = deformation_facets[i].normal;
    }
}

namespace {

bool perFacetSigmaColumnEnabled(
    const std::vector<MirrorFacet>& facets,
    double MirrorFacet::* value_member,
    bool MirrorFacet::* presence_member,
    const std::string& column_name)
{
    bool any_positive = false;
    for (const auto& facet : facets) {
        const double value = facet.*value_member;
        if (!std::isfinite(value) || value < 0.0) {
            throw std::runtime_error(
                column_name + " must be finite and >= 0 for facet " +
                std::to_string(facet.id));
        }
        any_positive = any_positive || value > 0.0;
    }
    if (!any_positive) {
        return false;
    }

    for (const auto& facet : facets) {
        const double value = facet.*value_member;
        const bool explicitly_present = facet.*presence_member;
        // Positive programmatically constructed values predate the presence
        // flags and remain explicit for backward compatibility. A zero must
        // carry the flag so that it can be distinguished from a missing cell.
        if (!explicitly_present && value == 0.0) {
            throw std::runtime_error(
                "per-facet " + column_name +
                " is active but missing facet " + std::to_string(facet.id) +
                "; provide an explicit 0 or a non-zero value for every facet");
        }
    }
    return true;
}

} // namespace

void applyFacetErrors(std::vector<MirrorFacet>& facets, const ErrorConfig& error) {
    const bool has_per_facet_misalignment =
        perFacetSigmaColumnEnabled(
            facets,
            &MirrorFacet::misalign_sigma_rad,
            &MirrorFacet::has_misalign_sigma_rad,
            "misalign_sigma_rad");
    if (error.facet_radial_position_sigma_m == 0.0 &&
        error.facet_normal_sigma_deg == 0.0 &&
        error.radius_of_curvature_sigma_m == 0.0 &&
        error.reflectivity_scale_sigma == 0.0 &&
        !has_per_facet_misalignment) {
        return;
    }

    // A positive per-facet column activates the whole column. In that mode zero
    // is an explicit zero for that facet and the uniform sigma is not added.
    // If the column is absent or all zeros, fall back to errors.cfg.
    const double normal_sigma_rad = has_per_facet_misalignment
        ? 0.0
        : error.facet_normal_sigma_deg * DEG_TO_RAD;

    for (auto& facet : facets) {
        Vec3 ideal_normal = facet.normal.normalized();
        const auto facet_seed = deriveSeed(
            error.random_seed,
            static_cast<std::uint64_t>(
                static_cast<std::int64_t>(facet.id) + 1));
        auto normal_deviate = [&](RandomStage stage) {
            std::mt19937_64 rng(stageSeed(facet_seed, stage));
            std::normal_distribution<double> unit_normal(0.0, 1.0);
            return unit_normal(rng);
        };

        if (error.facet_radial_position_sigma_m > 0.0) {
            const double offset = error.facet_radial_position_sigma_m *
                normal_deviate(RandomStage::FacetPositionError);
            facet.center += ideal_normal * offset;
        }

        const double facet_normal_sigma_rad = has_per_facet_misalignment
            ? facet.misalign_sigma_rad
            : normal_sigma_rad;
        if (facet_normal_sigma_rad > 0.0) {
            std::mt19937_64 normal_rng(stageSeed(
                facet_seed, RandomStage::FacetNormalError));
            facet.normal = perturbVectorOnSphere(
                ideal_normal, facet_normal_sigma_rad, normal_rng);
        }

        if (error.radius_of_curvature_sigma_m > 0.0 &&
            std::isfinite(facet.radius_of_curvature)) {
            double radius = facet.radius_of_curvature +
                            error.radius_of_curvature_sigma_m *
                                normal_deviate(RandomStage::FacetRadiusError);
            facet.radius_of_curvature = std::max(1e-9, radius);
        }

        if (error.reflectivity_scale_sigma > 0.0) {
            double scale = 1.0 + error.reflectivity_scale_sigma *
                normal_deviate(RandomStage::FacetReflectivityError);
            facet.reflectivity_scale = std::clamp(
                facet.reflectivity_scale * scale, 0.0, 1.0);
        }
    }
}

void applyFacetEfficiencyScales(std::vector<MirrorFacet>& facets,
                                const std::map<std::string, std::string>& cfg)
{
    const double uniform_scale =
        getDouble(cfg, "efficiency.mirror_reflectivity_scale", 1.0);
    if (!std::isfinite(uniform_scale) ||
        uniform_scale < 0.0 || uniform_scale > 1.0) {
        throw std::runtime_error(
            "efficiency.mirror_reflectivity_scale must be finite and in [0,1]");
    }

    const std::string csv_path =
        getString(cfg, "efficiency.mirror_reflectivity_scale_csv", "");
    if (!isDisabledText(csv_path)) {
        std::ifstream ifs(csv_path);
        if (!ifs) {
            throw std::runtime_error(
                "failed to open per-facet reflectivity scale CSV: " + csv_path);
        }

        std::string line;
        if (!std::getline(ifs, line)) {
            throw std::runtime_error(
                "per-facet reflectivity scale CSV is empty: " + csv_path);
        }
        const auto header_cells = splitCsvCells(line);
        std::map<std::string, int> header;
        for (int i = 0; i < static_cast<int>(header_cells.size()); ++i) {
            header[lowerCopy(header_cells[i])] = i;
        }
        if (headerIndex(header, "id") < 0 ||
            headerIndex(header, "reflectivity_scale") < 0) {
            throw std::runtime_error(
                "per-facet reflectivity scale CSV requires id,reflectivity_scale");
        }

        std::map<int, double> scales;
        int line_no = 1;
        while (std::getline(ifs, line)) {
            ++line_no;
            if (trim(line).empty()) {
                continue;
            }
            const auto cells = splitCsvCells(line);
            const int id = std::stoi(getRequiredCell(cells, header, "id"));
            double scale = 0.0;
            const std::string text =
                getRequiredCell(cells, header, "reflectivity_scale");
            if (!parseDoubleStrict(text, scale) ||
                !std::isfinite(scale) || scale < 0.0 || scale > 1.0) {
                throw std::runtime_error(
                    "per-facet reflectivity scale CSV line " +
                    std::to_string(line_no) +
                    " requires reflectivity_scale in [0,1]");
            }
            if (!scales.emplace(id, scale).second) {
                throw std::runtime_error(
                    "per-facet reflectivity scale CSV has duplicate id " +
                    std::to_string(id));
            }
        }

        for (auto& facet : facets) {
            const auto it = scales.find(facet.id);
            if (it == scales.end()) {
                throw std::runtime_error(
                    "per-facet reflectivity scale CSV is missing facet id " +
                    std::to_string(facet.id));
            }
            facet.reflectivity_scale = it->second;
            facet.has_reflectivity_scale = true;
            scales.erase(it);
        }
        if (!scales.empty()) {
            throw std::runtime_error(
                "per-facet reflectivity scale CSV contains unknown facet id " +
                std::to_string(scales.begin()->first));
        }
        return;
    }

    // Legacy mirror CSV values remain valid. Missing/blank values use the
    // uniform efficiency value, whose program default is one.
    for (auto& facet : facets) {
        if (!facet.has_reflectivity_scale) {
            facet.reflectivity_scale = uniform_scale;
        }
    }
}

double effectiveReflectDirectionSigmaRad(const std::vector<MirrorFacet>& facets,
                                         const ErrorConfig& error)
{
    const bool has_positive_per_facet_roughness =
        perFacetSigmaColumnEnabled(
            facets,
            &MirrorFacet::roughness_sigma_rad,
            &MirrorFacet::has_roughness_sigma_rad,
            "roughness_sigma_rad");
    return has_positive_per_facet_roughness
        ? 0.0
        : error.reflect_direction_sigma_deg * DEG_TO_RAD;
}

OpticalEfficiencyConfig buildEfficiencyConfig(const std::map<std::string, std::string>& cfg) {
    OpticalEfficiencyConfig eff;
    eff.constant_scale = getDouble(cfg, "efficiency.constant_scale", eff.constant_scale);
    if (!std::isfinite(eff.constant_scale) || eff.constant_scale < 0.0) {
        throw std::runtime_error(
            "efficiency.constant_scale must be finite and >= 0");
    }

    eff.mirror_reflectivity = parseEfficiencyFactor(cfg, "efficiency.mirror_reflectivity");
    eff.filter_transmission = parseEfficiencyFactor(cfg, "efficiency.filter_transmission");
    eff.sipm_pde = parseFirstEfficiencyFactor(
        cfg,
        {"sipm.pde",
         "sipm.pe_conversion",
         "efficiency.sipm_pde",
         "electronics.pe_conversion",
         "electronics.pde"});
    eff.atmosphere_transmission = parseEfficiencyFactor(cfg, "efficiency.atmosphere_transmission");
    if (!eff.atmosphere_transmission.enabled) {
        eff.atmosphere_transmission = parseEfficiencyFactor(cfg, "atmosphere.transmission");
    }
    eff.use_funnel_acceptance =
        getBool(cfg, "efficiency.funnel_acceptance",
                getBool(cfg, "efficiency.use_funnel_acceptance", false));

    const std::string atmosphere_model =
        lowerCopy(trim(getString(cfg, "atmosphere.model", "none")));
    if (eff.atmosphere_transmission.enabled &&
        !isDisabledText(atmosphere_model) &&
        !getBool(cfg, "atmosphere.allow_additional_spectral_factor", false)) {
        throw std::runtime_error(
            "both atmosphere.model and a spectral atmosphere transmission "
            "factor are enabled; set "
            "atmosphere.allow_additional_spectral_factor=true only when the "
            "second factor is intentional");
    }

    const std::string collector =
        getString(cfg, "camera.collector", "none");
    if (eff.use_funnel_acceptance && !isDisabledText(collector) &&
        !getBool(cfg, "camera.allow_analytic_funnel_with_collector", false)) {
        throw std::runtime_error(
            "efficiency.funnel_acceptance and camera.collector would apply "
            "two angular collector models; disable one or explicitly set "
            "camera.allow_analytic_funnel_with_collector=true");
    }

    return eff;
}

PropagationConfig buildPropagationConfig(const std::map<std::string, std::string>& cfg) {
    PropagationConfig propagation;
    propagation.speed_of_light_m_per_ns =
        getDouble(cfg, "propagation.speed_of_light_m_per_ns",
                  propagation.speed_of_light_m_per_ns);
    if (!std::isfinite(propagation.speed_of_light_m_per_ns) ||
        propagation.speed_of_light_m_per_ns <= 0.0) {
        throw std::runtime_error("propagation.speed_of_light_m_per_ns must be finite and > 0");
    }
    return propagation;
}

} // namespace lact
