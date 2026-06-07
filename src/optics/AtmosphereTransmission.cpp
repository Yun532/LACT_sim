#include "optics/AtmosphereTransmission.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace {

std::string trimCopy(const std::string& s)
{
    const auto begin = s.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos) {
        return "";
    }
    const auto end = s.find_last_not_of(" \t\r\n");
    return s.substr(begin, end - begin + 1);
}

std::string lowerCopyLocal(std::string s)
{
    for (auto& c : s) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return s;
}

bool parseDoubleLocal(const std::string& text, double& out)
{
    try {
        std::size_t pos = 0;
        out = std::stod(text, &pos);
        return pos == text.size();
    } catch (...) {
        return false;
    }
}

bool disabledText(const std::string& text)
{
    const std::string value = lowerCopyLocal(trimCopy(text));
    return value.empty() || value == "none" || value == "off" ||
           value == "false" || value == "no";
}

double altitudeKmToLogCm(double altitude_km)
{
    if (!std::isfinite(altitude_km) || altitude_km <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::log10(altitude_km * 1.0e5);
}

std::string getStringLocal(const std::map<std::string, std::string>& cfg,
                           const std::string& key,
                           const std::string& fallback)
{
    auto it = cfg.find(key);
    return it == cfg.end() ? fallback : it->second;
}

double getDoubleLocal(const std::map<std::string, std::string>& cfg,
                      const std::string& key,
                      double fallback)
{
    auto it = cfg.find(key);
    if (it == cfg.end()) {
        return fallback;
    }
    std::size_t pos = 0;
    const double value = std::stod(it->second, &pos);
    if (pos != it->second.size()) {
        throw std::runtime_error("invalid numeric config value for " + key + ": " + it->second);
    }
    return value;
}

} // namespace

AtmosphereTransmission::AtmosphereTransmission(const AtmosphereTransmissionConfig& cfg)
    : cfg_(cfg)
{
    if (!cfg_.enabled) {
        return;
    }
    if (cfg_.model != "modtran_tau") {
        throw std::runtime_error("unsupported atmosphere.model: " + cfg_.model);
    }
    if (cfg_.tau_table_path.empty()) {
        throw std::runtime_error("atmosphere.tau_table is required for modtran_tau");
    }
    if (cfg_.min_cos_theta <= 0.0 || cfg_.min_cos_theta > 1.0) {
        throw std::runtime_error("atmosphere.min_cos_theta must be in (0, 1]");
    }
    loadModtranTauTable(cfg_.tau_table_path);
}

void AtmosphereTransmission::loadModtranTauTable(const std::string& path)
{
    std::ifstream ifs(path);
    if (!ifs) {
        throw std::runtime_error("failed to open atmosphere tau table: " + path);
    }

    std::vector<double> rows;
    std::string line;
    double table_h2_km = std::numeric_limits<double>::quiet_NaN();
    while (std::getline(ifs, line)) {
        const std::string trimmed = trimCopy(line);
        if (trimmed.empty()) {
            continue;
        }
        if (trimmed.rfind("# H2=", 0) == 0) {
            const auto h2_pos = trimmed.find("H2=");
            const auto h2_end = trimmed.find(',', h2_pos == std::string::npos ? 0 : h2_pos);
            if (h2_pos != std::string::npos && h2_end != std::string::npos) {
                const std::string h2_text =
                    trimCopy(trimmed.substr(h2_pos + 3, h2_end - (h2_pos + 3)));
                double h2_value = 0.0;
                if (parseDoubleLocal(h2_text, h2_value)) {
                    table_h2_km = h2_value;
                    altitudes_km_.push_back(table_h2_km);
                }
            }
            const auto h1_pos = trimmed.find("H1=");
            if (h1_pos == std::string::npos) {
                continue;
            }
            std::stringstream ss(trimmed.substr(h1_pos + 3));
            double altitude = 0.0;
            while (ss >> altitude) {
                altitudes_km_.push_back(altitude);
            }
            continue;
        }
        if (trimmed[0] == '#') {
            continue;
        }

        std::stringstream ss(trimmed);
        double wavelength = 0.0;
        if (!(ss >> wavelength)) {
            continue;
        }
        wavelengths_nm_.push_back(wavelength);
        if (std::isfinite(table_h2_km)) {
            rows.push_back(0.0);
        }
        double value = 0.0;
        while (ss >> value) {
            rows.push_back(value);
        }
    }

    if (wavelengths_nm_.empty() || altitudes_km_.empty()) {
        throw std::runtime_error("atmosphere tau table missing wavelength rows or H1 header: " + path);
    }
    const std::size_t expected = wavelengths_nm_.size() * altitudes_km_.size();
    if (rows.size() != expected) {
        std::ostringstream oss;
        oss << "atmosphere tau table has " << rows.size()
            << " tau values, expected " << expected;
        throw std::runtime_error(oss.str());
    }
    tau_ = std::move(rows);
    log_altitudes_cm_.reserve(altitudes_km_.size());
    for (double altitude_km : altitudes_km_) {
        const double log_altitude = altitudeKmToLogCm(altitude_km);
        if (!std::isfinite(log_altitude)) {
            throw std::runtime_error("atmosphere tau table has invalid altitude");
        }
        log_altitudes_cm_.push_back(log_altitude);
    }
}

AtmosphereTransmission::LookupIndex
AtmosphereTransmission::locate(const std::vector<double>& grid, double value) const
{
    LookupIndex idx;
    if (grid.empty() || !std::isfinite(value)) {
        idx.valid = false;
        return idx;
    }
    if (value < grid.front()) {
        if (cfg_.out_of_range == "zero") {
            idx.valid = false;
            return idx;
        }
        idx.lo = idx.hi = 0;
        return idx;
    }
    if (value > grid.back()) {
        if (cfg_.out_of_range == "zero") {
            idx.valid = false;
            return idx;
        }
        idx.lo = idx.hi = grid.size() - 1;
        return idx;
    }

    auto it = std::lower_bound(grid.begin(), grid.end(), value);
    if (it == grid.begin()) {
        idx.lo = idx.hi = 0;
        return idx;
    }
    if (it == grid.end()) {
        idx.lo = idx.hi = grid.size() - 1;
        return idx;
    }
    idx.hi = static_cast<std::size_t>(it - grid.begin());
    if (*it == value) {
        idx.lo = idx.hi = idx.hi;
        return idx;
    }
    idx.lo = idx.hi - 1;
    idx.frac = (value - grid[idx.lo]) / (grid[idx.hi] - grid[idx.lo]);
    return idx;
}

double AtmosphereTransmission::tauAt(std::size_t iw, std::size_t ih) const
{
    const double value = tau_[iw * altitudes_km_.size() + ih];
    if (!std::isfinite(value) || value >= 99999.0) {
        if (cfg_.invalid_tau == "zero") {
            return 0.0;
        }
        return std::numeric_limits<double>::infinity();
    }
    return value;
}

double AtmosphereTransmission::interpolateTau(double wavelength_nm, double altitude_km) const
{
    const LookupIndex wi = locate(wavelengths_nm_, wavelength_nm);
    const LookupIndex hi = locate(log_altitudes_cm_, altitudeKmToLogCm(altitude_km));
    if (!wi.valid || !hi.valid) {
        return std::numeric_limits<double>::infinity();
    }
    const double t00 = tauAt(wi.lo, hi.lo);
    const double t10 = tauAt(wi.hi, hi.lo);
    const double t01 = tauAt(wi.lo, hi.hi);
    const double t11 = tauAt(wi.hi, hi.hi);
    if (!std::isfinite(t00) || !std::isfinite(t10) ||
        !std::isfinite(t01) || !std::isfinite(t11)) {
        return std::numeric_limits<double>::infinity();
    }
    const double t0 = t00 + wi.frac * (t10 - t00);
    const double t1 = t01 + wi.frac * (t11 - t01);
    return t0 + hi.frac * (t1 - t0);
}

double AtmosphereTransmission::resolveEmissionAltitude(double emission_altitude_km) const
{
    if (std::isfinite(emission_altitude_km)) {
        return emission_altitude_km;
    }
    if (cfg_.missing_emission_altitude == "default") {
        return cfg_.default_emission_altitude_km;
    }
    throw std::runtime_error("atmosphere modtran_tau requires photon emission_altitude_km");
}

double AtmosphereTransmission::transmission(double wavelength_nm,
                                            double emission_altitude_km,
                                            const Vec3& global_direction) const
{
    if (!cfg_.enabled) {
        return 1.0;
    }
    const double altitude = resolveEmissionAltitude(emission_altitude_km);
    double tau = interpolateTau(wavelength_nm, altitude);
    if (!std::isfinite(tau)) {
        return 0.0;
    }
    if (cfg_.slant_correction) {
        double cos_theta = -global_direction.normalized().z;
        if (cos_theta <= 0.0) {
            return 0.0;
        }
        if (cos_theta < cfg_.min_cos_theta) {
            if (cfg_.large_angle_behavior == "error") {
                throw std::runtime_error("atmosphere secant correction exceeded min_cos_theta");
            }
            if (cfg_.large_angle_behavior == "opaque" ||
                cfg_.large_angle_behavior == "zero") {
                return 0.0;
            }
            cos_theta = cfg_.min_cos_theta;
        }
        tau /= cos_theta;
    }
    return std::exp(-tau);
}

AtmosphereTransmissionConfig buildAtmosphereTransmissionConfig(
    const std::map<std::string, std::string>& cfg)
{
    AtmosphereTransmissionConfig out;
    out.model = lowerCopyLocal(trimCopy(getStringLocal(cfg, "atmosphere.model", out.model)));
    if (disabledText(out.model)) {
        out.enabled = false;
        out.model = "none";
        return out;
    }
    out.enabled = true;
    out.tau_table_path = getStringLocal(cfg, "atmosphere.tau_table", "");
    out.ground_altitude_km =
        getDoubleLocal(cfg, "atmosphere.ground_altitude_km", out.ground_altitude_km);
    const std::string slant =
        lowerCopyLocal(trimCopy(getStringLocal(cfg, "atmosphere.slant_correction", "secant")));
    out.slant_correction = !(disabledText(slant) || slant == "none");
    out.min_cos_theta = getDoubleLocal(cfg, "atmosphere.min_cos_theta", out.min_cos_theta);
    out.large_angle_behavior =
        lowerCopyLocal(trimCopy(getStringLocal(cfg, "atmosphere.large_angle_behavior",
                                               out.large_angle_behavior)));
    out.out_of_range =
        lowerCopyLocal(trimCopy(getStringLocal(cfg, "atmosphere.out_of_range",
                                               out.out_of_range)));
    out.invalid_tau =
        lowerCopyLocal(trimCopy(getStringLocal(cfg, "atmosphere.invalid_tau", out.invalid_tau)));
    out.missing_emission_altitude =
        lowerCopyLocal(trimCopy(getStringLocal(cfg, "atmosphere.missing_emission_altitude",
                                               out.missing_emission_altitude)));
    out.default_emission_altitude_km =
        getDoubleLocal(cfg, "atmosphere.default_emission_altitude_km",
                       out.default_emission_altitude_km);
    return out;
}

std::string atmosphereTransmissionDescription(const AtmosphereTransmissionConfig& cfg)
{
    if (!cfg.enabled) {
        return "not set";
    }
    std::ostringstream oss;
    oss << cfg.model << ": " << cfg.tau_table_path
        << ", H2=" << cfg.ground_altitude_km << " km"
        << ", slant=" << (cfg.slant_correction ? "secant" : "none")
        << ", min_cos=" << cfg.min_cos_theta;
    return oss.str();
}
