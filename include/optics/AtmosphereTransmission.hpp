#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include "core/Vec3.hpp"

struct AtmosphereTransmissionConfig {
    bool enabled = false;
    std::string model = "none";
    std::string tau_table_path;
    double ground_altitude_km = 4.4;
    bool slant_correction = true;
    double min_cos_theta = 0.2;
    std::string large_angle_behavior = "clamp";
    std::string out_of_range = "clamp";
    std::string invalid_tau = "opaque";
    std::string missing_emission_altitude = "error";
    double default_emission_altitude_km = 10.0;
};

class AtmosphereTransmission {
public:
    AtmosphereTransmission() = default;
    explicit AtmosphereTransmission(const AtmosphereTransmissionConfig& cfg);

    bool enabled() const { return cfg_.enabled; }
    const AtmosphereTransmissionConfig& config() const { return cfg_; }

    double transmission(double wavelength_nm,
                        double emission_altitude_km,
                        const Vec3& global_direction) const;

private:
    struct LookupIndex {
        std::size_t lo = 0;
        std::size_t hi = 0;
        double frac = 0.0;
        bool valid = true;
    };

    AtmosphereTransmissionConfig cfg_;
    std::vector<double> wavelengths_nm_;
    std::vector<double> altitudes_km_;
    std::vector<double> log_altitudes_cm_;
    std::vector<double> tau_;

    void loadModtranTauTable(const std::string& path);
    LookupIndex locate(const std::vector<double>& grid, double value) const;
    double tauAt(std::size_t iw, std::size_t ih) const;
    double interpolateTau(double wavelength_nm, double altitude_km) const;
    double resolveEmissionAltitude(double emission_altitude_km) const;
};

AtmosphereTransmissionConfig buildAtmosphereTransmissionConfig(
    const std::map<std::string, std::string>& cfg);
std::string atmosphereTransmissionDescription(const AtmosphereTransmissionConfig& cfg);
