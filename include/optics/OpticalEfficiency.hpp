#pragma once
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

class EfficiencyCurve {
public:
    bool empty() const {
        return points_.empty();
    }

    void loadCsv(const std::string& path) {
        std::ifstream ifs(path);
        if (!ifs) {
            throw std::runtime_error("failed to open efficiency curve: " + path);
        }

        points_.clear();
        std::string line;
        while (std::getline(ifs, line)) {
            auto comment = line.find('#');
            if (comment != std::string::npos) {
                line = line.substr(0, comment);
            }
            if (line.empty()) {
                continue;
            }
            std::replace(line.begin(), line.end(), ',', ' ');
            std::stringstream ss(line);
            double wavelength_nm = 0.0;
            double efficiency = 0.0;
            if (!(ss >> wavelength_nm >> efficiency)) {
                continue;
            }
            points_.push_back({wavelength_nm, efficiency});
        }

        std::sort(points_.begin(), points_.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });

        if (points_.empty()) {
            throw std::runtime_error("efficiency curve has no numeric rows: " + path);
        }
    }

    double evaluate(double wavelength_nm) const {
        if (points_.empty() || wavelength_nm < points_.front().first ||
            wavelength_nm > points_.back().first) {
            return 0.0;
        }

        auto it = std::lower_bound(
            points_.begin(), points_.end(), wavelength_nm,
            [](const auto& p, double wl) { return p.first < wl; });

        if (it == points_.begin()) {
            return it->second;
        }
        if (it == points_.end()) {
            return points_.back().second;
        }
        if (it->first == wavelength_nm) {
            return it->second;
        }

        const auto& hi = *it;
        const auto& lo = *(it - 1);
        double frac = (wavelength_nm - lo.first) / (hi.first - lo.first);
        return lo.second + frac * (hi.second - lo.second);
    }

private:
    std::vector<std::pair<double, double>> points_;
};

struct EfficiencyFactorConfig {
    bool enabled = false;
    bool use_curve = false;
    double constant = 1.0;
    std::string csv_path;
};

struct OpticalEfficiencyConfig {
    bool use_funnel_acceptance = false;

    double constant_scale = 1.0;

    EfficiencyFactorConfig mirror_reflectivity;
    EfficiencyFactorConfig filter_transmission;
    EfficiencyFactorConfig sipm_pde;
    EfficiencyFactorConfig atmosphere_transmission;
};

class OpticalEfficiency {
public:
    OpticalEfficiency() = default;

    explicit OpticalEfficiency(const OpticalEfficiencyConfig& cfg)
        : cfg_(cfg)
    {
        if (cfg_.mirror_reflectivity.enabled && cfg_.mirror_reflectivity.use_curve) {
            mirror_reflectivity_.loadCsv(cfg_.mirror_reflectivity.csv_path);
        }
        if (cfg_.filter_transmission.enabled && cfg_.filter_transmission.use_curve) {
            filter_transmission_.loadCsv(cfg_.filter_transmission.csv_path);
        }
        if (cfg_.sipm_pde.enabled && cfg_.sipm_pde.use_curve) {
            sipm_pde_.loadCsv(cfg_.sipm_pde.csv_path);
        }
        if (cfg_.atmosphere_transmission.enabled && cfg_.atmosphere_transmission.use_curve) {
            atmosphere_transmission_.loadCsv(cfg_.atmosphere_transmission.csv_path);
        }
    }

    double mirrorReflectivity(double wavelength_nm) const {
        return evaluateFactor(cfg_.mirror_reflectivity, mirror_reflectivity_, wavelength_nm);
    }

    double telescopeTransmission(double wavelength_nm) const {
        return evaluateFactor(cfg_.filter_transmission, filter_transmission_, wavelength_nm);
    }

    double funnelAcceptance(double incidence_angle_rad, double /*wavelength_nm*/) const {
        if (!cfg_.use_funnel_acceptance) return 1.0;
        double c = std::cos(incidence_angle_rad);
        return std::max(0.0, c);
    }

    double cathodeEfficiency(double wavelength_nm) const {
        return evaluateFactor(cfg_.sipm_pde, sipm_pde_, wavelength_nm);
    }

    double atmosphereTransmission(double wavelength_nm) const {
        return evaluateFactor(cfg_.atmosphere_transmission,
                              atmosphere_transmission_,
                              wavelength_nm);
    }

    double total(double wavelength_nm, double incidence_angle_rad) const {
        return cfg_.constant_scale
             * mirrorReflectivity(wavelength_nm)
             * telescopeTransmission(wavelength_nm)
             * atmosphereTransmission(wavelength_nm)
             * funnelAcceptance(incidence_angle_rad, wavelength_nm)
             * cathodeEfficiency(wavelength_nm);
    }

private:
    static double evaluateFactor(const EfficiencyFactorConfig& cfg,
                                 const EfficiencyCurve& curve,
                                 double wavelength_nm)
    {
        if (!cfg.enabled) return 1.0;
        if (cfg.use_curve) return curve.evaluate(wavelength_nm);
        return cfg.constant;
    }

    OpticalEfficiencyConfig cfg_;
    EfficiencyCurve mirror_reflectivity_;
    EfficiencyCurve filter_transmission_;
    EfficiencyCurve sipm_pde_;
    EfficiencyCurve atmosphere_transmission_;
};
