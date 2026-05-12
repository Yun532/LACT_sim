#include "optics/OpticalEfficiency.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace {

bool check(bool cond, const std::string& msg)
{
    if (!cond) {
        std::cerr << msg << "\n";
        return false;
    }
    return cond;
}

bool near(double a, double b, double tol = 1e-12)
{
    return std::abs(a - b) < tol;
}

std::filesystem::path writeCurve(const std::string& name, const std::string& text)
{
    auto path = std::filesystem::temp_directory_path() / name;
    std::ofstream os(path);
    os << text;
    return path;
}

EfficiencyFactorConfig curveFactor(const std::filesystem::path& path)
{
    EfficiencyFactorConfig cfg;
    cfg.enabled = true;
    cfg.use_curve = true;
    cfg.csv_path = path.string();
    return cfg;
}

EfficiencyFactorConfig constantFactor(double value)
{
    EfficiencyFactorConfig cfg;
    cfg.enabled = true;
    cfg.use_curve = false;
    cfg.constant = value;
    return cfg;
}

} // namespace

int main()
{
    bool ok = true;

    const auto curve_path = writeCurve(
        "lact_efficiency_curve_test.csv",
        "wavelength_nm,efficiency\n"
        "300,0.10\n"
        "400,0.50\n"
        "500,0.90\n");

    EfficiencyCurve curve;
    curve.loadCsv(curve_path.string());
    ok &= check(near(curve.evaluate(300.0), 0.10), "exact first point");
    ok &= check(near(curve.evaluate(400.0), 0.50), "exact middle point");
    ok &= check(near(curve.evaluate(350.0), 0.30), "linear interpolation");
    ok &= check(near(curve.evaluate(250.0), 0.0), "below range should be zero");
    ok &= check(near(curve.evaluate(550.0), 0.0), "above range should be zero");

    const auto duplicate_path = writeCurve(
        "lact_efficiency_duplicate_test.csv",
        "wavelength_nm,efficiency\n"
        "400,0.20\n"
        "400,0.40\n"
        "500,0.80\n");
    EfficiencyCurve duplicate_curve;
    duplicate_curve.loadCsv(duplicate_path.string());
    ok &= check(near(duplicate_curve.evaluate(400.0), 0.30),
                "duplicate wavelengths should be averaged");
    ok &= check(near(duplicate_curve.evaluate(450.0), 0.55),
                "interpolation after duplicate merge");

    OpticalEfficiencyConfig cfg;
    cfg.constant_scale = 2.0;
    cfg.mirror_reflectivity = curveFactor(curve_path);
    cfg.filter_transmission = constantFactor(0.80);
    cfg.sipm_pde = constantFactor(0.25);
    cfg.atmosphere_transmission = constantFactor(0.90);
    OpticalEfficiency eff(cfg);
    ok &= check(near(eff.total(400.0, 0.0), 2.0 * 0.50 * 0.80 * 0.25 * 0.90),
                "total efficiency product");

    OpticalEfficiency ideal;
    ok &= check(near(ideal.total(400.0, 0.0), 1.0),
                "disabled efficiency factors should default to one");

    std::cout << "efficiency curve checks passed\n";
    return ok ? 0 : 1;
}
