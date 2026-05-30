#include "optics/AtmosphereTransmission.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

bool nearlyEqual(double a, double b, double eps = 1e-6)
{
    return std::abs(a - b) <= eps;
}

bool check(bool condition, const char* message)
{
    if (!condition) {
        std::cerr << "FAIL: " << message << "\n";
        return false;
    }
    return true;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "usage: test_atmosphere_transmission /path/to/modtran.dat\n";
        return 2;
    }

    AtmosphereTransmissionConfig cfg;
    cfg.enabled = true;
    cfg.model = "modtran_tau";
    cfg.tau_table_path = argv[1];
    cfg.slant_correction = false;
    AtmosphereTransmission atm(cfg);

    bool ok = true;
    const double vertical = atm.transmission(440.0, 100.0, {0.0, 0.0, -1.0});
    ok &= check(nearlyEqual(vertical, std::exp(-0.207110), 1e-8),
                "440 nm, 100 km vertical transmission");

    AtmosphereTransmissionConfig slant_cfg = cfg;
    slant_cfg.slant_correction = true;
    AtmosphereTransmission slant_atm(slant_cfg);
    const double slant = slant_atm.transmission(440.0, 100.0,
                                                {std::sqrt(0.75), 0.0, -0.5});
    ok &= check(nearlyEqual(slant, std::exp(-0.207110 / 0.5), 1e-8),
                "secant slant correction");

    const double invalid = atm.transmission(200.0, 100.0, {0.0, 0.0, -1.0});
    ok &= check(invalid == 0.0, "99999 tau sentinel is opaque");

    if (!ok) {
        return 1;
    }
    std::cout << "Atmosphere transmission checks passed\n";
    return 0;
}
