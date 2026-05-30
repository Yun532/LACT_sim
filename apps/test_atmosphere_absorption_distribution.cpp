#include "optics/AtmosphereTransmission.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

struct Row {
    double altitude_km;
    double before_weight;
    double after_weight;
    double theory_weight;
};

bool nearlyEqual(double a, double b, double eps)
{
    return std::abs(a - b) <= eps;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr << "usage: test_atmosphere_absorption_distribution "
                  << "/path/to/modtran.dat /path/to/output.csv\n";
        return 2;
    }

    AtmosphereTransmissionConfig cfg;
    cfg.enabled = true;
    cfg.model = "modtran_tau";
    cfg.tau_table_path = argv[1];
    cfg.slant_correction = false;
    AtmosphereTransmission atm(cfg);

    const double wavelength_nm = 440.0;
    const double before = 10000.0;
    const std::vector<double> heights_km{
        4.45, 4.6, 5.0, 6.5, 10.0, 20.0, 40.0, 100.0};

    std::vector<Row> rows;
    rows.reserve(heights_km.size());
    bool ok = true;
    for (double h : heights_km) {
        const double t = atm.transmission(wavelength_nm, h, {0.0, 0.0, -1.0});
        const double after = before * t;
        const double theory = before * t;
        rows.push_back({h, before, after, theory});
        ok = ok && nearlyEqual(after, theory, 1e-9);
    }

    const std::filesystem::path out_path(argv[2]);
    if (out_path.has_parent_path()) {
        std::filesystem::create_directories(out_path.parent_path());
    }
    std::ofstream out(out_path);
    if (!out) {
        throw std::runtime_error("failed to write " + out_path.string());
    }
    out << "wavelength_nm,altitude_km,before_weight,after_weight,"
        << "theory_weight,transmission\n";
    for (const Row& row : rows) {
        out << wavelength_nm << ','
            << row.altitude_km << ','
            << row.before_weight << ','
            << row.after_weight << ','
            << row.theory_weight << ','
            << row.after_weight / row.before_weight << '\n';
    }

    if (!ok) {
        std::cerr << "atmosphere absorption distribution mismatch\n";
        return 1;
    }
    std::cout << "Wrote atmosphere absorption distribution to "
              << out_path.string() << "\n";
    return 0;
}
