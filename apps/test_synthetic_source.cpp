#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "io/PhotonCsvSource.hpp"
#include "io/SyntheticPhotonSource.hpp"

namespace {

bool nearlyEqual(double a, double b, double tol = 1e-12) {
    return std::abs(a - b) <= tol;
}

bool check(bool condition, const std::string& label) {
    if (!condition) {
        std::cerr << "FAILED: " << label << "\n";
        return false;
    }
    return true;
}

bool expectInvalidConfig(const SyntheticPhotonConfig& cfg, const std::string& label) {
    try {
        SyntheticPhotonSource source(cfg);
    } catch (const std::runtime_error&) {
        return true;
    }
    std::cerr << "FAILED: " << label << "\n";
    return false;
}

bool checkParallelBeam() {
    SyntheticPhotonConfig cfg;
    cfg.mode = SyntheticMode::ParallelBeam;
    cfg.n_bunches = 8;
    cfg.multiplicity = 3.0;
    cfg.photon_weight = 0.5;
    cfg.wavelength_nm = 405.0;
    cfg.time_ns = 12.5;
    cfg.source_plane_z = 120.0;
    cfg.beam_radius_m = 0.5;
    cfg.beam_direction = {0.01, 0.0, -1.0};
    cfg.event_id = 7;
    cfg.telescope_id = 2;
    cfg.random_seed = 99;

    SyntheticPhotonSource source(cfg);
    PhotonBunch first;
    bool ok = check(source.next(first), "parallel source produces first bunch");

    Vec3 expected_dir = cfg.beam_direction.normalized();
    ok &= check(nearlyEqual(first.photon.dir.x, expected_dir.x), "parallel direction x");
    ok &= check(nearlyEqual(first.photon.dir.y, expected_dir.y), "parallel direction y");
    ok &= check(nearlyEqual(first.photon.dir.z, expected_dir.z), "parallel direction z");
    ok &= check(nearlyEqual(first.photon.pos.z, cfg.source_plane_z), "parallel source plane z");
    ok &= check(first.photon.pos.x * first.photon.pos.x +
                first.photon.pos.y * first.photon.pos.y <=
                cfg.beam_radius_m * cfg.beam_radius_m + 1e-12,
                "parallel source samples inside disk");
    ok &= check(nearlyEqual(first.photon.wavelength_nm, cfg.wavelength_nm), "parallel wavelength");
    ok &= check(nearlyEqual(first.photon.time_ns, cfg.time_ns), "parallel time");
    ok &= check(nearlyEqual(first.photon.weight, cfg.photon_weight), "parallel photon weight");
    ok &= check(nearlyEqual(first.multiplicity, cfg.multiplicity), "parallel multiplicity");
    ok &= check(first.event_id == cfg.event_id, "parallel event id");
    ok &= check(first.telescope_id == cfg.telescope_id, "parallel telescope id");
    ok &= check(first.source_bunch_index == 0,
                "parallel first random-stream index");

    int produced = 1;
    PhotonBunch bunch;
    while (source.next(bunch)) {
        ++produced;
        ok &= check(bunch.source_bunch_index ==
                        static_cast<std::uint64_t>(produced - 1),
                    "parallel random-stream indices are sequential");
        ok &= check(nearlyEqual(bunch.photon.pos.z, cfg.source_plane_z), "parallel all z");
        ok &= check(nearlyEqual(bunch.photon.dir.norm(), 1.0), "parallel all directions unit");
    }
    ok &= check(produced == cfg.n_bunches, "parallel source count");
    ok &= check(!source.next(bunch), "parallel source exhausted");

    source.reset();
    PhotonBunch replay;
    ok &= check(source.next(replay), "parallel reset produces bunch");
    ok &= check(replay.source_bunch_index == 0,
                "parallel reset restarts random-stream index");
    ok &= check(nearlyEqual(replay.photon.pos.x, first.photon.pos.x), "parallel reset x");
    ok &= check(nearlyEqual(replay.photon.pos.y, first.photon.pos.y), "parallel reset y");
    ok &= check(nearlyEqual(replay.photon.pos.z, first.photon.pos.z), "parallel reset z");

    return ok;
}

bool checkPointSource() {
    SyntheticPhotonConfig cfg;
    cfg.mode = SyntheticMode::PointSource;
    cfg.n_bunches = 12;
    cfg.multiplicity = 2.0;
    cfg.source_position = {1.0, -2.0, 50.0};
    cfg.aperture_z = 0.5;
    cfg.aperture_radius_m = 1.25;
    cfg.random_seed = 123;

    SyntheticPhotonSource source(cfg);

    bool ok = true;
    int produced = 0;
    PhotonBunch bunch;
    while (source.next(bunch)) {
        ++produced;
        ok &= check(bunch.source_bunch_index ==
                        static_cast<std::uint64_t>(produced - 1),
                    "point-source random-stream indices are sequential");
        ok &= check(nearlyEqual(bunch.photon.pos.x, cfg.source_position.x), "point source pos x");
        ok &= check(nearlyEqual(bunch.photon.pos.y, cfg.source_position.y), "point source pos y");
        ok &= check(nearlyEqual(bunch.photon.pos.z, cfg.source_position.z), "point source pos z");
        ok &= check(nearlyEqual(bunch.photon.dir.norm(), 1.0), "point source direction unit");

        double t = (cfg.aperture_z - bunch.photon.pos.z) / bunch.photon.dir.z;
        Vec3 target = bunch.photon.pos + bunch.photon.dir * t;
        double r2 = target.x * target.x + target.y * target.y;

        ok &= check(t > 0.0, "point source target is in front");
        ok &= check(nearlyEqual(target.z, cfg.aperture_z), "point source intersects aperture z");
        ok &= check(r2 <= cfg.aperture_radius_m * cfg.aperture_radius_m + 1e-12,
                    "point source target inside aperture");
    }

    ok &= check(produced == cfg.n_bunches, "point source count");
    return ok;
}

bool checkInvalidConfigs() {
    bool ok = true;

    SyntheticPhotonConfig cfg;
    cfg.n_bunches = -1;
    ok &= expectInvalidConfig(cfg, "negative n_bunches rejected");

    cfg = SyntheticPhotonConfig{};
    cfg.beam_direction = {0.0, 0.0, 0.0};
    ok &= expectInvalidConfig(cfg, "zero beam direction rejected");

    cfg = SyntheticPhotonConfig{};
    cfg.beam_radius_m = -1.0;
    ok &= expectInvalidConfig(cfg, "negative beam radius rejected");

    cfg = SyntheticPhotonConfig{};
    cfg.mode = SyntheticMode::PointSource;
    cfg.aperture_radius_m = -0.1;
    ok &= expectInvalidConfig(cfg, "negative aperture radius rejected");

    cfg = SyntheticPhotonConfig{};
    cfg.photon_weight = -1.0;
    ok &= expectInvalidConfig(cfg, "negative photon weight rejected");

    cfg = SyntheticPhotonConfig{};
    cfg.wavelength_nm = 0.0;
    ok &= expectInvalidConfig(cfg, "non-positive wavelength rejected");

    return ok;
}

bool checkPhotonCsvRowIndex() {
    const auto path =
        std::filesystem::temp_directory_path() / "lact_photon_csv_row_index.csv";
    {
        std::ofstream output(path);
        output << "x_m,y_m,z_m,dir_x,dir_y,dir_z\n"
               << "0,0,1,0,0,-1\n"
               << "1,0,1,0,0,-1\n";
    }

    PhotonCsvConfig cfg;
    cfg.csv_path = path.string();
    cfg.default_eventio_2d = true;
    PhotonCsvSource source(cfg);
    PhotonBunch first;
    PhotonBunch second;
    const bool have_first = source.next(first);
    const bool have_second = source.next(second);
    std::filesystem::remove(path);

    bool ok = check(have_first && have_second, "PhotonCsv produces two rows");
    ok &= check(first.source_bunch_index == 0,
                "PhotonCsv first implicit random-stream index");
    ok &= check(second.source_bunch_index == 1,
                "PhotonCsv second implicit random-stream index");
    ok &= check(first.eventio_2d && second.eventio_2d,
                "PhotonCsv applies configured default 2D line semantics");
    return ok;
}

bool checkInvalidPhotonCsvRows() {
    bool ok = true;
    const auto path =
        std::filesystem::temp_directory_path() / "lact_photon_csv_invalid.csv";
    const auto expectInvalidRow =
        [&path, &ok](const std::string& row, const char* label) {
            {
                std::ofstream output(path);
                output << "x_m,y_m,z_m,dir_x,dir_y,dir_z,"
                          "wavelength_nm,weight,multiplicity\n"
                       << row << "\n";
            }
            PhotonCsvConfig cfg;
            cfg.csv_path = path.string();
            try {
                PhotonCsvSource source(cfg);
                (void)source;
                std::cerr << label << "\n";
                ok = false;
            } catch (...) {
            }
        };

    expectInvalidRow("0,0,1,0,0,0,400,1,1",
                     "PhotonCsv accepted a zero direction");
    expectInvalidRow("0,0,1,0,0,-1,400,-1,1",
                     "PhotonCsv accepted a negative weight");
    expectInvalidRow("0,0,1,0,0,-1,400,1,-1",
                     "PhotonCsv accepted a negative multiplicity");
    expectInvalidRow("nan,0,1,0,0,-1,400,1,1",
                     "PhotonCsv accepted a non-finite position");
    std::filesystem::remove(path);
    return ok;
}

} // namespace

int main() {
    bool ok = true;
    ok &= checkParallelBeam();
    ok &= checkPointSource();
    ok &= checkInvalidConfigs();
    ok &= checkPhotonCsvRowIndex();
    ok &= checkInvalidPhotonCsvRows();

    if (ok) {
        std::cout << "Synthetic photon source checks passed\n";
        return 0;
    }

    return 1;
}
