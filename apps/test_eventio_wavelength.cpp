#include "io/EventIOPhotonSource.hpp"

#include <cmath>
#include <cstdint>
#include <iostream>

int main(int argc, char** argv)
{
    EventIOPhotonConfig cfg;
    cfg.missing_wavelength_model = "cherenkov";
    cfg.missing_wavelength_min_nm = 260.0;
    cfg.missing_wavelength_max_nm = 1000.0;
    cfg.missing_wavelength_seed = 17;

    PhotonBunch bunch;
    bunch.shower_event_id = 3;
    bunch.array_id = 4;
    bunch.telescope_id = 5;
    bunch.source_bunch_index = 6;
    bunch.raw_wavelength_nm = 0.0;

    const double first = resolveEventIOPhotonWavelength(bunch, 0, cfg);
    const double repeat = resolveEventIOPhotonWavelength(bunch, 0, cfg);
    const double second = resolveEventIOPhotonWavelength(bunch, 1, cfg);
    if (first < 260.0 || first > 1000.0 ||
        second < 260.0 || second > 1000.0) {
        std::cerr << "sampled wavelength outside configured range\n";
        return 1;
    }
    if (first != repeat) {
        std::cerr << "wavelength sampling is not reproducible\n";
        return 1;
    }
    if (first == second) {
        std::cerr << "represented photons reused one wavelength draw\n";
        return 1;
    }

    bunch.raw_wavelength_nm = 415.0;
    if (resolveEventIOPhotonWavelength(bunch, 10, cfg) != 415.0) {
        std::cerr << "positive EventIO wavelength changed\n";
        return 1;
    }

    bunch.raw_wavelength_nm = -1.0;
    if (resolveEventIOPhotonWavelength(bunch, 10, cfg) !=
        cfg.default_wavelength_nm) {
        std::cerr << "negative CEFFIC wavelength placeholder changed\n";
        return 1;
    }

    const double emission_altitude =
        eventIO3DEmissionAltitudeKm(4.4, 0.0, -1.0, 560000.0);
    if (std::abs(emission_altitude - 10.0) > 1.0e-12) {
        std::cerr << "3D EventIO emission altitude conversion is wrong\n";
        return 1;
    }

    if (argc == 2) {
        cfg.path = argv[1];
        cfg.event_id_mode = "event_array100";
        cfg.max_shower_events = 1;
        const auto metadata = readEventIOMetadata(cfg);
        if (!std::isfinite(metadata.observation_altitude_m)) {
            std::cerr << "real EventIO input has no observation altitude\n";
            return 1;
        }
        cfg.observation_altitude_km = metadata.observation_altitude_m * 1.0e-3;

        std::uint64_t positive = 0;
        std::uint64_t zero = 0;
        std::uint64_t negative = 0;
        std::uint64_t bunch_2d = 0;
        std::uint64_t bunch_3d = 0;
        std::uint64_t missing_3d_altitude = 0;
        const auto stats = streamEventIOPhotonBunches(
            cfg,
            [&](const PhotonBunch& item) {
                if (item.raw_wavelength_nm > 0.0) ++positive;
                else if (item.raw_wavelength_nm < 0.0) ++negative;
                else ++zero;
                if (item.eventio_2d) {
                    ++bunch_2d;
                } else {
                    ++bunch_3d;
                    if (!std::isfinite(item.emission_altitude_km)) {
                        ++missing_3d_altitude;
                    }
                }
            });
        if (stats.photon_bunches != positive + zero + negative ||
            stats.photon_bunches != bunch_2d + bunch_3d) {
            std::cerr << "real EventIO bunch accounting mismatch\n";
            return 1;
        }
        if (bunch_3d > 0 && missing_3d_altitude > 0) {
            std::cerr << "3D EventIO bunches are missing emission altitude\n";
            return 1;
        }
        std::cout << "observation_altitude_m=" << metadata.observation_altitude_m
                  << " bunches=" << stats.photon_bunches
                  << " wavelength_positive=" << positive
                  << " wavelength_zero=" << zero
                  << " wavelength_negative=" << negative
                  << " bunch_2d=" << bunch_2d
                  << " bunch_3d=" << bunch_3d
                  << " missing_3d_altitude=" << missing_3d_altitude << '\n';
    }

    return 0;
}
