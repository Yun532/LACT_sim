#include "io/EventIOPhotonSource.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

std::string lowerCopy(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

std::uint64_t mixSeed(std::uint64_t seed, std::uint64_t value)
{
    seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
    return seed;
}

double unitRandom01(std::uint64_t seed)
{
    // SplitMix64 gives a deterministic pseudo-random value from bunch identity,
    // independent of streaming order and filtering.
    seed += 0x9e3779b97f4a7c15ULL;
    seed = (seed ^ (seed >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    seed = (seed ^ (seed >> 27U)) * 0x94d049bb133111ebULL;
    seed ^= seed >> 31U;
    constexpr double scale = 1.0 / 9007199254740992.0; // 2^53
    return static_cast<double>(seed >> 11U) * scale;
}

double sampleMissingWavelength(int shower_event_id,
                               int array_id,
                               int telescope_id,
                               std::uint64_t source_bunch_index,
                               std::uint64_t photon_index,
                               const EventIOPhotonConfig& cfg)
{
    const std::string model = lowerCopy(cfg.missing_wavelength_model);
    if (model.empty() || model == "default" || model == "fixed" ||
        model == "constant" || model == "none" || model == "off") {
        return cfg.default_wavelength_nm;
    }

    const double lo = cfg.missing_wavelength_min_nm;
    const double hi = cfg.missing_wavelength_max_nm;
    if (!std::isfinite(lo) || !std::isfinite(hi) || lo <= 0.0 || hi <= lo) {
        throw std::runtime_error(
            "source.missing_wavelength_min_nm/max_nm must be finite with 0 < min < max");
    }

    std::uint64_t seed = cfg.missing_wavelength_seed;
    seed = mixSeed(seed, static_cast<std::uint64_t>(shower_event_id + 1000003));
    seed = mixSeed(seed, static_cast<std::uint64_t>(array_id + 100003));
    seed = mixSeed(seed, static_cast<std::uint64_t>(telescope_id + 10007));
    seed = mixSeed(seed, source_bunch_index + 1ULL);
    seed = mixSeed(seed, photon_index + 1ULL);
    const double u = std::clamp(unitRandom01(seed), 1e-16, 1.0 - 1e-16);

    if (model == "uniform") {
        return lo + u * (hi - lo);
    }
    if (model == "cherenkov" || model == "cherenkov_1_over_lambda2" ||
        model == "1_over_lambda2" || model == "one_over_lambda2") {
        const double inv_lo = 1.0 / lo;
        const double inv_hi = 1.0 / hi;
        return 1.0 / (inv_lo - u * (inv_lo - inv_hi));
    }

    throw std::runtime_error("unsupported source.missing_wavelength_model: " +
                             cfg.missing_wavelength_model);
}

} // namespace

double resolveEventIOPhotonWavelength(const PhotonBunch& bunch,
                                      std::uint64_t photon_index,
                                      const EventIOPhotonConfig& cfg)
{
    if (bunch.raw_wavelength_nm > 0.0) {
        return bunch.raw_wavelength_nm;
    }
    if (bunch.raw_wavelength_nm < 0.0) {
        // CEFFIC bunches no longer carry a physical wavelength. The
        // placeholder is used only by wavelength-agnostic geometry code.
        return cfg.default_wavelength_nm;
    }
    return sampleMissingWavelength(bunch.shower_event_id,
                                   bunch.array_id,
                                   bunch.telescope_id,
                                   bunch.source_bunch_index,
                                   photon_index,
                                   cfg);
}

double eventIO3DEmissionAltitudeKm(double observation_altitude_km,
                                   double bunch_z_cm,
                                   double direction_cosine_z,
                                   double emission_distance_cm)
{
    if (!std::isfinite(observation_altitude_km) ||
        !std::isfinite(bunch_z_cm) ||
        !std::isfinite(direction_cosine_z) ||
        !std::isfinite(emission_distance_cm)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return observation_altitude_km +
           (bunch_z_cm - direction_cosine_z * emission_distance_cm) * 1.0e-5;
}
