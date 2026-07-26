#include "app/PhotonResponseSampler.hpp"

#include <cmath>
#include <iostream>
#include <set>

using namespace lact;

namespace {

bool check(bool condition, const char* message)
{
    if (!condition) {
        std::cerr << message << '\n';
    }
    return condition;
}

PhotonBunch makeBunch(double multiplicity)
{
    PhotonBunch bunch;
    bunch.multiplicity = multiplicity;
    bunch.event_id = 101;
    bunch.shower_event_id = 1;
    bunch.array_id = 1;
    bunch.telescope_id = 4;
    bunch.source_bunch_index = 17;
    bunch.raw_wavelength_nm = 0.0;
    bunch.photon.weight = 1.0;
    return bunch;
}

} // namespace

int main()
{
    bool ok = true;
    EventIOPhotonConfig eventio;
    eventio.missing_wavelength_min_nm = 260.0;
    eventio.missing_wavelength_max_nm = 1000.0;
    eventio.missing_wavelength_seed = 246813579ULL;

    PhotonResponseConfig expectation_cfg;
    expectation_cfg.mode = PhotonResponseMode::Expectation;
    PhotonResponseSampler expectation(expectation_cfg, eventio);
    auto bunch = makeBunch(3.25);
    ok &= check(expectation.candidateCount(bunch) == 1,
                "expectation mode should keep one weighted ray");
    auto expected = expectation.candidate(bunch, 0);
    ok &= check(std::abs(expected.photon.weight - 3.25) < 1.0e-12,
                "expectation mode lost bunch multiplicity");
    ok &= check(expected.stream_id != 0 &&
                    expected.photon.random_stream_id != 0,
                "expectation mode should assign non-zero response and photon streams");
    const auto expected_repeat = expectation.candidate(bunch, 0);
    ok &= check(expected_repeat.stream_id == expected.stream_id &&
                    expected_repeat.photon.random_stream_id ==
                        expected.photon.random_stream_id,
                "expectation mode streams should be reproducible");
    auto translated_bunch = bunch;
    translated_bunch.photon.pos = {123.0, -456.0, 789.0};
    translated_bunch.photon.dir = {0.1, 0.2, -0.97};
    const auto translated_expected = expectation.candidate(translated_bunch, 0);
    ok &= check(translated_expected.stream_id == expected.stream_id,
                "photon identity stream must not depend on coordinates");
    auto expected_pre = expectation.applyPreGeometry(expected, 0.4, 0.2);
    ok &= check(expected_pre.survives &&
                std::abs(expected.photon.weight - 1.3) < 1.0e-12,
                "expectation atmosphere weight is incorrect");

    PhotonResponseConfig stochastic_cfg;
    stochastic_cfg.mode = PhotonResponseMode::StochasticPe;
    stochastic_cfg.seed = 42;
    PhotonResponseSampler stochastic(stochastic_cfg, eventio);
    ok &= check(stochastic.candidateCount(bunch) == 4,
                "fractional stochastic bunch should produce a remainder candidate");
    std::set<double> wavelengths;
    std::set<std::uint64_t> streams;
    for (std::uint64_t i = 0; i < 4; ++i) {
        const auto candidate = stochastic.candidate(bunch, i);
        wavelengths.insert(candidate.photon.wavelength_nm);
        streams.insert(candidate.stream_id);
        const double wanted_fraction = i == 3 ? 0.25 : 1.0;
        ok &= check(std::abs(candidate.represented_fraction - wanted_fraction) < 1.0e-12,
                    "represented photon fraction is incorrect");
    }
    ok &= check(wavelengths.size() == 4,
                "missing wavelengths were not sampled independently");
    ok &= check(streams.size() == 4,
                "represented photons reused a random stream");
    ok &= check(stochastic.candidate(bunch, 2).stream_id ==
                    stochastic.candidate(bunch, 2).stream_id,
                "stochastic stream is not reproducible");

    auto certain = stochastic.candidate(bunch, 0);
    auto certain_pre = stochastic.applyPreGeometry(certain, 1.0, 1.0);
    ok &= check(certain_pre.survives &&
                certain.photon.optical_efficiency_preapplied,
                "unit pre-geometry probability should always survive");
    OpticalSurfaceHit hit;
    hit.relative_efficiency = 1.0;
    ok &= check(stochastic.acceptPostGeometry(certain, hit) &&
                hit.relative_efficiency == 1.0,
                "unit post-geometry probability should always survive");

    auto impossible = stochastic.candidate(bunch, 0);
    auto impossible_pre = stochastic.applyPreGeometry(impossible, 0.0, 1.0);
    ok &= check(!impossible_pre.survives,
                "zero pre-geometry probability should never survive");
    try {
        auto invalid_probability = stochastic.candidate(bunch, 0);
        OpticalSurfaceHit invalid_hit;
        invalid_hit.relative_efficiency = 1.01;
        (void)stochastic.acceptPostGeometry(
            invalid_probability, invalid_hit);
        ok = false;
        std::cerr << "post-geometry probability above one was silently clamped\n";
    } catch (...) {
    }

    std::uint64_t pre_survivors = 0;
    std::uint64_t detected = 0;
    constexpr std::uint64_t trials = 20000;
    for (std::uint64_t i = 0; i < trials; ++i) {
        auto trial_bunch = makeBunch(1.0);
        trial_bunch.source_bunch_index = i;
        auto trial = stochastic.candidate(trial_bunch, 0);
        const auto pre = stochastic.applyPreGeometry(trial, 0.5, 0.4);
        if (!pre.survives) {
            continue;
        }
        ++pre_survivors;
        OpticalSurfaceHit trial_hit;
        trial_hit.relative_efficiency = 0.25;
        if (stochastic.acceptPostGeometry(trial, trial_hit)) {
            ++detected;
        }
    }
    const double pre_fraction = static_cast<double>(pre_survivors) / trials;
    const double detected_fraction = static_cast<double>(detected) / trials;
    ok &= check(std::abs(pre_fraction - 0.20) < 0.015,
                "pre-geometry Bernoulli sampling has the wrong mean");
    ok &= check(std::abs(detected_fraction - 0.05) < 0.008,
                "combined stochastic detector response has the wrong mean");

    auto ceffic_bunch = makeBunch(0.25);
    ceffic_bunch.raw_wavelength_nm = -1.0;
    ceffic_bunch.photon.optical_efficiency_preapplied = true;
    auto ceffic = stochastic.candidate(ceffic_bunch, 0);
    auto ceffic_pre = stochastic.applyPreGeometry(ceffic, 0.0, 0.0);
    ok &= check(ceffic_pre.survives &&
                std::abs(ceffic.remaining_probability - 0.25) < 1.0e-12,
                "CEFFIC bunch should bypass wavelength response but retain remainder probability");

    const auto default_config = buildPhotonResponseConfig({});
    ok &= check(default_config.stochastic(),
                "unspecified response mode should default to stochastic_pe");
    const auto explicit_expectation = buildPhotonResponseConfig({
        {"response.mode", "expectation"}});
    ok &= check(!explicit_expectation.stochastic(),
                "explicit expectation response mode was not preserved");
    const auto parsed = buildPhotonResponseConfig({
        {"response.mode", "stochastic"}, {"response.seed", "99"}});
    ok &= check(parsed.stochastic() && parsed.seed == 99,
                "response configuration parsing failed");
    try {
        (void)buildPhotonResponseConfig({{"response.mode", "invalid"}});
        ok = false;
        std::cerr << "invalid response mode was accepted\n";
    } catch (...) {
    }
    try {
        auto invalid_bunch = makeBunch(-1.0);
        (void)stochastic.candidateCount(invalid_bunch);
        ok = false;
        std::cerr << "negative stochastic multiplicity was accepted\n";
    } catch (...) {
    }
    try {
        auto invalid_bunch = makeBunch(-1.0);
        (void)expectation.candidateCount(invalid_bunch);
        ok = false;
        std::cerr << "negative expectation multiplicity was accepted\n";
    } catch (...) {
    }

    return ok ? 0 : 1;
}
