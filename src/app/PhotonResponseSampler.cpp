#include "app/PhotonResponseSampler.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace lact {
namespace {

constexpr std::uint64_t kPreGeometryStage = 1;
constexpr std::uint64_t kPostGeometryStage = 2;

std::uint64_t mixSeed(std::uint64_t seed, std::uint64_t value)
{
    seed ^= value + 0x9e3779b97f4a7c15ULL +
            (seed << 6U) + (seed >> 2U);
    return seed;
}

std::string lowerCopy(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

std::string configValue(const std::map<std::string, std::string>& cfg,
                        const std::string& key,
                        const std::string& fallback)
{
    const auto found = cfg.find(key);
    return found == cfg.end() ? fallback : found->second;
}

std::uint64_t parseSeed(const std::string& text)
{
    if (text.empty() || text.front() == '-') {
        throw std::runtime_error("response.seed must be an unsigned integer");
    }
    std::size_t consumed = 0;
    const auto value = std::stoull(text, &consumed);
    if (consumed != text.size()) {
        throw std::runtime_error("response.seed must be an unsigned integer");
    }
    return value;
}

double probability(double value)
{
    if (!std::isfinite(value)) {
        throw std::runtime_error("photon response probability must be finite");
    }
    return std::clamp(value, 0.0, 1.0);
}

} // namespace

const char* PhotonResponseConfig::modeName() const
{
    return stochastic() ? "stochastic_pe" : "expectation";
}

PhotonResponseConfig buildPhotonResponseConfig(
    const std::map<std::string, std::string>& cfg)
{
    PhotonResponseConfig response;
    std::string mode = lowerCopy(configValue(
        cfg, "response.mode",
        configValue(cfg, "detector.response_mode", response.modeName())));
    if (mode == "stochastic" || mode == "photon" || mode == "discrete_pe") {
        mode = "stochastic_pe";
    } else if (mode == "expected" || mode == "weighted") {
        mode = "expectation";
    }
    if (mode == "stochastic_pe") {
        response.mode = PhotonResponseMode::StochasticPe;
    } else if (mode != "expectation") {
        throw std::runtime_error(
            "response.mode must be expectation or stochastic_pe");
    }
    response.seed = parseSeed(configValue(
        cfg, "response.seed", std::to_string(response.seed)));
    return response;
}

std::uint64_t photonResponseStreamId(const PhotonResponseConfig& cfg,
                                     const PhotonBunch& bunch,
                                     std::uint64_t represented_index)
{
    std::uint64_t seed = cfg.seed;
    seed = mixSeed(seed,
                   static_cast<std::uint64_t>(bunch.shower_event_id + 1000003));
    seed = mixSeed(seed, static_cast<std::uint64_t>(bunch.array_id + 100003));
    seed = mixSeed(seed, static_cast<std::uint64_t>(bunch.telescope_id + 10007));
    seed = mixSeed(seed, bunch.source_bunch_index + 1ULL);
    return mixSeed(seed, represented_index + 1ULL);
}

double photonResponseUniform01(std::uint64_t stream_id,
                               std::uint64_t stage)
{
    std::uint64_t seed = mixSeed(stream_id, stage);
    seed += 0x9e3779b97f4a7c15ULL;
    seed = (seed ^ (seed >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    seed = (seed ^ (seed >> 27U)) * 0x94d049bb133111ebULL;
    seed ^= seed >> 31U;
    constexpr double scale = 1.0 / 9007199254740992.0;
    return static_cast<double>(seed >> 11U) * scale;
}

PhotonResponseSampler::PhotonResponseSampler(PhotonResponseConfig response,
                                             EventIOPhotonConfig eventio)
    : response_(response), eventio_(std::move(eventio))
{
}

std::uint64_t PhotonResponseSampler::candidateCount(
    const PhotonBunch& bunch) const
{
    if (!response_.stochastic()) {
        return 1;
    }
    if (!std::isfinite(bunch.multiplicity) || bunch.multiplicity < 0.0) {
        throw std::runtime_error(
            "response.mode=stochastic_pe requires a finite, non-negative "
            "EventIO bunch multiplicity");
    }
    if (bunch.multiplicity == 0.0) {
        return 0;
    }
    const double count = std::ceil(bunch.multiplicity);
    if (count > static_cast<double>(std::numeric_limits<std::uint64_t>::max())) {
        throw std::runtime_error("EventIO bunch multiplicity is too large to expand");
    }
    return static_cast<std::uint64_t>(count);
}

PhotonCandidate PhotonResponseSampler::candidate(
    const PhotonBunch& bunch,
    std::uint64_t represented_index) const
{
    const std::uint64_t count = candidateCount(bunch);
    if (represented_index >= count) {
        throw std::out_of_range("represented photon index exceeds bunch multiplicity");
    }

    PhotonCandidate out;
    out.photon = bunch.photon;
    out.represented_index = represented_index;
    out.stochastic = response_.stochastic();
    if (!out.stochastic) {
        out.represented_fraction = bunch.multiplicity;
        out.photon.weight *= bunch.multiplicity;
        return out;
    }

    if (std::abs(out.photon.weight - 1.0) > 1.0e-12) {
        throw std::runtime_error(
            "response.mode=stochastic_pe requires source.photon_weight=1 "
            "so accepted photoelectrons remain discrete");
    }
    out.represented_fraction = std::min(
        1.0, bunch.multiplicity - static_cast<double>(represented_index));
    out.remaining_probability = out.represented_fraction;
    out.stream_id = photonResponseStreamId(
        response_, bunch, represented_index);
    out.photon.random_stream_id = out.stream_id;
    out.photon.wavelength_nm = resolveEventIOPhotonWavelength(
        bunch, represented_index, eventio_);
    return out;
}

PreGeometryResponse PhotonResponseSampler::applyPreGeometry(
    PhotonCandidate& candidate,
    double atmosphere_transmission,
    double wavelength_detection_probability) const
{
    PreGeometryResponse result;
    result.expected_weight_before_atmosphere =
        candidate.stochastic
            ? candidate.photon.weight * candidate.represented_fraction
            : candidate.photon.weight;

    const double atmosphere = candidate.photon.optical_efficiency_preapplied
                                  ? 1.0
                                  : probability(atmosphere_transmission);
    result.expected_weight_after_atmosphere =
        result.expected_weight_before_atmosphere * atmosphere;

    if (!candidate.stochastic) {
        if (!candidate.photon.optical_efficiency_preapplied) {
            candidate.photon.weight *= atmosphere;
        }
        result.survives = candidate.photon.weight > 0.0;
        return result;
    }

    if (candidate.photon.optical_efficiency_preapplied) {
        return result;
    }

    const double pretrace_probability = probability(
        candidate.remaining_probability * atmosphere *
        probability(wavelength_detection_probability));
    result.survives = photonResponseUniform01(
        candidate.stream_id, kPreGeometryStage) < pretrace_probability;
    if (result.survives) {
        candidate.remaining_probability = 1.0;
        candidate.photon.optical_efficiency_preapplied = true;
    }
    return result;
}

bool PhotonResponseSampler::acceptPostGeometry(
    const PhotonCandidate& candidate,
    OpticalSurfaceHit& hit) const
{
    if (!candidate.stochastic) {
        return hit.relative_efficiency > 0.0;
    }
    const double detection_probability = probability(
        hit.relative_efficiency * candidate.remaining_probability);
    const bool accepted = photonResponseUniform01(
        candidate.stream_id, kPostGeometryStage) < detection_probability;
    hit.accepted = accepted;
    if (accepted) {
        hit.relative_efficiency = 1.0;
    }
    return accepted;
}

} // namespace lact
