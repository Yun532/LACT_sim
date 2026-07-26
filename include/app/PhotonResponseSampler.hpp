#pragma once

#include <cstdint>
#include <map>
#include <string>

#include "core/PhotonBunch.hpp"
#include "io/EventIOPhotonSource.hpp"
#include "optics/OpticalSurfaceHit.hpp"

namespace lact {

enum class PhotonResponseMode {
    Expectation,
    StochasticPe,
};

struct PhotonResponseConfig {
    PhotonResponseMode mode = PhotonResponseMode::StochasticPe;
    std::uint64_t seed = 1357913579ULL;

    bool stochastic() const { return mode == PhotonResponseMode::StochasticPe; }
    const char* modeName() const;
};

PhotonResponseConfig buildPhotonResponseConfig(
    const std::map<std::string, std::string>& cfg);

struct PhotonCandidate {
    Photon photon;
    std::uint64_t represented_index = 0;
    std::uint64_t stream_id = 0;
    double represented_fraction = 1.0;
    double remaining_probability = 1.0;
    bool stochastic = false;
};

struct PreGeometryResponse {
    bool survives = true;
    double expected_weight_before_atmosphere = 0.0;
    double expected_weight_after_atmosphere = 0.0;
};

class PhotonResponseSampler {
public:
    PhotonResponseSampler(PhotonResponseConfig response,
                          EventIOPhotonConfig eventio);

    const PhotonResponseConfig& config() const { return response_; }

    std::uint64_t candidateCount(const PhotonBunch& bunch) const;
    PhotonCandidate candidate(const PhotonBunch& bunch,
                              std::uint64_t represented_index) const;

    PreGeometryResponse applyPreGeometry(PhotonCandidate& candidate,
                                         double atmosphere_transmission,
                                         double wavelength_detection_probability) const;

    bool acceptPostGeometry(const PhotonCandidate& candidate,
                            OpticalSurfaceHit& hit) const;

private:
    PhotonResponseConfig response_;
    EventIOPhotonConfig eventio_;
};

std::uint64_t photonResponseStreamId(const PhotonResponseConfig& cfg,
                                     const PhotonBunch& bunch,
                                     std::uint64_t represented_index);
std::uint64_t photonIdentityStreamId(const PhotonBunch& bunch,
                                     std::uint64_t represented_index);

double photonResponseUniform01(std::uint64_t stream_id,
                               std::uint64_t stage);

} // namespace lact
