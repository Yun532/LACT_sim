#include "io/SyntheticPhotonSource.hpp"
#include <cmath>
#include <stdexcept>

namespace {
constexpr double PI = 3.14159265358979323846;
}

SyntheticPhotonSource::SyntheticPhotonSource(const SyntheticPhotonConfig& cfg)
    : cfg_(cfg), index_(0), rng_(cfg.random_seed) {
    if (cfg_.n_bunches < 0) {
        throw std::runtime_error("SyntheticPhotonSource: n_bunches must be >= 0");
    }
    if (!std::isfinite(cfg_.multiplicity) || cfg_.multiplicity < 0.0) {
        throw std::runtime_error("SyntheticPhotonSource: multiplicity must be finite and >= 0");
    }
    if (!std::isfinite(cfg_.photon_weight) || cfg_.photon_weight < 0.0) {
        throw std::runtime_error("SyntheticPhotonSource: photon_weight must be finite and >= 0");
    }
    if (!std::isfinite(cfg_.wavelength_nm) || cfg_.wavelength_nm <= 0.0) {
        throw std::runtime_error("SyntheticPhotonSource: wavelength_nm must be finite and > 0");
    }
    if (!std::isfinite(cfg_.time_ns)) {
        throw std::runtime_error("SyntheticPhotonSource: time_ns must be finite");
    }
    if (!std::isfinite(cfg_.source_plane_z)) {
        throw std::runtime_error("SyntheticPhotonSource: source_plane_z must be finite");
    }
    if (!std::isfinite(cfg_.beam_radius_m) || cfg_.beam_radius_m < 0.0) {
        throw std::runtime_error("SyntheticPhotonSource: beam_radius_m must be finite and >= 0");
    }
    if (!std::isfinite(cfg_.beam_direction.x) ||
        !std::isfinite(cfg_.beam_direction.y) ||
        !std::isfinite(cfg_.beam_direction.z) ||
        cfg_.beam_direction.norm() <= 0.0) {
        throw std::runtime_error("SyntheticPhotonSource: beam_direction must be finite and non-zero");
    }
    if (!std::isfinite(cfg_.source_position.x) ||
        !std::isfinite(cfg_.source_position.y) ||
        !std::isfinite(cfg_.source_position.z)) {
        throw std::runtime_error("SyntheticPhotonSource: source_position must be finite");
    }
    if (!std::isfinite(cfg_.aperture_z)) {
        throw std::runtime_error("SyntheticPhotonSource: aperture_z must be finite");
    }
    if (!std::isfinite(cfg_.aperture_radius_m) || cfg_.aperture_radius_m < 0.0) {
        throw std::runtime_error("SyntheticPhotonSource: aperture_radius_m must be finite and >= 0");
    }
}

void SyntheticPhotonSource::reset() {
    index_ = 0;
    rng_.seed(cfg_.random_seed);
}

bool SyntheticPhotonSource::next(PhotonBunch& out) {
    if (index_ >= cfg_.n_bunches) {
        return false;
    }

    switch (cfg_.mode) {
        case SyntheticMode::ParallelBeam:
            out = makeParallelBeamBunch();
            break;
        case SyntheticMode::PointSource:
            out = makePointSourceBunch();
            break;
        default:
            throw std::runtime_error("SyntheticPhotonSource: unsupported synthetic mode");
    }

    // Preserve a stable, unique identity for every generated bunch.  Optical
    // stochastic stages derive their per-photon streams from this index; if it
    // remains at the PhotonBunch default (zero), every synthetic photon that
    // hits the same facet receives the same roughness deflection.
    out.source_bunch_index = static_cast<std::uint64_t>(index_);
    ++index_;
    return true;
}

Vec3 SyntheticPhotonSource::sampleDisk(double radius) {
    // Uniform in area
    double u1 = uni_(rng_);
    double u2 = uni_(rng_);

    double r = radius * std::sqrt(u1);
    double phi = 2.0 * PI * u2;

    return {r * std::cos(phi), r * std::sin(phi), 0.0};
}

PhotonBunch SyntheticPhotonSource::makeParallelBeamBunch() {
    PhotonBunch b;
    b.event_id = cfg_.event_id;
    b.telescope_id = cfg_.telescope_id;
    b.multiplicity = cfg_.multiplicity;

    Vec3 xy = sampleDisk(cfg_.beam_radius_m);

    b.photon.pos = {xy.x, xy.y, cfg_.source_plane_z};
    b.photon.dir = cfg_.beam_direction.normalized();
    b.photon.wavelength_nm = cfg_.wavelength_nm;
    b.photon.time_ns = cfg_.time_ns;
    b.photon.weight = cfg_.photon_weight;

    return b;
}

PhotonBunch SyntheticPhotonSource::makePointSourceBunch() {
    PhotonBunch b;
    b.event_id = cfg_.event_id;
    b.telescope_id = cfg_.telescope_id;
    b.multiplicity = cfg_.multiplicity;

    Vec3 target_xy = sampleDisk(cfg_.aperture_radius_m);
    Vec3 target{target_xy.x, target_xy.y, cfg_.aperture_z};

    Vec3 dir = (target - cfg_.source_position).normalized();

    b.photon.pos = cfg_.source_position;
    b.photon.dir = dir;
    b.photon.wavelength_nm = cfg_.wavelength_nm;
    b.photon.time_ns = cfg_.time_ns;
    b.photon.weight = cfg_.photon_weight;

    return b;
}
