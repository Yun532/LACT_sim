#pragma once
#include <cstdint>
#include <random>
#include "io/PhotonSource.hpp"

enum class SyntheticMode {
    ParallelBeam, // parallel rays entering from a plane
    PointSource   // rays emitted from one point toward a telescope aperture
};

struct SyntheticPhotonConfig {
    SyntheticMode mode = SyntheticMode::ParallelBeam;

    int n_bunches = 1000;
    double multiplicity = 1.0;

    double wavelength_nm = 400.0;
    double time_ns = 0.0;
    double photon_weight = 1.0;

    // For ParallelBeam:
    // rays start from z = source_plane_z and are sampled uniformly in a disk of beam_radius_m
    double source_plane_z = 100.0;
    double beam_radius_m = 1.0;
    Vec3 beam_direction{0.0, 0.0, -1.0};

    // For PointSource:
    Vec3 source_position{0.0, 0.0, 100.0};
    double aperture_z = 0.0;
    double aperture_radius_m = 1.0;

    // Common
    int event_id = 0;
    int telescope_id = 0;
    std::uint64_t random_seed = 123456789ULL;
};

class SyntheticPhotonSource : public PhotonSource {
public:
    explicit SyntheticPhotonSource(const SyntheticPhotonConfig& cfg);

    bool next(PhotonBunch& out) override;
    void reset() override;

private:
    SyntheticPhotonConfig cfg_;
    int index_ = 0;
    std::mt19937_64 rng_;
    std::uniform_real_distribution<double> uni_{0.0, 1.0};

    Vec3 sampleDisk(double radius);
    PhotonBunch makeParallelBeamBunch();
    PhotonBunch makePointSourceBunch();
};
