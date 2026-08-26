#pragma once
#include <cstdint>
#include <limits>

#include "core/Photon.hpp"

enum class PhotonOrigin {
    Cherenkov,
    Nsb,
    Dark,
};

struct PhotonBunch {
    Photon photon;
    double multiplicity = 1.0;   // how many photons this bunch represents
    int event_id = 0;            // compatibility/output event identifier
    int shower_event_id = 0;     // original CORSIKA shower identifier
    int array_id = 0;            // original MC_TELARRAY/core-reuse identifier
    int telescope_id = -1;
    std::uint64_t source_bunch_index = 0;
    double raw_wavelength_nm = 400.0;
    PhotonOrigin origin = PhotonOrigin::Cherenkov;
    bool eventio_2d = false;     // true when the original EventIO block has no explicit z/cz
    double emission_altitude_km = std::numeric_limits<double>::quiet_NaN();
    bool has_emitter = false;    // true when CORSIKA IACT STORE-EMITTER data was attached
    double emitter_mass_gev = std::numeric_limits<double>::quiet_NaN();
    double emitter_charge = std::numeric_limits<double>::quiet_NaN();
    double emitter_energy_gev = std::numeric_limits<double>::quiet_NaN();
    double emitter_time_ns = std::numeric_limits<double>::quiet_NaN();
};
