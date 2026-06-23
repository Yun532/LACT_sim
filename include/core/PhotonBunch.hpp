#pragma once
#include <limits>

#include "core/Photon.hpp"

struct PhotonBunch {
    Photon photon;
    double multiplicity = 1.0;   // how many photons this bunch represents
    int event_id = 0;            // reserved for later
    int telescope_id = -1;       // reserved for later, useful for EventIO/CORSIKA adapters
    bool eventio_2d = false;     // true when the original EventIO block has no explicit z/cz
    double emission_altitude_km = std::numeric_limits<double>::quiet_NaN();
    bool has_emitter = false;    // true when CORSIKA IACT STORE-EMITTER data was attached
    double emitter_mass_gev = std::numeric_limits<double>::quiet_NaN();
    double emitter_charge = std::numeric_limits<double>::quiet_NaN();
    double emitter_energy_gev = std::numeric_limits<double>::quiet_NaN();
    double emitter_time_ns = std::numeric_limits<double>::quiet_NaN();
};
