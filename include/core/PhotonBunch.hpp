#pragma once
#include "core/Photon.hpp"

struct PhotonBunch {
    Photon photon;
    double multiplicity = 1.0;   // how many photons this bunch represents
    int event_id = 0;            // reserved for later
    int telescope_id = -1;       // reserved for later, useful for EventIO/CORSIKA adapters
    bool eventio_2d = false;     // true when the original EventIO block has no explicit z/cz
};
