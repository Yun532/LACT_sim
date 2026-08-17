#pragma once

// EventIO input adaptation: frame conversion for incoming bunches and
// the wavelength/atmosphere metadata carried on the CORSIKA input card.

#include <map>
#include <optional>
#include <string>

#include "app/OpticalSimCommon.hpp"
#include "io/CorsikaTraceStats.hpp"

namespace lact {

PhotonBunch transformEventIOBunchToTraceFrame(
    const PhotonBunch& input,
    const TelescopeConfig& telescope_cfg,
    const EventIOMetadata& metadata,
    const SourceRuntimeConfig& source_runtime_cfg);

std::optional<WavelengthRange> wavelengthRangeFromInputCard(
    const EventIOMetadata& metadata);

bool hasExplicitMissingWavelengthRange(const std::map<std::string, std::string>& cfg);

void applyEventIOWavelengthMetadata(EventIOPhotonConfig& eventio_cfg,
                                    const EventIOMetadata& metadata,
                                    const std::map<std::string, std::string>& cfg);

void applyEventIOAtmosphereMetadata(EventIOPhotonConfig& eventio_cfg,
                                    const EventIOMetadata& metadata);

} // namespace lact
