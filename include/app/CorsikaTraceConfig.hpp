#pragma once

// Builders that turn the flat key/value config into the typed
// output, waveform, profiling and histogram settings.

#include <map>
#include <string>

#include "app/OpticalSimCommon.hpp"
#include "io/CorsikaTraceStats.hpp"

namespace lact {

CorsikaTraceOutputConfig buildCorsikaTraceOutputConfig(
    const std::map<std::string, std::string>& cfg);

std::string normalizeWaveformTimeReference(std::string value);

WaveformOutputConfig buildWaveformOutputConfig(
    const std::map<std::string, std::string>& cfg);

CollectorDebugConfig buildCollectorDebugConfig(
    const std::map<std::string, std::string>& cfg);

ProfileConfig buildProfileConfig(const std::map<std::string, std::string>& cfg);

AtmosphereHistogramConfig buildAtmosphereHistogramConfig(
    const std::map<std::string, std::string>& cfg);

} // namespace lact
