#pragma once

#include <cstdint>
#include <map>
#include <vector>

#include <hdf5.h>

#include "app/OpticalSimCommon.hpp"
#include "io/CorsikaTraceOutputTypes.hpp"
#include "io/CameraElectronicsEvent.hpp"

namespace lact {

struct Hdf5WaveformImage {
    std::int32_t image_index = 0;
    int event_id = 0;
    int telescope_id = 0;
    double time_first_ns = 0.0;
    double time_mean_ns = 0.0;
};

void writeHdf5Waveforms(
    hid_t file,
    const CorsikaTraceOutputConfig& output,
    const WaveformOutputConfig& waveform,
    const ElectronicsConfig& electronics,
    const NsbConfig& nsb,
    const std::vector<std::int32_t>& pixel_id_axis,
    const std::vector<Hdf5WaveformImage>& images,
    const std::map<WaveformKey, WaveformPixelAccumulator>& waveforms,
    const std::vector<RawWaveformHit>& raw_hits,
    const CameraElectronicsEventMap& electronics_events);

} // namespace lact
