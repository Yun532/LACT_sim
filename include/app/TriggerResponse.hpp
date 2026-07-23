#pragma once

#include <cstddef>
#include <functional>
#include <limits>
#include <vector>

#include "app/OpticalSimCommon.hpp"

namespace lact {

struct BinnedPeTriggerDecision {
    bool triggered = false;
    int n_pixels_above_threshold = 0;
    std::size_t window_start_bin = 0;
    double trigger_time_ns = 0.0;
};

struct TelescopeTriggerTime {
    int telescope_id = -1;
    // Raw local camera-trigger time. Kept for output/debug compatibility.
    double trigger_time_ns = 0.0;
    // Time used by the array coincidence test. NaN falls back to the raw
    // trigger time, preserving the previous interface for callers that do
    // not request geometric correction.
    double coincidence_time_ns = std::numeric_limits<double>::quiet_NaN();
    double geometric_delay_ns = 0.0;
};

struct ArrayTriggerDecision {
    bool triggered = false;
    std::vector<int> coincident_telescope_ids;
    double coincidence_start_time_ns =
        std::numeric_limits<double>::quiet_NaN();
    double coincidence_end_time_ns =
        std::numeric_limits<double>::quiet_NaN();
};

BinnedPeTriggerDecision evaluateBinnedPeTrigger(
    std::size_t n_pixels,
    std::size_t n_bins,
    double bin_width_ns,
    double first_bin_center_time_ns,
    const TriggerConfig& trigger,
    const std::function<double(std::size_t, std::size_t)>& pe_at);

ArrayTriggerDecision evaluateArrayTrigger(
    std::vector<TelescopeTriggerTime> telescope_triggers,
    const TriggerConfig& trigger);

} // namespace lact
