#include "app/TriggerResponse.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lact {
namespace {

double coincidenceTime(const TelescopeTriggerTime& item)
{
    return std::isfinite(item.coincidence_time_ns)
        ? item.coincidence_time_ns
        : item.trigger_time_ns;
}

} // namespace

BinnedPeTriggerDecision evaluateBinnedPeTrigger(
    std::size_t n_pixels,
    std::size_t n_bins,
    double bin_width_ns,
    double first_bin_center_time_ns,
    const TriggerConfig& trigger,
    const std::function<double(std::size_t, std::size_t)>& pe_at)
{
    BinnedPeTriggerDecision out;
    if (!trigger.enabled || n_pixels == 0 || n_bins == 0 || !pe_at) {
        return out;
    }
    if (!(bin_width_ns > 0.0) || !std::isfinite(bin_width_ns)) {
        throw std::runtime_error("trigger bin width must be finite and > 0");
    }

    const std::size_t window_bins =
        trigger.camera_coincidence_window_ns > 0.0
            ? std::min(n_bins, std::max<std::size_t>(
                                   1,
                                   static_cast<std::size_t>(std::ceil(
                                       trigger.camera_coincidence_window_ns /
                                       bin_width_ns))))
            : n_bins;

    std::vector<int> pixels_above_threshold(n_bins, 0);
    for (std::size_t col = 0; col < n_pixels; ++col) {
        double window_pe = 0.0;
        for (std::size_t bin = 0; bin < window_bins; ++bin) {
            window_pe += pe_at(col, bin);
        }
        if (window_pe >= trigger.pixel_threshold_pe) {
            ++pixels_above_threshold[0];
        }
        for (std::size_t start = 1; start < n_bins; ++start) {
            window_pe -= pe_at(col, start - 1);
            const std::size_t add_bin = start + window_bins - 1;
            if (add_bin < n_bins) {
                window_pe += pe_at(col, add_bin);
            }
            if (window_pe >= trigger.pixel_threshold_pe) {
                ++pixels_above_threshold[start];
            }
        }
    }

    for (std::size_t start = 0; start < n_bins; ++start) {
        if (pixels_above_threshold[start] > out.n_pixels_above_threshold) {
            out.n_pixels_above_threshold = pixels_above_threshold[start];
            out.window_start_bin = start;
        }
    }
    out.triggered =
        out.n_pixels_above_threshold >= trigger.camera_multiplicity;
    const std::size_t center_offset = std::min(
        n_bins - 1 - out.window_start_bin, window_bins / 2);
    out.trigger_time_ns = first_bin_center_time_ns +
                          static_cast<double>(out.window_start_bin + center_offset) *
                              bin_width_ns;
    return out;
}

ArrayTriggerDecision evaluateArrayTrigger(
    std::vector<TelescopeTriggerTime> telescope_triggers,
    const TriggerConfig& trigger)
{
    ArrayTriggerDecision out;
    if (!trigger.enabled || telescope_triggers.empty()) {
        return out;
    }
    telescope_triggers.erase(
        std::remove_if(telescope_triggers.begin(), telescope_triggers.end(),
                       [](const TelescopeTriggerTime& item) {
                           return !std::isfinite(coincidenceTime(item));
                       }),
        telescope_triggers.end());
    if (telescope_triggers.empty()) {
        return out;
    }

    std::sort(telescope_triggers.begin(), telescope_triggers.end(),
              [](const TelescopeTriggerTime& a, const TelescopeTriggerTime& b) {
                  const double a_time = coincidenceTime(a);
                  const double b_time = coincidenceTime(b);
                  if (a_time != b_time) {
                      return a_time < b_time;
                  }
                  return a.telescope_id < b.telescope_id;
              });

    std::size_t best_begin = 0;
    std::size_t best_end = 0;
    if (trigger.array_coincidence_window_ns <= 0.0) {
        best_end = telescope_triggers.size() - 1;
    } else {
        std::size_t begin = 0;
        for (std::size_t end = 0; end < telescope_triggers.size(); ++end) {
            while (begin < end &&
                   coincidenceTime(telescope_triggers[end]) -
                           coincidenceTime(telescope_triggers[begin]) >
                       trigger.array_coincidence_window_ns) {
                ++begin;
            }
            if (end - begin > best_end - best_begin) {
                best_begin = begin;
                best_end = end;
            }
        }
    }

    const std::size_t count = best_end - best_begin + 1;
    out.triggered = static_cast<int>(count) >= trigger.array_multiplicity;
    if (!out.triggered) {
        return out;
    }
    out.coincidence_start_time_ns =
        coincidenceTime(telescope_triggers[best_begin]);
    out.coincidence_end_time_ns =
        coincidenceTime(telescope_triggers[best_end]);
    out.coincident_telescope_ids.reserve(count);
    for (std::size_t i = best_begin; i <= best_end; ++i) {
        out.coincident_telescope_ids.push_back(telescope_triggers[i].telescope_id);
    }
    std::sort(out.coincident_telescope_ids.begin(),
              out.coincident_telescope_ids.end());
    return out;
}

} // namespace lact
