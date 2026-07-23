#include "app/TriggerResponse.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using namespace lact;

namespace {

bool check(bool condition, const char* message)
{
    if (!condition) {
        std::cerr << message << '\n';
    }
    return condition;
}

} // namespace

int main()
{
    bool ok = true;
    TriggerConfig trigger;
    trigger.enabled = true;
    trigger.pixel_threshold_pe = 3.0;
    trigger.camera_multiplicity = 2;
    trigger.camera_coincidence_window_ns = 2.0;
    trigger.array_multiplicity = 2;
    trigger.array_coincidence_window_ns = 10.0;

    const std::vector<std::vector<double>> pe{
        {0.0, 2.0, 2.0, 0.0},
        {0.0, 0.0, 3.0, 0.0},
        {3.0, 0.0, 0.0, 0.0},
    };
    const auto camera = evaluateBinnedPeTrigger(
        pe.size(), pe.front().size(), 1.0, 100.5, trigger,
        [&](std::size_t pixel, std::size_t bin) { return pe[pixel][bin]; });
    ok &= check(camera.triggered, "camera trigger was missed");
    ok &= check(camera.n_pixels_above_threshold == 2,
                "camera multiplicity is incorrect");
    ok &= check(camera.window_start_bin == 1,
                "camera trigger selected the wrong time window");
    ok &= check(std::abs(camera.trigger_time_ns - 102.5) < 1.0e-12,
                "camera trigger time is incorrect");

    const auto array = evaluateArrayTrigger(
        {{1, 100.0}, {2, 107.0}, {3, 130.0}}, trigger);
    ok &= check(array.triggered, "array coincidence was missed");
    ok &= check(array.coincident_telescope_ids == std::vector<int>({1, 2}),
                "array coincidence selected the wrong telescopes");

    const auto corrected_array = evaluateArrayTrigger(
        {{1, 100.0, 100.0}, {2, 250.0, 104.0}, {3, 130.0, 130.0}},
        trigger);
    ok &= check(corrected_array.triggered &&
                    corrected_array.coincident_telescope_ids ==
                        std::vector<int>({1, 2}),
                "array trigger did not use corrected coincidence times");
    ok &= check(std::abs(corrected_array.coincidence_start_time_ns - 100.0) < 1e-12 &&
                    std::abs(corrected_array.coincidence_end_time_ns - 104.0) < 1e-12,
                "array trigger did not report the corrected coincidence interval");

    trigger.array_multiplicity = 3;
    const auto separated = evaluateArrayTrigger(
        {{1, 100.0}, {2, 107.0}, {3, 130.0}}, trigger);
    ok &= check(!separated.triggered,
                "time-separated telescopes incorrectly triggered the array");

    trigger.array_coincidence_window_ns = 0.0;
    const auto no_window = evaluateArrayTrigger(
        {{1, 100.0}, {2, 107.0}, {3, 130.0}}, trigger);
    ok &= check(no_window.triggered &&
                    no_window.coincident_telescope_ids.size() == 3,
                "zero array window should retain legacy count-only behavior");
    return ok ? 0 : 1;
}
