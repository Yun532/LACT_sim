#include "app/TriggerResponse.hpp"
#include "io/CorsikaTraceOutputTypes.hpp"

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
    ok &= check(camera.first_trigger_window_start_bin == 1,
                "camera trigger selected the wrong first threshold window");
    ok &= check(std::abs(camera.trigger_time_ns - 102.5) < 1.0e-12,
                "camera trigger time is incorrect");
    ok &= check(std::abs(camera.first_trigger_time_ns - 102.5) < 1.0e-12 &&
                    std::abs(camera.max_multiplicity_time_ns - 102.5) < 1.0e-12,
                "camera trigger diagnostic times are incorrect");

    const std::vector<std::vector<double>> growing_pe{
        {3.0, 0.0, 3.0, 0.0},
        {3.0, 0.0, 3.0, 0.0},
        {0.0, 0.0, 3.0, 0.0},
    };
    trigger.camera_coincidence_window_ns = 1.0;
    const auto growing_camera = evaluateBinnedPeTrigger(
        growing_pe.size(), growing_pe.front().size(), 1.0, 100.5, trigger,
        [&](std::size_t pixel, std::size_t bin) {
            return growing_pe[pixel][bin];
        });
    ok &= check(growing_camera.triggered,
                "growing camera trigger was missed");
    ok &= check(growing_camera.first_trigger_window_start_bin == 0 &&
                    growing_camera.window_start_bin == 2,
                "first and maximum-multiplicity windows were not separated");
    ok &= check(std::abs(growing_camera.trigger_time_ns - 100.5) < 1.0e-12 &&
                    std::abs(growing_camera.first_trigger_time_ns - 100.5) <
                        1.0e-12 &&
                    std::abs(growing_camera.max_multiplicity_time_ns - 102.5) <
                        1.0e-12,
                "first and maximum-multiplicity times are incorrect");
    ok &= check(growing_camera.n_pixels_above_threshold == 3,
                "maximum camera multiplicity was not retained");
    trigger.camera_coincidence_window_ns = 2.0;

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

    NsbConfig nsb;
    nsb.enabled = true;
    nsb.rate_pe_per_ns_per_pixel = 1.5;
    nsb.seed = 1234567;
    WaveformOutputConfig waveform;
    waveform.enabled = true;
    waveform.source = "pe";
    waveform.time_bin_width_ns = 1.0;
    constexpr std::size_t nsb_pixels = 8;
    constexpr std::size_t nsb_bins = 6;
    const auto generate_nsb = [&](int event_id) {
        std::vector<double> cells(nsb_pixels * nsb_bins, 0.0);
        generateTimeBinnedNsbPe(
            nsb, waveform, event_id, 2, nsb_pixels, nsb_bins,
            [&](std::size_t pixel, std::size_t bin, float value) {
                cells[pixel * nsb_bins + bin] += value;
            });
        return cells;
    };
    const auto nsb_first = generate_nsb(101);
    const auto nsb_repeat = generate_nsb(101);
    const auto nsb_other_event = generate_nsb(102);
    ok &= check(nsb_first == nsb_repeat,
                "time-binned NSB realization is not reproducible");
    ok &= check(nsb_first != nsb_other_event,
                "time-binned NSB realization did not separate events");

    trigger.pixel_threshold_pe = 3.0;
    trigger.camera_multiplicity = 2;
    trigger.camera_coincidence_window_ns = 2.0;
    const auto nsb_trigger = evaluateBinnedPeTrigger(
        nsb_pixels, nsb_bins, waveform.time_bin_width_ns, 0.5, trigger,
        [&](std::size_t pixel, std::size_t bin) {
            return nsb_first[pixel * nsb_bins + bin];
        });
    const auto nsb_trigger_repeat = evaluateBinnedPeTrigger(
        nsb_pixels, nsb_bins, waveform.time_bin_width_ns, 0.5, trigger,
        [&](std::size_t pixel, std::size_t bin) {
            return nsb_repeat[pixel * nsb_bins + bin];
        });
    ok &= check(nsb_trigger.triggered == nsb_trigger_repeat.triggered &&
                    nsb_trigger.n_pixels_above_threshold ==
                        nsb_trigger_repeat.n_pixels_above_threshold &&
                    nsb_trigger.trigger_time_ns ==
                        nsb_trigger_repeat.trigger_time_ns,
                "trigger did not reproduce the saved NSB realization");
    return ok ? 0 : 1;
}
