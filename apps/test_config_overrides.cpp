#include "app/OpticalSimCommon.hpp"
#include "app/CorsikaTraceConfig.hpp"

#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

using namespace lact;

namespace {

bool check(bool condition, const std::string& message)
{
    if (!condition) std::cerr << message << '\n';
    return condition;
}

template <std::size_t N>
ConfigCommandLine parse(const char* const (&values)[N])
{
    char* argv[N];
    for (std::size_t i = 0; i < N; ++i) {
        argv[i] = const_cast<char*>(values[i]);
    }
    return parseConfigCommandLine(static_cast<int>(N), argv);
}

} // namespace

int main()
{
    bool ok = true;
    const char* args[] = {
        "run_corsika_trace",
        "production.cfg",
        "-C", "nsb.enabled=true",
        "-Celectronics.microcell.saturation_enabled=false",
        "--set=output.lact_root_path=batch/run_17.root",
        "--set", "nsb.seed=17",
        "-C", "nsb.seed=18",
    };
    const auto parsed = parse(args);
    ok &= check(parsed.positional == std::vector<std::string>{"production.cfg"},
                "configuration positional argument was not preserved");
    ok &= check(parsed.overrides.size() == 5,
                "not all repeated configuration overrides were parsed");

    std::map<std::string, std::string> cfg{
        {"nsb.enabled", "false"},
        {"nsb.seed", "1"},
    };
    applyConfigOverrides(cfg, parsed.overrides);
    ok &= check(cfg["nsb.enabled"] == "true", "boolean override failed");
    ok &= check(cfg["nsb.seed"] == "18",
                "later repeated override must win");
    ok &= check(cfg["electronics.microcell.saturation_enabled"] == "false",
                "compact -Ckey=value syntax failed");
    ok &= check(cfg["output.lact_root_path"] == "batch/run_17.root",
                "--set=key=value syntax failed");

    const auto default_source = buildSourceRuntimeConfig({});
    ok &= check(default_source.event_id_mode == "event",
                "default EventIO identity must remain backward compatible");
    const auto array_source = buildSourceRuntimeConfig({
        {"source.event_id_mode", "event_array100"},
    });
    ok &= check(array_source.event_id_mode == "event_array100",
                "array/core identity mode must remain explicitly selectable");

    try {
        (void)buildWaveformOutputConfig({
            {"waveform.enabled", "true"},
            {"waveform.source", "pe"},
            {"waveform.time_window_start_ns", "0"},
            {"waveform.time_window_end_ns", "10"},
            {"waveform.time_bin_width_ns", "4"},
        });
        std::cerr << "non-integral waveform sampling grid was accepted\n";
        ok = false;
    } catch (const std::runtime_error&) {
    }

    try {
        const char* invalid[] = {"program", "config.cfg", "-C", "missing_equal"};
        (void)parse(invalid);
        std::cerr << "override without '=' was accepted\n";
        ok = false;
    } catch (...) {
    }
    try {
        const char* invalid[] = {"program", "config.cfg", "--unknown"};
        (void)parse(invalid);
        std::cerr << "unknown option was accepted\n";
        ok = false;
    } catch (...) {
    }

    return ok ? 0 : 1;
}
