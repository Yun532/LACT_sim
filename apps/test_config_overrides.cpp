#include "app/OpticalSimCommon.hpp"
#include "app/CorsikaTraceConfig.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
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

    const auto original_cwd = std::filesystem::current_path();
    const auto path_test_root = std::filesystem::temp_directory_path() /
        ("lact_config_path_" + std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        const auto project = path_test_root / "project";
        const auto config_dir = project / "configs" / "examples";
        const auto asset = project / "configs" / "assets" / "input.csv";
        const auto cwd = path_test_root / "cwd";
        std::filesystem::create_directories(config_dir);
        std::filesystem::create_directories(asset.parent_path());
        std::filesystem::create_directories(cwd / "cwd-only");
        std::ofstream(config_dir / "main.cfg") << "source.mode=PhotonCsv\n";
        std::ofstream(asset) << "x,y\n0,0\n";
        std::ofstream(cwd / "cwd-only" / "input.csv") << "x,y\n0,0\n";

        const auto root_relative = resolveRelativePath(
            (config_dir / "main.cfg").string(),
            "configs/assets/input.csv");
        ok &= check(std::filesystem::equivalent(root_relative, asset),
                    "repository-root-style config path fallback failed");

        const auto owner_relative = resolveRelativePath(
            (config_dir / "main.cfg").string(),
            "../assets/input.csv");
        ok &= check(std::filesystem::equivalent(owner_relative, asset),
                    "config-directory-relative path resolution changed");

        const auto absolute = resolveRelativePath(
            (config_dir / "main.cfg").string(), asset.string());
        ok &= check(std::filesystem::equivalent(absolute, asset),
                    "absolute input path was prefixed or rewritten");

        std::filesystem::current_path(cwd);
        const auto cwd_relative = resolveRelativePath(
            (config_dir / "main.cfg").string(), "cwd-only/input.csv");
        ok &= check(std::filesystem::equivalent(
                        cwd_relative, cwd / "cwd-only" / "input.csv"),
                    "working-directory path fallback failed");

        std::filesystem::create_directories(config_dir / "ambiguous");
        std::filesystem::create_directories(cwd / "ambiguous");
        std::ofstream(config_dir / "ambiguous" / "input.csv") << "owner\n";
        std::ofstream(cwd / "ambiguous" / "input.csv") << "cwd\n";
        try {
            (void)resolveRelativePath(
                (config_dir / "main.cfg").string(), "ambiguous/input.csv");
            std::cerr << "ambiguous config input path was accepted\n";
            ok = false;
        } catch (const std::runtime_error&) {
        }
    } catch (const std::exception& error) {
        std::cerr << "config path test setup failed: " << error.what() << '\n';
        ok = false;
    }
    std::filesystem::current_path(original_cwd);
    std::error_code cleanup_error;
    std::filesystem::remove_all(path_test_root, cleanup_error);

    return ok ? 0 : 1;
}
