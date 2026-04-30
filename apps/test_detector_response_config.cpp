#include "app/OpticalSimCommon.hpp"

#include <cmath>
#include <iostream>
#include <map>
#include <string>

using namespace lact;

namespace {

bool check(bool cond, const std::string& msg)
{
    if (!cond) {
        std::cerr << msg << "\n";
        return false;
    }
    return cond;
}

bool expectInvalid(const std::map<std::string, std::string>& cfg,
                   const std::string& label)
{
    try {
        (void)buildNsbConfig(cfg);
        (void)buildTriggerConfig(cfg);
    } catch (...) {
        return true;
    }
    std::cerr << "expected invalid config: " << label << "\n";
    return false;
}

} // namespace

int main()
{
    bool ok = true;

    auto default_nsb = buildNsbConfig({});
    ok &= check(!default_nsb.enabled, "NSB should default to disabled");
    ok &= check(default_nsb.model == "constant_rate", "NSB default model");
    ok &= check(std::abs(default_nsb.window_ns - 16.0) < 1e-12, "NSB default window");

    std::map<std::string, std::string> nsb_cfg{
        {"nsb.enabled", "true"},
        {"nsb.model", "constant_rate"},
        {"nsb.rate_pe_per_ns_per_pixel", "0.05"},
        {"nsb.window_ns", "20"},
        {"nsb.seed", "7"},
    };
    auto nsb = buildNsbConfig(nsb_cfg);
    ok &= check(nsb.enabled, "NSB enabled parsed");
    ok &= check(std::abs(nsb.rate_pe_per_ns_per_pixel - 0.05) < 1e-12,
                "NSB rate parsed");
    ok &= check(std::abs(nsb.window_ns - 20.0) < 1e-12, "NSB window parsed");
    ok &= check(nsb.seed == 7, "NSB seed parsed");

    auto trigger_default = buildTriggerConfig({});
    ok &= check(!trigger_default.enabled, "trigger should default to disabled");
    ok &= check(std::abs(trigger_default.pixel_threshold_pe - 5.0) < 1e-12,
                "trigger default threshold");

    std::map<std::string, std::string> trigger_cfg{
        {"trigger.enabled", "true"},
        {"trigger.pixel_threshold_pe", "3.5"},
        {"trigger.camera_multiplicity", "4"},
        {"trigger.array_multiplicity", "2"},
        {"trigger.coincidence_window_ns", "30"},
    };
    auto trigger = buildTriggerConfig(trigger_cfg);
    ok &= check(trigger.enabled, "trigger enabled parsed");
    ok &= check(std::abs(trigger.pixel_threshold_pe - 3.5) < 1e-12,
                "trigger threshold parsed");
    ok &= check(trigger.camera_multiplicity == 4, "camera multiplicity parsed");
    ok &= check(trigger.array_multiplicity == 2, "array multiplicity parsed");
    ok &= check(std::abs(trigger.coincidence_window_ns - 30.0) < 1e-12,
                "coincidence window parsed");

    ok &= expectInvalid({{"nsb.rate_pe_per_ns_per_pixel", "-1"}}, "negative NSB rate");
    ok &= expectInvalid({{"trigger.camera_multiplicity", "0"}}, "zero camera multiplicity");

    std::cout << "detector response config checks passed\n";
    return ok ? 0 : 1;
}
