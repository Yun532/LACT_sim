#include "app/OpticalSimCommon.hpp"

#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

int main()
{
    try {
        std::map<std::string, std::string> values{
            {"electronics.sipm.channels_per_pixel", "8"},
            {"electronics.sipm.microcells_per_channel", "33792"},
            {"electronics.sipm.saturation_enabled", "true"},
            {"electronics.sipm.saturation_model", "hard_no_recovery"},
        };
        const auto config = lact::buildElectronicsConfig(values);
        const lact::ElectronicsResponse response(config);
        const double cells = 8.0 * 33792.0;
        const double primary = cells;
        const double expected = cells * (1.0 - std::exp(-1.0));
        const double actual = response.saturatedPe(primary);
        if (std::abs(actual - expected) > 1e-9 * cells) {
            throw std::runtime_error("hard_no_recovery saturation mismatch");
        }
        if (response.saturatedPe(0.0) != 0.0) {
            throw std::runtime_error("zero p.e. saturation mismatch");
        }
        const std::vector<double> primary_waveform{
            0.05 * cells, 0.15 * cells, 0.30 * cells, 0.50 * cells};
        const auto fired_waveform =
            response.saturatedWaveform(primary_waveform);
        double primary_sum = 0.0;
        double fired_sum = 0.0;
        for (std::size_t i = 0; i < primary_waveform.size(); ++i) {
            primary_sum += primary_waveform[i];
            fired_sum += fired_waveform[i];
            if (fired_waveform[i] < 0.0 ||
                fired_waveform[i] > primary_waveform[i]) {
                throw std::runtime_error(
                    "time-binned saturation produced an invalid fired p.e. sample");
            }
        }
        if (std::abs(fired_sum - response.saturatedPe(primary_sum)) >
            1e-9 * cells) {
            throw std::runtime_error(
                "fired waveform integral must match integrated saturation");
        }

        values["electronics.sipm.saturation_enabled"] = "false";
        const lact::ElectronicsResponse linear(
            lact::buildElectronicsConfig(values));
        if (linear.saturatedPe(primary) != primary) {
            throw std::runtime_error("disabled saturation must be linear");
        }
        if (linear.saturatedWaveform(primary_waveform) != primary_waveform) {
            throw std::runtime_error(
                "disabled saturation must preserve every waveform sample");
        }
        std::cout << "electronics saturation checks passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
