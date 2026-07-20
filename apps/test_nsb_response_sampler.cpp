#include "app/NsbResponseSampler.hpp"
#include "app/OpticalSimCommon.hpp"

#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using namespace lact;

namespace {

bool check(bool condition, const std::string& message)
{
    if (!condition) {
        std::cerr << message << "\n";
        return false;
    }
    return true;
}

bool sameRealization(const NsbRealization& left, const NsbRealization& right)
{
    if (left.total_pe != right.total_pe ||
        left.samples.size() != right.samples.size() ||
        left.integrated_pe_by_pixel != right.integrated_pe_by_pixel) {
        return false;
    }
    for (std::size_t i = 0; i < left.samples.size(); ++i) {
        if (left.samples[i].pixel_col != right.samples[i].pixel_col ||
            left.samples[i].time_bin != right.samples[i].time_bin ||
            left.samples[i].pe != right.samples[i].pe) {
            return false;
        }
    }
    return true;
}

} // namespace

int main()
{
    bool ok = true;
    const std::vector<int> pixel_ids{10, 20, 30, 40};

    NsbConfig uniform;
    uniform.enabled = true;
    uniform.rate_pe_per_ns_per_pixel = 0.2;
    uniform.seed = 1234567;
    const auto first = generateNsbRealization(
        uniform, 17, 3, pixel_ids, 50, 1.0);
    const auto repeated = generateNsbRealization(
        uniform, 17, 3, pixel_ids, 50, 1.0);
    const auto other_event = generateNsbRealization(
        uniform, 18, 3, pixel_ids, 50, 1.0);
    ok &= check(sameRealization(first, repeated),
                "same NSB key must reproduce exactly");
    ok &= check(!sameRealization(first, other_event),
                "different event IDs should produce different NSB samples");
    ok &= check(std::abs(first.expected_total_pe - 40.0) < 1e-12,
                "uniform expected total p.e. is wrong");
    const double integrated_total = std::accumulate(
        first.integrated_pe_by_pixel.begin(),
        first.integrated_pe_by_pixel.end(), 0.0);
    ok &= check(integrated_total == static_cast<double>(first.total_pe),
                "integrated image must equal sparse waveform total");

    NsbConfig scaled = uniform;
    scaled.spatial_model = "pixel_scale";
    scaled.pixel_relative_scale = {
        {{-1, 10}, 1.0},
        {{-1, 20}, 2.0},
        {{-1, 30}, 0.0},
        {{-1, 40}, 1.0},
        {{3, 40}, 4.0},
    };
    const auto rates = resolveNsbPixelRates(scaled, 3, pixel_ids);
    ok &= check(rates == std::vector<double>({0.2, 0.4, 0.0, 0.8}),
                "pixel scale or telescope override resolution is wrong");
    const auto scaled_realization = generateNsbRealization(
        scaled, 17, 3, pixel_ids, 50, 1.0);
    ok &= check(std::abs(scaled_realization.expected_total_pe - 70.0) < 1e-12,
                "scaled expected total p.e. is wrong");
    ok &= check(scaled_realization.integrated_pe_by_pixel[2] == 0.0f,
                "zero-rate pixel received NSB p.e.");

    const std::vector<int> statistical_pixels{0, 1};
    NsbConfig statistical = uniform;
    statistical.rate_pe_per_ns_per_pixel = 0.1;
    constexpr int trials = 5000;
    double sum0 = 0.0;
    double sum1 = 0.0;
    double sum00 = 0.0;
    double sum11 = 0.0;
    double sum01 = 0.0;
    for (int event_id = 0; event_id < trials; ++event_id) {
        const auto sample = generateNsbRealization(
            statistical, event_id, 0, statistical_pixels, 10, 1.0);
        const double x0 = sample.integrated_pe_by_pixel[0];
        const double x1 = sample.integrated_pe_by_pixel[1];
        sum0 += x0;
        sum1 += x1;
        sum00 += x0 * x0;
        sum11 += x1 * x1;
        sum01 += x0 * x1;
    }
    const double mean0 = sum0 / trials;
    const double mean1 = sum1 / trials;
    const double variance0 = sum00 / trials - mean0 * mean0;
    const double variance1 = sum11 / trials - mean1 * mean1;
    const double covariance = sum01 / trials - mean0 * mean1;
    ok &= check(std::abs(mean0 - 1.0) < 0.06 && std::abs(mean1 - 1.0) < 0.06,
                "uniform NSB pixel means are inconsistent with Poisson expectation");
    ok &= check(std::abs(variance0 - 1.0) < 0.10 &&
                    std::abs(variance1 - 1.0) < 0.10,
                "uniform NSB pixel variances are inconsistent with Poisson expectation");
    ok &= check(std::abs(covariance) < 0.06,
                "uniform NSB pixel covariance is too large");

    try {
        NsbConfig incomplete = scaled;
        incomplete.pixel_relative_scale.erase({-1, 20});
        (void)resolveNsbPixelRates(incomplete, 3, pixel_ids);
        std::cerr << "missing pixel scale should fail\n";
        ok = false;
    } catch (...) {
    }

    if (ok) {
        std::cout << "NSB response sampler checks passed\n";
    }
    return ok ? 0 : 1;
}
