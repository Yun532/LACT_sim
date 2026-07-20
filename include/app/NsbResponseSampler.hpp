#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace lact {

struct NsbConfig;

struct NsbPeSample {
    std::size_t pixel_col = 0;
    std::size_t time_bin = 0;
    float pe = 0.0f;
};

struct NsbRealization {
    std::size_t n_pixels = 0;
    std::size_t n_bins = 0;
    double bin_width_ns = 0.0;
    double expected_total_pe = 0.0;
    std::uint64_t total_pe = 0;
    std::vector<double> rate_pe_per_ns_by_pixel;
    std::vector<float> integrated_pe_by_pixel;
    std::vector<NsbPeSample> samples;
};

std::vector<double> resolveNsbPixelRates(
    const NsbConfig& nsb,
    int telescope_id,
    const std::vector<int>& pixel_ids);

NsbRealization generateNsbRealization(
    const NsbConfig& nsb,
    int event_id,
    int telescope_id,
    const std::vector<int>& pixel_ids,
    std::size_t n_bins,
    double bin_width_ns);

} // namespace lact
