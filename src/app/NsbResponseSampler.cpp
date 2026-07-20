#include "app/NsbResponseSampler.hpp"

#include "app/OpticalSimCommon.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>
#include <utility>

namespace lact {
namespace {

std::uint64_t mixSeed(std::uint64_t seed, std::uint64_t value)
{
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    value ^= value >> 31U;
    seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
    return seed;
}

std::mt19937_64 makeRng(const NsbConfig& nsb, int event_id, int telescope_id)
{
    std::uint64_t seed = nsb.seed;
    seed = mixSeed(seed, static_cast<std::uint64_t>(
                             static_cast<std::int64_t>(event_id) + 0x80000000LL));
    seed = mixSeed(seed, static_cast<std::uint64_t>(telescope_id + 1));
    return std::mt19937_64(seed);
}

} // namespace

std::vector<double> resolveNsbPixelRates(
    const NsbConfig& nsb,
    int telescope_id,
    const std::vector<int>& pixel_ids)
{
    std::vector<double> rates(pixel_ids.size(), nsb.rate_pe_per_ns_per_pixel);
    if (nsb.spatial_model == "uniform") {
        return rates;
    }
    if (nsb.spatial_model != "pixel_scale") {
        throw std::runtime_error("unsupported NSB spatial model: " + nsb.spatial_model);
    }

    for (std::size_t col = 0; col < pixel_ids.size(); ++col) {
        const int pixel_id = pixel_ids[col];
        const auto exact = nsb.pixel_relative_scale.find({telescope_id, pixel_id});
        const auto fallback = nsb.pixel_relative_scale.find({-1, pixel_id});
        if (exact != nsb.pixel_relative_scale.end()) {
            rates[col] *= exact->second;
        } else if (fallback != nsb.pixel_relative_scale.end()) {
            rates[col] *= fallback->second;
        } else {
            throw std::runtime_error(
                "NSB pixel_scale has no row for telescope " +
                std::to_string(telescope_id) + ", pixel " +
                std::to_string(pixel_id));
        }
    }
    return rates;
}

NsbRealization generateNsbRealization(
    const NsbConfig& nsb,
    int event_id,
    int telescope_id,
    const std::vector<int>& pixel_ids,
    std::size_t n_bins,
    double bin_width_ns)
{
    NsbRealization out;
    out.n_pixels = pixel_ids.size();
    out.n_bins = n_bins;
    out.bin_width_ns = bin_width_ns;
    out.integrated_pe_by_pixel.assign(pixel_ids.size(), 0.0f);

    if (!nsb.enabled || pixel_ids.empty() || n_bins == 0 ||
        nsb.rate_pe_per_ns_per_pixel <= 0.0) {
        return out;
    }
    if (!(bin_width_ns > 0.0) || !std::isfinite(bin_width_ns)) {
        throw std::runtime_error("NSB bin width must be finite and > 0");
    }

    out.rate_pe_per_ns_by_pixel =
        resolveNsbPixelRates(nsb, telescope_id, pixel_ids);
    double rate_sum = 0.0;
    bool uniform_rates = true;
    double first_rate = out.rate_pe_per_ns_by_pixel.front();
    for (double rate : out.rate_pe_per_ns_by_pixel) {
        if (!std::isfinite(rate) || rate < 0.0) {
            throw std::runtime_error("resolved NSB pixel rate must be finite and >= 0");
        }
        rate_sum += rate;
        uniform_rates = uniform_rates && rate == first_rate;
    }
    out.expected_total_pe =
        rate_sum * bin_width_ns * static_cast<double>(n_bins);
    if (!(out.expected_total_pe > 0.0)) {
        return out;
    }
    if (!std::isfinite(out.expected_total_pe)) {
        throw std::runtime_error("total expected NSB p.e. must be finite");
    }

    auto rng = makeRng(nsb, event_id, telescope_id);
    std::poisson_distribution<unsigned long long> total_poisson(
        out.expected_total_pe);
    out.total_pe = static_cast<std::uint64_t>(total_poisson(rng));
    if (out.total_pe == 0) {
        return out;
    }

    std::uniform_int_distribution<std::size_t> time_bin_dist(0, n_bins - 1);
    std::uniform_int_distribution<std::size_t> uniform_pixel_dist(
        0, pixel_ids.size() - 1);
    std::discrete_distribution<std::size_t> weighted_pixel_dist(
        out.rate_pe_per_ns_by_pixel.begin(), out.rate_pe_per_ns_by_pixel.end());
    std::map<std::size_t, std::uint64_t> pe_by_cell;
    for (std::uint64_t i = 0; i < out.total_pe; ++i) {
        const std::size_t col = uniform_rates
            ? uniform_pixel_dist(rng)
            : weighted_pixel_dist(rng);
        const std::size_t bin = time_bin_dist(rng);
        ++pe_by_cell[col * n_bins + bin];
    }

    out.samples.reserve(pe_by_cell.size());
    for (const auto& item : pe_by_cell) {
        const std::size_t col = item.first / n_bins;
        const std::size_t bin = item.first % n_bins;
        if (item.second > static_cast<std::uint64_t>(
                              std::numeric_limits<std::uint32_t>::max())) {
            throw std::runtime_error("NSB cell count exceeds supported range");
        }
        const float pe = static_cast<float>(item.second);
        out.samples.push_back({col, bin, pe});
        out.integrated_pe_by_pixel[col] += pe;
    }
    return out;
}

} // namespace lact
