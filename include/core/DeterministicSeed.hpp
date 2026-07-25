#pragma once

#include <cstdint>

namespace lact {

// Stable integer-only seed derivation. Keep these constants and operations
// versioned: unlike std::hash<double>, this has defined cross-platform bits.
inline std::uint64_t splitMix64(std::uint64_t value)
{
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31U);
}

inline std::uint64_t deriveSeed(std::uint64_t seed, std::uint64_t value)
{
    return splitMix64(seed ^ splitMix64(value));
}

enum class RandomStage : std::uint64_t {
    MissingWavelength = 1,
    PreGeometryAcceptance = 2,
    MirrorRoughness = 3,
    PostGeometryAcceptance = 4,
    FacetPositionError = 5,
    FacetNormalError = 6,
    FacetRadiusError = 7,
    FacetReflectivityError = 8,
};

inline std::uint64_t stageSeed(std::uint64_t seed, RandomStage stage)
{
    return deriveSeed(seed, static_cast<std::uint64_t>(stage));
}

} // namespace lact
