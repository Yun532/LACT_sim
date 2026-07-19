#include "io/EventIOPhotonSource.hpp"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

extern "C" {
#include "initial.h"
#include "io_basic.h"
#include "io_hess.h"
#include "mc_tel.h"
FILE* fileopen(const char* fname, const char* mode);
int fileclose(FILE* f);
int read_simtel_mc_shower(IO_BUFFER* iobuf, MCShower* mcs);
}

namespace {

constexpr double RAD_TO_DEG = 180.0 / 3.14159265358979323846;

struct IoBufferGuard {
    IO_BUFFER* ptr = nullptr;
    ~IoBufferGuard() {
        if (ptr) {
            free_io_buffer(ptr);
        }
    }
};

struct FileGuard {
    FILE* ptr = nullptr;
    ~FileGuard() {
        if (ptr) {
            fileclose(ptr);
        }
    }
};

double downwardDirZ(double cx, double cy) {
    double z2 = 1.0 - cx * cx - cy * cy;
    if (z2 < 0.0) {
        z2 = 0.0;
    }
    return -std::sqrt(z2);
}

std::string lowerCopy(std::string s) {
    for (auto& c : s) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return s;
}

int outputEventId(int shower_event_id, int array_id, const EventIOPhotonConfig& cfg) {
    const std::string mode = lowerCopy(cfg.event_id_mode);
    if (mode == "event_array100" || mode == "runid") {
        return shower_event_id * 100 + array_id;
    }
    if (mode == "event") {
        return shower_event_id;
    }
    throw std::runtime_error("unsupported source.event_id_mode: " + cfg.event_id_mode);
}

bool keepRow(int event_id, int shower_event_id, int telescope_id, const EventIOPhotonConfig& cfg) {
    if (cfg.filter_shower_event_id && shower_event_id != cfg.selected_shower_event_id) {
        return false;
    }
    if (cfg.filter_event_id && event_id != cfg.selected_event_id) {
        return false;
    }
    if (!cfg.selected_event_ids.empty() &&
        !std::binary_search(cfg.selected_event_ids.begin(),
                            cfg.selected_event_ids.end(), event_id)) {
        return false;
    }
    if (cfg.filter_telescope_id && telescope_id != cfg.selected_telescope_id) {
        return false;
    }
    return true;
}

std::uint64_t mixSeed(std::uint64_t seed, std::uint64_t value)
{
    seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
    return seed;
}

double unitRandom01(std::uint64_t seed)
{
    // SplitMix64 gives a deterministic pseudo-random value from bunch identity,
    // independent of streaming order and filtering.
    seed += 0x9e3779b97f4a7c15ULL;
    seed = (seed ^ (seed >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    seed = (seed ^ (seed >> 27U)) * 0x94d049bb133111ebULL;
    seed = seed ^ (seed >> 31U);
    constexpr double denom = 1.0 / 9007199254740992.0; // 2^53
    return static_cast<double>(seed >> 11U) * denom;
}

double sampleMissingWavelength(int shower_event_id,
                               int array_id,
                               int telescope_id,
                               std::uint64_t source_bunch_index,
                               std::uint64_t photon_index,
                               const EventIOPhotonConfig& cfg)
{
    const std::string model = lowerCopy(cfg.missing_wavelength_model);
    if (model.empty() || model == "default" || model == "fixed" ||
        model == "constant" || model == "none" || model == "off") {
        return cfg.default_wavelength_nm;
    }

    const double lo = cfg.missing_wavelength_min_nm;
    const double hi = cfg.missing_wavelength_max_nm;
    if (!std::isfinite(lo) || !std::isfinite(hi) || lo <= 0.0 || hi <= lo) {
        throw std::runtime_error(
            "source.missing_wavelength_min_nm/max_nm must be finite with 0 < min < max");
    }

    std::uint64_t seed = cfg.missing_wavelength_seed;
    seed = mixSeed(seed, static_cast<std::uint64_t>(shower_event_id + 1000003));
    seed = mixSeed(seed, static_cast<std::uint64_t>(array_id + 100003));
    seed = mixSeed(seed, static_cast<std::uint64_t>(telescope_id + 10007));
    seed = mixSeed(seed, source_bunch_index + 1ULL);
    seed = mixSeed(seed, photon_index + 1ULL);
    const double u = std::min(1.0 - 1e-16, std::max(1e-16, unitRandom01(seed)));

    if (model == "uniform") {
        return lo + u * (hi - lo);
    }
    if (model == "cherenkov" || model == "cherenkov_1_over_lambda2" ||
        model == "1_over_lambda2" || model == "one_over_lambda2") {
        const double inv_lo = 1.0 / lo;
        const double inv_hi = 1.0 / hi;
        return 1.0 / (inv_lo - u * (inv_lo - inv_hi));
    }

    throw std::runtime_error("unsupported source.missing_wavelength_model: " +
                             cfg.missing_wavelength_model);
}

using EventIdentity = std::pair<int, int>;

void registerOutputEventId(std::map<int, EventIdentity>& identities,
                           int event_id,
                           int shower_event_id,
                           int array_id,
                           const EventIOPhotonConfig& cfg)
{
    const std::string mode = lowerCopy(cfg.event_id_mode);
    if (mode != "event_array100" && mode != "runid") {
        return;
    }
    const EventIdentity identity{shower_event_id, array_id};
    const auto [it, inserted] = identities.emplace(event_id, identity);
    if (!inserted && it->second != identity) {
        throw std::runtime_error(
            "source.event_id_mode produced a collision for event_id=" +
            std::to_string(event_id) + ": (shower_event_id,array_id)=(" +
            std::to_string(it->second.first) + "," +
            std::to_string(it->second.second) + ") and (" +
            std::to_string(shower_event_id) + "," +
            std::to_string(array_id) + ")");
    }
}

void buildOutputEventIdentityMap(EventIOMetadata& metadata,
                                 const EventIOPhotonConfig& cfg)
{
    metadata.output_event_identity.clear();
    const std::string mode = lowerCopy(cfg.event_id_mode);
    if (mode != "event_array100" && mode != "runid") {
        return;
    }
    for (const auto& item : metadata.array_offsets_by_shower) {
        const int shower_event_id = item.first;
        const std::size_t n_arrays = std::max(item.second.x_m.size(),
                                              item.second.y_m.size());
        for (std::size_t array_index = 0; array_index < n_arrays; ++array_index) {
            const int array_id = static_cast<int>(array_index);
            const int event_id = outputEventId(shower_event_id, array_id, cfg);
            registerOutputEventId(metadata.output_event_identity,
                                  event_id,
                                  shower_event_id,
                                  array_id,
                                  cfg);
        }
    }
}

int selectedShowerEventId(const EventIOPhotonConfig& cfg) {
    if (!cfg.filter_event_id && cfg.selected_event_ids.empty()) {
        if (cfg.filter_shower_event_id) {
            return cfg.selected_shower_event_id;
        }
        return cfg.default_event_id;
    }
    const int selected_event = cfg.filter_event_id
                                   ? cfg.selected_event_id
                                   : cfg.selected_event_ids.front();
    const std::string mode = lowerCopy(cfg.event_id_mode);
    if (mode == "event_array100" || mode == "runid") {
        return selected_event / 100;
    }
    return selected_event;
}

int selectedArrayId(const EventIOPhotonConfig& cfg) {
    if (!cfg.filter_event_id && cfg.selected_event_ids.empty()) {
        return 0;
    }
    const int selected_event = cfg.filter_event_id
                                   ? cfg.selected_event_id
                                   : cfg.selected_event_ids.front();
    const std::string mode = lowerCopy(cfg.event_id_mode);
    if (mode == "event_array100" || mode == "runid") {
        return selected_event % 100;
    }
    return 0;
}

double interpolateAtmosphereThickness(const EventIOMetadata& metadata,
                                      double altitude_m)
{
    if (!std::isfinite(altitude_m) || metadata.atmosphere.size() < 2) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const auto& table = metadata.atmosphere;
    if (altitude_m <= table.front().altitude_m) {
        return table.front().thickness_g_cm2;
    }
    if (altitude_m >= table.back().altitude_m) {
        return table.back().thickness_g_cm2;
    }
    for (std::size_t i = 1; i < table.size(); ++i) {
        const auto& lo = table[i - 1];
        const auto& hi = table[i];
        if (altitude_m <= hi.altitude_m) {
            const double frac =
                (altitude_m - lo.altitude_m) / (hi.altitude_m - lo.altitude_m);
            return lo.thickness_g_cm2 +
                   frac * (hi.thickness_g_cm2 - lo.thickness_g_cm2);
        }
    }
    return std::numeric_limits<double>::quiet_NaN();
}

double interpolateAtmosphereHeight(const EventIOMetadata& metadata,
                                   double thickness_g_cm2)
{
    if (!std::isfinite(thickness_g_cm2) || metadata.atmosphere.size() < 2) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const auto& table = metadata.atmosphere;
    for (std::size_t i = 1; i < table.size(); ++i) {
        const auto& lower_alt = table[i - 1];
        const auto& upper_alt = table[i];
        const double thick_hi = lower_alt.thickness_g_cm2;
        const double thick_lo = upper_alt.thickness_g_cm2;
        if (thickness_g_cm2 <= thick_hi && thickness_g_cm2 >= thick_lo) {
            const double frac =
                (thickness_g_cm2 - thick_hi) / (thick_lo - thick_hi);
            return lower_alt.altitude_m +
                   frac * (upper_alt.altitude_m - lower_alt.altitude_m);
        }
    }
    if (thickness_g_cm2 > table.front().thickness_g_cm2) {
        return table.front().altitude_m;
    }
    return table.back().altitude_m;
}

void fillDerivedAtmosphericTruth(EventIOMetadata& metadata, EventIOEventHeader& event)
{
    if (!std::isfinite(event.starting_grammage_g_cm2) &&
        std::isfinite(event.h_first_int_m)) {
        event.starting_grammage_g_cm2 =
            interpolateAtmosphereThickness(metadata, event.h_first_int_m);
    }
    if (!std::isfinite(event.h_max_m) && std::isfinite(event.x_max_g_cm2)) {
        event.h_max_m = interpolateAtmosphereHeight(metadata, event.x_max_g_cm2);
    }
}

void fillDerivedAtmosphericTruth(EventIOMetadata& metadata)
{
    for (auto& event : metadata.events) {
        fillDerivedAtmosphericTruth(metadata, event);
    }
    if (metadata.selected_event) {
        auto it = std::find_if(
            metadata.events.begin(), metadata.events.end(),
            [&metadata](const EventIOEventHeader& event) {
                return event.shower_event_id ==
                       metadata.selected_event->shower_event_id;
            });
        if (it != metadata.events.end()) {
            metadata.selected_event = *it;
        }
    }
}

int readCurrentEventId(IO_BUFFER* iobuf, int item_type, int& current_event_id) {
    real data[273];
    int rc = read_tel_block(iobuf, item_type, data, 273);
    if (rc < 0) {
        return rc;
    }
    if (item_type == IO_TYPE_MC_EVTH) {
        current_event_id = static_cast<int>(Nint(data[1]));
    }
    return 0;
}

EventIOEventHeader& eventHeaderForShower(EventIOMetadata& metadata, int shower_event_id)
{
    auto it = std::find_if(
        metadata.events.begin(), metadata.events.end(),
        [shower_event_id](const EventIOEventHeader& event) {
            return event.shower_event_id == shower_event_id;
        });
    if (it != metadata.events.end()) {
        return *it;
    }
    metadata.events.push_back(EventIOEventHeader{});
    metadata.events.back().shower_event_id = shower_event_id;
    return metadata.events.back();
}

int readEventHeaderMetadata(IO_BUFFER* iobuf,
                            int item_type,
                            int selected_shower_event_id,
                            int& current_event_id,
                            EventIOMetadata& metadata)
{
    real data[273];
    int rc = read_tel_block(iobuf, item_type, data, 273);
    if (rc < 0) {
        return rc;
    }
    if (item_type == IO_TYPE_MC_EVTH) {
        current_event_id = static_cast<int>(Nint(data[1]));
        EventIOEventHeader& event = eventHeaderForShower(metadata, current_event_id);
        event.shower_event_id = current_event_id;
        event.primary_type = static_cast<int>(Nint(data[2]));
        event.energy_gev = data[3];
        event.theta_deg = data[10] * RAD_TO_DEG;
        event.phi_deg = data[11] * RAD_TO_DEG;
        event.altitude_deg = 90.0 - event.theta_deg;
        event.azimuth_north_to_east_deg = (data[92] - data[11] + M_PI) * RAD_TO_DEG;
        event.core_x_m = data[98] * 0.01;
        event.core_y_m = data[118] * 0.01;
        event.array_rotation_deg = data[92] * RAD_TO_DEG;
        if (std::isfinite(data[6]) && data[6] != 0.0) {
            event.h_first_int_m = std::fabs(data[6]) * 0.01;
        }
        fillDerivedAtmosphericTruth(metadata, event);
        if (current_event_id == selected_shower_event_id) {
            metadata.selected_event = event;
        }
    } else if (item_type == IO_TYPE_MC_EVTE) {
        current_event_id = static_cast<int>(Nint(data[1]));
        EventIOEventHeader& event = eventHeaderForShower(metadata, current_event_id);
        event.ground_gammas = data[2];
        event.ground_electrons = data[3];
        event.ground_hadrons = data[4];
        event.ground_muons = data[5];
        if (current_event_id == selected_shower_event_id) {
            metadata.selected_event = event;
        }
    }
    return 0;
}

int readLongitudinalMetadata(IO_BUFFER* iobuf,
                             int selected_shower_event_id,
                             EventIOMetadata& metadata)
{
    constexpr int kMaxProfiles = 16;
    constexpr int kMaxDepthBins = 10000;

    int event_id = 0;
    int profile_type = 0;
    int n_profiles = 0;
    int n_depth_bins = 0;
    double depth_step_g_cm2 = 0.0;
    std::vector<double> data(kMaxProfiles * kMaxDepthBins, 0.0);
    const int rc = read_shower_longitudinal(iobuf, &event_id, &profile_type,
                                            data.data(), kMaxDepthBins,
                                            &n_profiles, &n_depth_bins,
                                            &depth_step_g_cm2, kMaxProfiles);
    if (rc < 0) {
        return rc;
    }
    if (profile_type != 1 || n_profiles <= 0 || n_depth_bins <= 0 ||
        depth_step_g_cm2 <= 0.0) {
        return 0;
    }

    EventIOEventHeader& event = eventHeaderForShower(metadata, event_id);
    event.shower_event_id = event_id;

    int max_bin = -1;
    double max_value = -1.0;
    for (int bin = 0; bin < n_depth_bins; ++bin) {
        double value = 0.0;
        if (n_profiles > 2) {
            value = data[1 * kMaxDepthBins + bin] + data[2 * kMaxDepthBins + bin];
        } else {
            for (int profile = 0; profile < n_profiles; ++profile) {
                value += data[profile * kMaxDepthBins + bin];
            }
        }
        if (value > max_value) {
            max_value = value;
            max_bin = bin;
        }
    }
    if (max_bin >= 0 && max_value > 0.0) {
        event.x_max_g_cm2 = (static_cast<double>(max_bin) + 0.5) * depth_step_g_cm2;
        fillDerivedAtmosphericTruth(metadata, event);
        if (event_id == selected_shower_event_id) {
            metadata.selected_event = event;
        }
    }
    return 0;
}

int readAtmosphereMetadata(IO_BUFFER* iobuf, EventIOMetadata& metadata)
{
    AtmProf atmosphere{};
    const int rc = read_atmprof(iobuf, &atmosphere);
    if (rc < 0) {
        return rc;
    }

    if (std::isfinite(atmosphere.obslev)) {
        metadata.observation_altitude_m = atmosphere.obslev * 0.01;
    }

    metadata.atmosphere.clear();
    metadata.atmosphere.reserve(atmosphere.n_alt);
    for (unsigned i = 0; i < atmosphere.n_alt; ++i) {
        metadata.atmosphere.push_back(EventIOAtmosphereSample{
            atmosphere.alt_km[i] * 1000.0,
            atmosphere.thick[i],
        });
    }
    std::sort(metadata.atmosphere.begin(), metadata.atmosphere.end(),
              [](const auto& a, const auto& b) {
                  return a.altitude_m < b.altitude_m;
              });
    fillDerivedAtmosphericTruth(metadata);

    if (atmosphere.atmprof_fname != nullptr) free(atmosphere.atmprof_fname);
    if (atmosphere.alt_km != nullptr) free(atmosphere.alt_km);
    if (atmosphere.rho != nullptr) free(atmosphere.rho);
    if (atmosphere.thick != nullptr) free(atmosphere.thick);
    if (atmosphere.refidx_m1 != nullptr) free(atmosphere.refidx_m1);
    return 0;
}

void freeMcShowerProfileContent(MCShower& shower)
{
    const int n_profiles = std::min(shower.num_profiles, H_MAX_PROFILE);
    for (int i = 0; i < n_profiles; ++i) {
        if (shower.profile[i].content != nullptr) {
            free(shower.profile[i].content);
            shower.profile[i].content = nullptr;
            shower.profile[i].max_steps = 0;
        }
    }
}

int readSimtelMcShowerMetadata(IO_BUFFER* iobuf,
                               int selected_shower_event_id,
                               EventIOMetadata& metadata)
{
    MCShower shower{};
    int rc = read_simtel_mc_shower(iobuf, &shower);
    if (rc < 0) {
        freeMcShowerProfileContent(shower);
        return rc;
    }

    EventIOEventHeader& event = eventHeaderForShower(metadata, shower.shower_num);
    event.shower_event_id = shower.shower_num;
    event.primary_type = shower.primary_id;
    event.energy_gev = shower.energy * 1000.0;
    event.altitude_deg = shower.altitude * RAD_TO_DEG;
    event.theta_deg = 90.0 - event.altitude_deg;
    event.azimuth_north_to_east_deg = shower.azimuth * RAD_TO_DEG;
    event.h_first_int_m = shower.h_first_int;
    event.x_max_g_cm2 = shower.xmax;
    event.h_max_m = shower.hmax;
    event.starting_grammage_g_cm2 = shower.depth_start;
    event.has_simtel_mc_shower = true;
    if (shower.shower_num == selected_shower_event_id) {
        metadata.selected_event = event;
    }

    freeMcShowerProfileContent(shower);
    return 0;
}

void freeLinkedStrings(struct linked_string& list) {
    struct linked_string* xl = &list;
    while (xl != nullptr) {
        struct linked_string* next = xl->next;
        if (xl->text != nullptr) {
            free(xl->text);
            xl->text = nullptr;
        }
        if (xl != &list) {
            free(xl);
        }
        xl = next;
    }
    list.next = nullptr;
}

int readInputLinesMetadata(IO_BUFFER* iobuf, EventIOMetadata& metadata) {
    struct linked_string lines;
    lines.text = nullptr;
    lines.next = nullptr;
    int rc = read_input_lines(iobuf, &lines);
    if (rc < 0) {
        return rc;
    }
    for (struct linked_string* xl = &lines; xl != nullptr; xl = xl->next) {
        if (xl->text != nullptr) {
            metadata.input_lines.emplace_back(xl->text);
        }
    }
    freeLinkedStrings(lines);
    return 0;
}

int readTelescopePositionsMetadata(IO_BUFFER* iobuf, EventIOMetadata& metadata) {
    constexpr int MAX_TEL = 10000;
    int ntel = 0;
    std::vector<double> x(MAX_TEL), y(MAX_TEL), z(MAX_TEL), r(MAX_TEL);
    int rc = read_tel_pos(iobuf, MAX_TEL, &ntel, x.data(), y.data(), z.data(), r.data());
    if (rc < 0) {
        return rc;
    }
    metadata.telescopes.clear();
    for (int i = 0; i < ntel; ++i) {
        EventIOTelescopePosition tel;
        tel.telescope_id = i;
        tel.x_m = x[static_cast<std::size_t>(i)] * 0.01;
        tel.y_m = y[static_cast<std::size_t>(i)] * 0.01;
        tel.z_m = z[static_cast<std::size_t>(i)] * 0.01;
        tel.radius_m = r[static_cast<std::size_t>(i)] * 0.01;
        metadata.telescopes.push_back(tel);
    }
    return 0;
}

int readArrayOffsetsMetadata(IO_BUFFER* iobuf,
                             int current_event_id,
                             int selected_shower_event_id,
                             EventIOMetadata& metadata)
{
    constexpr int MAX_ARRAY = 10000;
    int narray = 0;
    double toff = 0.0;
    std::vector<double> x(MAX_ARRAY), y(MAX_ARRAY), w(MAX_ARRAY);
    int rc = read_tel_offset_w(iobuf, MAX_ARRAY, &narray, &toff, x.data(), y.data(), w.data());
    if (rc < 0) {
        return rc;
    }
    EventIOArrayOffsets offsets;
    offsets.time_offset_ns = toff;
    offsets.x_m.reserve(static_cast<std::size_t>(narray));
    offsets.y_m.reserve(static_cast<std::size_t>(narray));
    offsets.weight.reserve(static_cast<std::size_t>(narray));
    for (int i = 0; i < narray; ++i) {
        offsets.x_m.push_back(x[static_cast<std::size_t>(i)] * 0.01);
        offsets.y_m.push_back(y[static_cast<std::size_t>(i)] * 0.01);
        // MC_TELOFF follows CORSIKA's centimetre coordinate convention.  This
        // is an event/core-sampling area weight, not a per-photon efficiency.
        offsets.weight.push_back(w[static_cast<std::size_t>(i)] * 1.0e-4);
        offsets.has_explicit_weights = offsets.has_explicit_weights ||
                                       w[static_cast<std::size_t>(i)] != 0.0;
    }
    metadata.array_offsets_by_shower[current_event_id] = offsets;
    if (current_event_id == selected_shower_event_id) {
        metadata.selected_event_offsets = offsets;
    }
    return 0;
}

PhotonBunch makeBunch(const struct bunch& b,
                      int shower_event_id,
                      int array_id,
                      int event_id,
                      int telescope_id,
                      std::uint64_t source_bunch_index,
                      const EventIOPhotonConfig& cfg)
{
    PhotonBunch out;
    out.photon.pos = {b.x * 0.01, b.y * 0.01, 0.0};
    out.photon.dir = {b.cx, b.cy, downwardDirZ(b.cx, b.cy)};
    out.photon.normalizeDirection();
    out.photon.time_ns = b.ctime;
    out.photon.weight = cfg.default_weight;
    out.photon.optical_efficiency_preapplied = b.lambda < 0.0;
    out.multiplicity = b.photons * cfg.default_multiplicity;
    out.event_id = event_id;
    out.shower_event_id = shower_event_id;
    out.array_id = array_id;
    out.telescope_id = telescope_id;
    out.source_bunch_index = source_bunch_index;
    out.raw_wavelength_nm = b.lambda;
    out.photon.wavelength_nm = resolveEventIOPhotonWavelength(out, 0, cfg);
    out.eventio_2d = true;
    if (std::isfinite(b.zem) && b.zem > 0.0) {
        out.emission_altitude_km = b.zem * 1.0e-5;
    }
    return out;
}

bool isEmitterRecord(const struct bunch& b) {
    return b.lambda >= 9000.0;
}

bool isEmitterRecord(const struct bunch3d& b) {
    return b.lambda >= 9000.0;
}

void attachEmitter(PhotonBunch& out, const struct bunch& emitter) {
    out.has_emitter = true;
    out.emitter_mass_gev = emitter.cx;
    out.emitter_charge = emitter.cy;
    out.emitter_energy_gev = emitter.photons;
    out.emitter_time_ns = emitter.zem;
}

void attachEmitter(PhotonBunch& out, const struct bunch3d& emitter) {
    out.has_emitter = true;
    out.emitter_mass_gev = emitter.cx;
    out.emitter_charge = emitter.cy;
    out.emitter_energy_gev = emitter.photons;
    out.emitter_time_ns = emitter.dist;
}

PhotonBunch makeBunch3d(const struct bunch3d& b,
                        int shower_event_id,
                        int array_id,
                        int event_id,
                        int telescope_id,
                        std::uint64_t source_bunch_index,
                        const EventIOPhotonConfig& cfg)
{
    PhotonBunch out;
    out.photon.pos = {b.x * 0.01, b.y * 0.01, b.z * 0.01};
    out.photon.dir = {b.cx, b.cy, b.cz};
    out.photon.normalizeDirection();
    out.photon.time_ns = b.ctime;
    out.photon.weight = cfg.default_weight;
    out.photon.optical_efficiency_preapplied = b.lambda < 0.0;
    out.multiplicity = b.photons * cfg.default_multiplicity;
    out.event_id = event_id;
    out.shower_event_id = shower_event_id;
    out.array_id = array_id;
    out.telescope_id = telescope_id;
    out.source_bunch_index = source_bunch_index;
    out.raw_wavelength_nm = b.lambda;
    out.photon.wavelength_nm = resolveEventIOPhotonWavelength(out, 0, cfg);
    out.eventio_2d = false;
    out.emission_altitude_km = eventIO3DEmissionAltitudeKm(
        cfg.observation_altitude_km, b.z, b.cz, b.dist);
    return out;
}

int readPhotonBlock(IO_BUFFER* iobuf,
                    int current_event_id,
                    const EventIOPhotonConfig& cfg,
                    const EventIOPhotonCallback& on_bunch,
                    std::map<int, EventIdentity>& event_identities,
                    std::size_t& emitted,
                    std::size_t& emitted_2d)
{
    int array_id = 0;
    int telescope_id = 0;
    int nbunches = 0;
    double photons = 0.0;
    IO_ITEM_HEADER item_header;
    item_header.type = IO_TYPE_MC_PHOTONS;
    int rc = get_item_begin(iobuf, &item_header);
    if (rc < 0) {
        return rc;
    }

    const int version_group = static_cast<int>(item_header.version / 1000);
    if (version_group != 0 && version_group != 1) {
        get_item_end(iobuf, &item_header);
        return -1;
    }

    array_id = get_short(iobuf);
    telescope_id = get_short(iobuf);
    photons = get_real(iobuf);
    nbunches = get_long(iobuf);

    // CORSIKA may store non-photon ground particles with the same block type.
    // They are useful for shower diagnostics but not for optical tracing.
    const bool particle_block = (array_id == 999 && telescope_id == 999);
    if (particle_block) {
        return get_item_end(iobuf, &item_header);
    }

    const int event_id = outputEventId(current_event_id, array_id, cfg);
    if (!keepRow(event_id, current_event_id, telescope_id, cfg)) {
        return get_item_end(iobuf, &item_header);
    }
    registerOutputEventId(event_identities, event_id, current_event_id, array_id, cfg);

    bool have_pending = false;
    struct bunch pending{};
    std::uint64_t pending_index = 0;
    auto emit_pending = [&]() {
        if (!have_pending) {
            return;
        }
        on_bunch(makeBunch(pending, current_event_id, array_id, event_id,
                           telescope_id, pending_index, cfg));
        ++emitted;
        ++emitted_2d;
        have_pending = false;
    };
    auto emit_with_emitter = [&](const struct bunch& emitter) {
        PhotonBunch out = makeBunch(pending, current_event_id, array_id, event_id,
                                    telescope_id, pending_index, cfg);
        attachEmitter(out, emitter);
        on_bunch(out);
        ++emitted;
        ++emitted_2d;
        have_pending = false;
    };

    for (int i = 0; i < nbunches; ++i) {
        struct bunch b;
        if (version_group == 0) {
            b.x = get_real(iobuf);
            b.y = get_real(iobuf);
            b.cx = get_real(iobuf);
            b.cy = get_real(iobuf);
            b.ctime = get_real(iobuf);
            b.zem = get_real(iobuf);
            b.photons = get_real(iobuf);
            b.lambda = get_real(iobuf);
        } else {
            b.x = 0.1 * get_short(iobuf);
            b.y = 0.1 * get_short(iobuf);
            b.cx = get_short(iobuf) / 30000.0;
            if (b.cx > 1.0) b.cx = 1.0;
            if (b.cx < -1.0) b.cx = -1.0;
            b.cy = get_short(iobuf) / 30000.0;
            if (b.cy > 1.0) b.cy = 1.0;
            if (b.cy < -1.0) b.cy = -1.0;
            b.ctime = 0.1 * get_short(iobuf);
            b.zem = std::pow(10.0, 0.001 * get_short(iobuf));
            b.photons = 0.01 * get_short(iobuf);
            b.lambda = get_short(iobuf);
        }

        if (isEmitterRecord(b)) {
            if (have_pending) {
                if (cfg.read_emitter_info) {
                    emit_with_emitter(b);
                } else {
                    emit_pending();
                }
            }
            continue;
        }

        emit_pending();
        pending = b;
        pending_index = static_cast<std::uint64_t>(i);
        have_pending = true;
    }
    emit_pending();

    return get_item_end(iobuf, &item_header);
}

int readPhoton3dBlock(IO_BUFFER* iobuf,
                      int current_event_id,
                      const EventIOPhotonConfig& cfg,
                      const EventIOPhotonCallback& on_bunch,
                      std::map<int, EventIdentity>& event_identities,
                      std::size_t& emitted,
                      std::size_t& emitted_3d)
{
    int array_id = 0;
    int telescope_id = 0;
    int nbunches = 0;
    double photons = 0.0;
    int rc = read_tel_photons3d(iobuf, 0, &array_id, &telescope_id, &photons,
                                nullptr, &nbunches);
    if (rc != -10) {
        return rc;
    }

    std::vector<struct bunch3d> bunches(static_cast<std::size_t>(nbunches));
    rc = read_tel_photons3d(iobuf, nbunches, &array_id, &telescope_id, &photons,
                            bunches.data(), &nbunches);
    if (rc < 0) {
        return rc;
    }

    const int event_id = outputEventId(current_event_id, array_id, cfg);
    if (!keepRow(event_id, current_event_id, telescope_id, cfg)) {
        return 0;
    }
    registerOutputEventId(event_identities, event_id, current_event_id, array_id, cfg);
    bool have_pending = false;
    struct bunch3d pending{};
    std::uint64_t pending_index = 0;
    auto emit_pending = [&]() {
        if (!have_pending) {
            return;
        }
        on_bunch(makeBunch3d(pending, current_event_id, array_id, event_id,
                             telescope_id, pending_index, cfg));
        ++emitted;
        ++emitted_3d;
        have_pending = false;
    };
    auto emit_with_emitter = [&](const struct bunch3d& emitter) {
        PhotonBunch out = makeBunch3d(pending, current_event_id, array_id,
                                      event_id, telescope_id, pending_index, cfg);
        attachEmitter(out, emitter);
        on_bunch(out);
        ++emitted;
        ++emitted_3d;
        have_pending = false;
    };
    for (int i = 0; i < nbunches; ++i) {
        const auto& b = bunches[static_cast<std::size_t>(i)];
        if (isEmitterRecord(b)) {
            if (have_pending) {
                if (cfg.read_emitter_info) {
                    emit_with_emitter(b);
                } else {
                    emit_pending();
                }
            }
            continue;
        }
        emit_pending();
        pending = b;
        pending_index = static_cast<std::uint64_t>(i);
        have_pending = true;
    }
    emit_pending();
    return 0;
}

int readTelArray(IO_BUFFER* iobuf,
                 int current_event_id,
                 const EventIOPhotonConfig& cfg,
                 const EventIOPhotonCallback& on_bunch,
                 std::map<int, EventIdentity>& event_identities,
                 EventIOStreamStats& stats)
{
    IO_ITEM_HEADER array_header;
    int array_id = 0;
    int rc = begin_read_tel_array(iobuf, &array_header, &array_id);
    if (rc < 0) {
        return rc;
    }
    int type = 0;
    while ((type = next_subitem_type(iobuf)) > 0) {
        if (type == IO_TYPE_MC_PHOTONS) {
            rc = readPhotonBlock(iobuf, current_event_id, cfg, on_bunch,
                                 event_identities,
                                 stats.photon_bunches, stats.photon_bunches_2d);
        } else if (type == IO_TYPE_MC_PHOTONS3D) {
            rc = readPhoton3dBlock(iobuf, current_event_id, cfg, on_bunch,
                                   event_identities,
                                   stats.photon_bunches, stats.photon_bunches_3d);
        } else {
            rc = skip_subitem(iobuf);
        }
        if (rc < 0) {
            get_item_end(iobuf, &array_header);
            return rc;
        }
    }
    return end_read_tel_array(iobuf, &array_header);
}

void printMetadataBrief(const EventIOMetadata& metadata) {
    std::cerr << "EventIOPhotonSource: metadata shower_events="
              << metadata.events.size()
              << " telescopes=" << metadata.telescopes.size()
              << " input_card_lines=" << metadata.input_lines.size();
    if (!metadata.events.empty()) {
        int min_event = std::numeric_limits<int>::max();
        int max_event = std::numeric_limits<int>::min();
        for (const auto& event : metadata.events) {
            min_event = std::min(min_event, event.shower_event_id);
            max_event = std::max(max_event, event.shower_event_id);
        }
        std::cerr << " shower_event_range=" << min_event << ".." << max_event;
    }
    std::cerr << "\n";
}

} // namespace

double resolveEventIOPhotonWavelength(const PhotonBunch& bunch,
                                      std::uint64_t photon_index,
                                      const EventIOPhotonConfig& cfg)
{
    if (bunch.raw_wavelength_nm > 0.0) {
        return bunch.raw_wavelength_nm;
    }
    if (bunch.raw_wavelength_nm < 0.0) {
        // CEFFIC bunches no longer carry a physical wavelength.  The
        // placeholder is used only by wavelength-agnostic geometry code.
        return cfg.default_wavelength_nm;
    }
    return sampleMissingWavelength(bunch.shower_event_id,
                                   bunch.array_id,
                                   bunch.telescope_id,
                                   bunch.source_bunch_index,
                                   photon_index,
                                   cfg);
}

double eventIO3DEmissionAltitudeKm(double observation_altitude_km,
                                   double bunch_z_cm,
                                   double direction_cosine_z,
                                   double emission_distance_cm)
{
    if (!std::isfinite(observation_altitude_km) ||
        !std::isfinite(bunch_z_cm) ||
        !std::isfinite(direction_cosine_z) ||
        !std::isfinite(emission_distance_cm)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return observation_altitude_km +
           (bunch_z_cm - direction_cosine_z * emission_distance_cm) * 1.0e-5;
}

std::optional<EventIOTelescopePosition>
EventIOMetadata::telescopeById(int telescope_id) const {
    for (const auto& tel : telescopes) {
        if (tel.telescope_id == telescope_id) {
            return tel;
        }
    }
    return std::nullopt;
}

std::optional<EventIOArrayOffsets>
EventIOMetadata::arrayOffsetsForShower(int shower_event_id) const {
    const auto it = array_offsets_by_shower.find(shower_event_id);
    if (it == array_offsets_by_shower.end()) {
        return std::nullopt;
    }
    return it->second;
}

EventIOMetadata readEventIOMetadata(const EventIOPhotonConfig& cfg) {
    if (cfg.path.empty()) {
        throw std::runtime_error("readEventIOMetadata: source.eventio_path is required");
    }

    IoBufferGuard iobuf;
    iobuf.ptr = allocate_io_buffer(5000000L);
    if (!iobuf.ptr) {
        throw std::runtime_error("readEventIOMetadata: failed to allocate IO buffer");
    }
    if (iobuf.ptr->max_length < 1000000000L) {
        iobuf.ptr->max_length = 1000000000L;
    }

    FileGuard input;
    input.ptr = fileopen(cfg.path.c_str(), READ_BINARY);
    if (!input.ptr) {
        throw std::runtime_error("readEventIOMetadata: failed to open " + cfg.path);
    }
    iobuf.ptr->input_file = input.ptr;

    EventIOMetadata metadata;
    metadata.selected_array_id = selectedArrayId(cfg);
    const int selected_shower_event_id = selectedShowerEventId(cfg);
    int current_event_id = cfg.default_event_id;
    IO_ITEM_HEADER item_header;
    while (find_io_block(iobuf.ptr, &item_header) == 0) {
        if (read_io_block(iobuf.ptr, &item_header) != 0) {
            break;
        }
        int rc = 0;
        switch (static_cast<int>(item_header.type)) {
            case IO_TYPE_MC_INPUTCFG:
                rc = readInputLinesMetadata(iobuf.ptr, metadata);
                break;
            case IO_TYPE_MC_ATMPROF:
                rc = readAtmosphereMetadata(iobuf.ptr, metadata);
                break;
            case IO_TYPE_MC_TELPOS:
                rc = readTelescopePositionsMetadata(iobuf.ptr, metadata);
                break;
            case IO_TYPE_MC_EVTH:
            case IO_TYPE_MC_RUNH:
            case IO_TYPE_MC_EVTE:
            case IO_TYPE_MC_RUNE:
                rc = readEventHeaderMetadata(iobuf.ptr, static_cast<int>(item_header.type),
                                             selected_shower_event_id, current_event_id,
                                             metadata);
                break;
            case IO_TYPE_SIMTEL_MC_SHOWER:
                rc = readSimtelMcShowerMetadata(iobuf.ptr, selected_shower_event_id,
                                                metadata);
                break;
            case IO_TYPE_MC_LONGI:
                rc = readLongitudinalMetadata(iobuf.ptr, selected_shower_event_id,
                                              metadata);
                break;
            case IO_TYPE_MC_TELOFF:
                rc = readArrayOffsetsMetadata(iobuf.ptr, current_event_id,
                                              selected_shower_event_id, metadata);
                break;
            default:
                rc = 0;
                break;
        }
        if (rc < 0) {
            throw std::runtime_error("readEventIOMetadata: failed while reading block type " +
                                     std::to_string(item_header.type));
        }
    }
    iobuf.ptr->input_file = nullptr;
    buildOutputEventIdentityMap(metadata, cfg);
    return metadata;
}

EventIOStreamStats streamEventIOPhotonBunches(
    const EventIOPhotonConfig& cfg,
    const EventIOPhotonCallback& on_bunch,
    const EventIOProgressCallback& on_progress)
{
    if (cfg.path.empty()) {
        throw std::runtime_error("streamEventIOPhotonBunches: source.eventio_path is required");
    }
    if (!on_bunch) {
        throw std::runtime_error("streamEventIOPhotonBunches: on_bunch callback is required");
    }

    IoBufferGuard iobuf;
    iobuf.ptr = allocate_io_buffer(5000000L);
    if (!iobuf.ptr) {
        throw std::runtime_error("streamEventIOPhotonBunches: failed to allocate IO buffer");
    }
    if (iobuf.ptr->max_length < 1000000000L) {
        iobuf.ptr->max_length = 1000000000L;
    }

    FileGuard input;
    input.ptr = fileopen(cfg.path.c_str(), READ_BINARY);
    if (!input.ptr) {
        throw std::runtime_error("streamEventIOPhotonBunches: failed to open " + cfg.path);
    }
    iobuf.ptr->input_file = input.ptr;

    int current_event_id = cfg.default_event_id;
    IO_ITEM_HEADER item_header;
    EventIOStreamStats stats;
    std::set<int> streamed_shower_events;
    std::map<int, EventIdentity> event_identities;
    bool stop_after_current_block = false;
    std::size_t next_report_rows = 1000000;
    const auto load_start = std::chrono::steady_clock::now();
    while (find_io_block(iobuf.ptr, &item_header) == 0) {
        if (stop_after_current_block) {
            break;
        }
        if (read_io_block(iobuf.ptr, &item_header) != 0) {
            break;
        }

        int rc = 0;
        switch (static_cast<int>(item_header.type)) {
            case IO_TYPE_MC_EVTH:
            case IO_TYPE_MC_RUNH:
            case IO_TYPE_MC_EVTE:
            case IO_TYPE_MC_RUNE:
                rc = readCurrentEventId(iobuf.ptr, static_cast<int>(item_header.type),
                                        current_event_id);
                if (rc >= 0 && static_cast<int>(item_header.type) == IO_TYPE_MC_EVTH) {
                    if (!cfg.filter_shower_event_id ||
                        current_event_id == cfg.selected_shower_event_id) {
                        streamed_shower_events.insert(current_event_id);
                        if (cfg.max_shower_events > 0 &&
                            static_cast<int>(streamed_shower_events.size()) >
                                cfg.max_shower_events) {
                            stop_after_current_block = true;
                        }
                    }
                }
                break;
            case IO_TYPE_MC_TELARRAY:
                if (stop_after_current_block) {
                    rc = 0;
                    break;
                }
                rc = readTelArray(iobuf.ptr, current_event_id, cfg, on_bunch,
                                  event_identities,
                                  stats);
                break;
            case IO_TYPE_MC_PHOTONS:
                if (stop_after_current_block) {
                    rc = 0;
                    break;
                }
                rc = readPhotonBlock(iobuf.ptr, current_event_id, cfg, on_bunch,
                                     event_identities,
                                     stats.photon_bunches, stats.photon_bunches_2d);
                break;
            case IO_TYPE_MC_PHOTONS3D:
                if (stop_after_current_block) {
                    rc = 0;
                    break;
                }
                rc = readPhoton3dBlock(iobuf.ptr, current_event_id, cfg, on_bunch,
                                       event_identities,
                                       stats.photon_bunches, stats.photon_bunches_3d);
                break;
            default:
                rc = 0;
                break;
        }
        if (rc < 0) {
            throw std::runtime_error("streamEventIOPhotonBunches: failed while reading block type " +
                                     std::to_string(item_header.type));
        }
        if (on_progress && stats.photon_bunches >= next_report_rows) {
            const auto now = std::chrono::steady_clock::now();
            const double elapsed_s =
                std::chrono::duration<double>(now - load_start).count();
            on_progress({stats.photon_bunches,
                         stats.photon_bunches_2d,
                         stats.photon_bunches_3d,
                         current_event_id,
                         elapsed_s,
                         false});
            next_report_rows += 1000000;
        }
    }

    iobuf.ptr->input_file = nullptr;
    if (stats.photon_bunches == 0) {
        throw std::runtime_error("streamEventIOPhotonBunches: no photon bunches matched filters in " +
                                 cfg.path);
    }
    if (on_progress) {
        const auto load_done = std::chrono::steady_clock::now();
        const double elapsed_s =
            std::chrono::duration<double>(load_done - load_start).count();
        on_progress({stats.photon_bunches,
                     stats.photon_bunches_2d,
                     stats.photon_bunches_3d,
                     current_event_id,
                     elapsed_s,
                     true});
    }
    return stats;
}

EventIOPhotonSource::EventIOPhotonSource(const EventIOPhotonConfig& cfg)
    : cfg_(cfg)
{
    if (cfg_.path.empty()) {
        throw std::runtime_error("EventIOPhotonSource: source.eventio_path is required");
    }
    load();
}

void EventIOPhotonSource::reset() {
    index_ = 0;
}

bool EventIOPhotonSource::next(PhotonBunch& out) {
    if (index_ >= rows_.size()) {
        return false;
    }
    out = rows_[index_++];
    return true;
}

void EventIOPhotonSource::load() {
    std::cerr << "EventIOPhotonSource: reading metadata from " << cfg_.path << "\n";
    metadata_ = readEventIOMetadata(cfg_);
    printMetadataBrief(metadata_);
    std::cerr << "EventIOPhotonSource: loading photon bunches "
              << "(preloading all matching bunches before tracing)\n";
    auto stats = streamEventIOPhotonBunches(
        cfg_,
        [this](const PhotonBunch& bunch) {
            rows_.push_back(bunch);
        },
        [](const EventIOStreamProgress& progress) {
            std::cerr << "EventIOPhotonSource: "
                      << (progress.final ? "loaded total photon_bunches="
                                         : "loaded photon_bunches=")
                      << progress.photon_bunches
                      << " current_shower_event=" << progress.current_shower_event
                      << " elapsed_s=" << progress.elapsed_s << "\n";
        });
    (void)stats;
}
