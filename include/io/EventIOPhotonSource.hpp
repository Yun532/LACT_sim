#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "io/PhotonSource.hpp"

struct EventIOPhotonConfig {
    std::string path;
    bool local_telescope_frame = true;
    std::string event_id_mode = "event";
    double default_wavelength_nm = 400.0;
    std::string missing_wavelength_model = "cherenkov";
    double missing_wavelength_min_nm = 260.0;
    double missing_wavelength_max_nm = 1000.0;
    std::uint64_t missing_wavelength_seed = 246813579ULL;
    double observation_altitude_km = std::numeric_limits<double>::quiet_NaN();
    double default_time_ns = 0.0;
    double default_weight = 1.0;
    double default_multiplicity = 1.0;
    int default_event_id = 0;
    int default_telescope_id = 0;
    bool filter_telescope_id = false;
    int selected_telescope_id = 0;
    bool filter_event_id = false;
    int selected_event_id = 0;
    bool filter_shower_event_id = false;
    int selected_shower_event_id = 0;
    int max_shower_events = 0;
    bool read_emitter_info = false;
};

struct EventIOTelescopePosition {
    int telescope_id = 0;
    double x_m = 0.0;
    double y_m = 0.0;
    double z_m = 0.0;
    double radius_m = 0.0;
};

struct EventIOArrayOffsets {
    double time_offset_ns = 0.0;
    std::vector<double> x_m;
    std::vector<double> y_m;
    std::vector<double> weight;
    bool has_explicit_weights = false;
};

struct EventIOEventHeader {
    int shower_event_id = 0;
    int primary_type = 0;
    double energy_gev = 0.0;
    double theta_deg = 0.0;
    double phi_deg = 0.0;
    double altitude_deg = std::numeric_limits<double>::quiet_NaN();
    double azimuth_north_to_east_deg = 0.0;
    double core_x_m = 0.0;
    double core_y_m = 0.0;
    double array_rotation_deg = 0.0;
    double h_first_int_m = std::numeric_limits<double>::quiet_NaN();
    double x_max_g_cm2 = std::numeric_limits<double>::quiet_NaN();
    double h_max_m = std::numeric_limits<double>::quiet_NaN();
    double starting_grammage_g_cm2 = std::numeric_limits<double>::quiet_NaN();
    double ground_gammas = std::numeric_limits<double>::quiet_NaN();
    double ground_electrons = std::numeric_limits<double>::quiet_NaN();
    double ground_hadrons = std::numeric_limits<double>::quiet_NaN();
    double ground_muons = std::numeric_limits<double>::quiet_NaN();
    bool has_simtel_mc_shower = false;
};

struct EventIOAtmosphereSample {
    double altitude_m = 0.0;
    double thickness_g_cm2 = 0.0;
};

struct EventIOMetadata {
    std::vector<std::string> input_lines;
    std::vector<EventIOTelescopePosition> telescopes;
    std::vector<EventIOEventHeader> events;
    std::vector<EventIOAtmosphereSample> atmosphere;
    std::optional<EventIOEventHeader> selected_event;
    std::optional<EventIOArrayOffsets> selected_event_offsets;
    std::map<int, EventIOArrayOffsets> array_offsets_by_shower;
    std::map<int, std::pair<int, int>> output_event_identity;
    int selected_array_id = 0;
    double observation_altitude_m = std::numeric_limits<double>::quiet_NaN();

    std::optional<EventIOTelescopePosition> telescopeById(int telescope_id) const;
    std::optional<EventIOArrayOffsets> arrayOffsetsForShower(int shower_event_id) const;
};

double resolveEventIOPhotonWavelength(const PhotonBunch& bunch,
                                      std::uint64_t photon_index,
                                      const EventIOPhotonConfig& cfg);

double eventIO3DEmissionAltitudeKm(double observation_altitude_km,
                                   double bunch_z_cm,
                                   double direction_cosine_z,
                                   double emission_distance_cm);

EventIOMetadata readEventIOMetadata(const EventIOPhotonConfig& cfg);

struct EventIOStreamProgress {
    std::size_t photon_bunches = 0;
    std::size_t photon_bunches_2d = 0;
    std::size_t photon_bunches_3d = 0;
    int current_shower_event = 0;
    double elapsed_s = 0.0;
    bool final = false;
};

struct EventIOStreamStats {
    std::size_t photon_bunches = 0;
    std::size_t photon_bunches_2d = 0;
    std::size_t photon_bunches_3d = 0;
};

using EventIOPhotonCallback = std::function<void(const PhotonBunch&)>;
using EventIOProgressCallback = std::function<void(const EventIOStreamProgress&)>;

EventIOStreamStats streamEventIOPhotonBunches(
    const EventIOPhotonConfig& cfg,
    const EventIOPhotonCallback& on_bunch,
    const EventIOProgressCallback& on_progress = EventIOProgressCallback{});

class EventIOPhotonSource : public PhotonSource {
public:
    explicit EventIOPhotonSource(const EventIOPhotonConfig& cfg);

    bool next(PhotonBunch& out) override;
    void reset() override;

    std::size_t size() const { return rows_.size(); }
    bool localTelescopeFrame() const { return cfg_.local_telescope_frame; }
    const EventIOMetadata& metadata() const { return metadata_; }

private:
    EventIOPhotonConfig cfg_;
    std::size_t index_ = 0;
    std::vector<PhotonBunch> rows_;
    EventIOMetadata metadata_;

    void load();
};
