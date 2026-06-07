#pragma once

#include <cstddef>
#include <functional>
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
};

struct EventIOEventHeader {
    int shower_event_id = 0;
    int primary_type = 0;
    double energy_gev = 0.0;
    double theta_deg = 0.0;
    double phi_deg = 0.0;
    double azimuth_north_to_east_deg = 0.0;
    double core_x_m = 0.0;
    double core_y_m = 0.0;
    double array_rotation_deg = 0.0;
};

struct EventIOMetadata {
    std::vector<std::string> input_lines;
    std::vector<EventIOTelescopePosition> telescopes;
    std::vector<EventIOEventHeader> events;
    std::optional<EventIOEventHeader> selected_event;
    std::optional<EventIOArrayOffsets> selected_event_offsets;
    std::map<int, EventIOArrayOffsets> array_offsets_by_shower;
    int selected_array_id = 0;

    std::optional<EventIOTelescopePosition> telescopeById(int telescope_id) const;
    std::optional<EventIOArrayOffsets> arrayOffsetsForShower(int shower_event_id) const;
};

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
