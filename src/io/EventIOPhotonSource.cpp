#include "io/EventIOPhotonSource.hpp"

#include <cmath>
#include <chrono>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

extern "C" {
#include "initial.h"
#include "io_basic.h"
#include "mc_tel.h"
FILE* fileopen(const char* fname, const char* mode);
int fileclose(FILE* f);
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
    if (cfg.filter_telescope_id && telescope_id != cfg.selected_telescope_id) {
        return false;
    }
    return true;
}

int selectedShowerEventId(const EventIOPhotonConfig& cfg) {
    if (!cfg.filter_event_id) {
        if (cfg.filter_shower_event_id) {
            return cfg.selected_shower_event_id;
        }
        return cfg.default_event_id;
    }
    const std::string mode = lowerCopy(cfg.event_id_mode);
    if (mode == "event_array100" || mode == "runid") {
        return cfg.selected_event_id / 100;
    }
    return cfg.selected_event_id;
}

int selectedArrayId(const EventIOPhotonConfig& cfg) {
    if (!cfg.filter_event_id) {
        return 0;
    }
    const std::string mode = lowerCopy(cfg.event_id_mode);
    if (mode == "event_array100" || mode == "runid") {
        return cfg.selected_event_id % 100;
    }
    return 0;
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
        EventIOEventHeader event;
        event.shower_event_id = current_event_id;
        event.primary_type = static_cast<int>(Nint(data[2]));
        event.energy_gev = data[3];
        event.theta_deg = data[10] * RAD_TO_DEG;
        event.phi_deg = data[11] * RAD_TO_DEG;
        event.azimuth_north_to_east_deg = (data[92] - data[11] + M_PI) * RAD_TO_DEG;
        event.core_x_m = data[98] * 0.01;
        event.core_y_m = data[118] * 0.01;
        event.array_rotation_deg = data[92] * RAD_TO_DEG;
        metadata.events.push_back(event);
        if (current_event_id == selected_shower_event_id) {
            metadata.selected_event = event;
        }
    }
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
        offsets.weight.push_back(w[static_cast<std::size_t>(i)]);
    }
    metadata.array_offsets_by_shower[current_event_id] = offsets;
    if (current_event_id == selected_shower_event_id) {
        metadata.selected_event_offsets = offsets;
    }
    return 0;
}

PhotonBunch makeBunch(const struct bunch& b,
                      int event_id,
                      int telescope_id,
                      const EventIOPhotonConfig& cfg)
{
    PhotonBunch out;
    out.photon.pos = {b.x * 0.01, b.y * 0.01, 0.0};
    out.photon.dir = {b.cx, b.cy, downwardDirZ(b.cx, b.cy)};
    out.photon.normalizeDirection();
    out.photon.time_ns = b.ctime;
    out.photon.wavelength_nm = b.lambda > 0.0 ? b.lambda : cfg.default_wavelength_nm;
    out.photon.weight = cfg.default_weight;
    out.multiplicity = b.photons * cfg.default_multiplicity;
    out.event_id = event_id;
    out.telescope_id = telescope_id;
    out.eventio_2d = true;
    return out;
}

PhotonBunch makeBunch3d(const struct bunch3d& b,
                        int event_id,
                        int telescope_id,
                        const EventIOPhotonConfig& cfg)
{
    PhotonBunch out;
    out.photon.pos = {b.x * 0.01, b.y * 0.01, b.z * 0.01};
    out.photon.dir = {b.cx, b.cy, b.cz};
    out.photon.normalizeDirection();
    out.photon.time_ns = b.ctime;
    out.photon.wavelength_nm = b.lambda > 0.0 ? b.lambda : cfg.default_wavelength_nm;
    out.photon.weight = cfg.default_weight;
    out.multiplicity = b.photons * cfg.default_multiplicity;
    out.event_id = event_id;
    out.telescope_id = telescope_id;
    out.eventio_2d = false;
    return out;
}

int readPhotonBlock(IO_BUFFER* iobuf,
                    int current_event_id,
                    const EventIOPhotonConfig& cfg,
                    const EventIOPhotonCallback& on_bunch,
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
        on_bunch(makeBunch(b, event_id, telescope_id, cfg));
        ++emitted;
        ++emitted_2d;
    }

    return get_item_end(iobuf, &item_header);
}

int readPhoton3dBlock(IO_BUFFER* iobuf,
                      int current_event_id,
                      const EventIOPhotonConfig& cfg,
                      const EventIOPhotonCallback& on_bunch,
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
    for (int i = 0; i < nbunches; ++i) {
        on_bunch(makeBunch3d(bunches[static_cast<std::size_t>(i)], event_id, telescope_id, cfg));
        ++emitted;
        ++emitted_3d;
    }
    return 0;
}

int readTelArray(IO_BUFFER* iobuf,
                 int current_event_id,
                 const EventIOPhotonConfig& cfg,
                 const EventIOPhotonCallback& on_bunch,
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
                                 stats.photon_bunches, stats.photon_bunches_2d);
        } else if (type == IO_TYPE_MC_PHOTONS3D) {
            rc = readPhoton3dBlock(iobuf, current_event_id, cfg, on_bunch,
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
                                  stats);
                break;
            case IO_TYPE_MC_PHOTONS:
                if (stop_after_current_block) {
                    rc = 0;
                    break;
                }
                rc = readPhotonBlock(iobuf.ptr, current_event_id, cfg, on_bunch,
                                     stats.photon_bunches, stats.photon_bunches_2d);
                break;
            case IO_TYPE_MC_PHOTONS3D:
                if (stop_after_current_block) {
                    rc = 0;
                    break;
                }
                rc = readPhoton3dBlock(iobuf.ptr, current_event_id, cfg, on_bunch,
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
